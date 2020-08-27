% Get signals from patches of sbx TRANS files
% Set how many patches should be extracted
% It's possible to load in a mask, and ROIs denoted by the mask will be
% excluded from signal retrieval
%
%
% Leander de kraker
% 2019-7-10
%


%% Load file
[fn, pn] = uigetfile('.dat');
[sbxt, dim, freq] = transmemap([pn fn]);

removeROIs = false;
% the pretty visualisation takes too long with many patches, so there is a
% faster visualisation for above maxnPlot amount of patches
maxnPlot = 150;

%% The number of rows and column that will get extracted from the sbx file
nX = 12; % across the width 
nY = 10; % across the height
n = nY*nX; % update the number of patches

fprintf('will divide the sbx file into %d x %d = %d patches\n', nX, nY, n);
colors = jet(n);

%% OR The number of signals that will get extracted from the sbx file
n = 50;
fprintf('requested number of patches: %d\n',n)

% Calculate how many rows and columns need to be made for this many patches
nX = ceil(sqrt(n));
if (nX-1) * nX >= n
    nY = nX-1;
else
    nY = nX; 
end
n = nY*nX; % update the number of patches

fprintf('will divide the sbx file into %d x %d = %d patches\n', nX, nY, n);
colors = jet(n);


%% Remove existing ROIs using a mask?

load([pn fn(1:end-9) 'SPSIG.mat'],'Mask')
Mask = Mask';
removeROIs = true;

%% Divide the sbx file into n patches

% Extracting a frame from the transposed data file (may take some time)
fprintf('extracting an example frame...\n')
tic
frameExample = sbxt.Data.y(1,:);
frameExample = double(reshape(frameExample, [dim(2), dim(3)]));
if toc>60; fprintf('finally '); end
fprintf('done with example frame after %.1fsec\n',toc)

limsx = floor(linspace(1, dim(2), nX+1));
limsy = floor(linspace(1, dim(3), nY+1));

% preallocate the signal
signal = zeros(dim(1),n);
maskPatch = zeros(dim(2),dim(3));

% Extract the signals and plot them
figure;
subplot(2,1,1)
if n<maxnPlot
    rainbowFrame = repmat({frameExample}, [n,1]); % frameExample colored according to the patches
end

PP = struct('A',[], 'P',[], 'Con',struct('x',[],'y',[]), 'Cnt',n);
i = 1;
tic
for x = 1:nX
    for y = 1:nY
        maskPatch(limsx(x):limsx(x+1)-1, limsy(y):limsy(y+1)-1) = i;
        
        if removeROIs
            maskPatch(Mask>0) = 0; % remove the existing ROIs
        end
        if n<maxnPlot
            rainbowFrame{i}(maskPatch~=i) = 0; % Set pixels that are not part of this one as 0
        end
        
        % Get which indexes correspond to the ith patch of the data
        idx = find(maskPatch(:)==i);
        signal(:,i) = mean(sbxt.Data.y(:,idx),2); % Extract the signal from these indexes

        plot(signal(1:1500,i),'color', [colors(i,:) 0.5])
        hold on
        
        % Get a contour
        PP.A(i) = length(idx);
        con = bwboundaries(maskPatch==i);
        biggest = 1;
        biggestSize = 0;
        for c = 1:length(con)
            if biggestSize<length(con{c})
                biggest = c;
                biggestSize = length(con{c});
            end
        end
        PP.P(1,i) = mean(con{biggest}(:,1));
        PP.P(2,i) = mean(con{biggest}(:,2));
        PP.Con(i).x = con{biggest}(:,1);
        PP.Con(i).y = con{biggest}(:,2);
        
        fprintf('done with patch %d/%d, elapsed=%.0f sec\n', i, n, toc)
        i = i + 1;
    end
end
xlabel('time (frame)')
ylabel('raw fluorescence')
title(sprintf('signals from the %d (%d x %d) patches', n, nX, nY))

subplot(2,1,2) % plot where the signals originate from
if n<maxnPlot
    rainbowFrameRGB = CreateRGB2(rainbowFrame, colors);
    imagesc(permute(rainbowFrameRGB, [2 1 3])*2)
    title('example frame colored same as extracted signal')
else
    imagesc(log(double(frameExample))')
    colormap(gray)
    title('example frame with colors same as extracted signal')
end
hold on
for i = 1:PP.Cnt
    plot(PP.Con(i).x, PP.Con(i).y, 'color', colors(i,:))
end


%% Create df/f

xas = (1:dim(1))./freq;
window = round(5000*freq/15); %adapted for framerate
if length(xas) < window
    window = length(xas)-1; %otherwise you will get inf
end

% baseline estimation due to Pnevmatikakis et al
sigDff = basecorrect(signal, window);


%% Create Simple Estimate of Activity by Leander (SEAL)
% SEAL takes the non-negative derivative of slightly smoothed raw data,
% which quite closely, but more nuancetly resembles deconvolved data
% .... Also it's like 100000x faster than deconvolution
noiseFloor = 0;
% filter parameters to reduce noise
windowSize = 2;

seal = imgaussfilt(sigDff, [windowSize, 0.01]);
seal = diff(seal) - noiseFloor;
seal(end+1, :) = 0; % fill at end, because derivative is 1 shorter then original
seal(seal<0) = 0;

%% Creating a fake decon signal based on seal

repressed = zeros(1,size(sigDff,2)); % how much does the signal need to be repressed?
decay = 1.2; % how much does the repression decay per timestep?
decon = zeros(size(seal)); % the fake decon signal
threshold = prctile(sigDff,80); % the threshold is different per ROI
for t = 1:size(seal,1)
    spiking = sigDff(t,:)-repressed > threshold;
    repressed = repressed / decay; % repression decays
    repressed = repressed + spiking.*0.5; % neurons that spiked need to be repressed
    decon(t,:) = spiking;
end
fprintf('done creating fake decon signal\n')


%% Create deconvolved data
% Leander did not include deconvolution yet because the signal has too many
% patches it will take too long to deconvolve

%% Plot the signals
if exist('decon','var')
    nplots = 4;
else
    nplots = 3;
end

len = 4000;
x = xas(1:len);
h = gobjects(nplots,1);
figure
h(1) = subplot(nplots,1,1);
imagesc(x, 1:n, signal(1:len,:)')
title('raw signal')
h(2) = subplot(nplots,1,2);
imagesc(x, 1:n, sigDff(1:len,:)')
title('Df/f')
h(3) = subplot(nplots,1,3);
imagesc(x, 1:n, seal(1:len,:)')
title('Simple Estimate of Activity by Leander (SEAL)')
xlabel('time (sec)')
ylabel('patch (n)')
if exist('decon','var')
    h(4) = subplot(nplots,1,4);
    imagesc(x, 1:n, decon(1:len,:)')
    title('deconvolved HAA no code fuck that')
end
linkaxes(h,'xy')

%% Save the signals als a file which can be used in extract

savename = [fn(1:end-9), sprintf('signals-%d-Patches_SPSIG',n)];

sigraw = signal;
sig = sigDff;
den = sig; % den and decon NEED to exist for extract not to crash :(
Mask = maskPatch';
if ~exist('BImg','var')
    BImg = frameExample;
end
BImg = BImg';

save(savename, 'freq', 'sigraw', 'sig', 'seal', 'den', 'decon', 'Mask', 'PP','BImg')
if n<maxnPlot
    save(savename,'rainbowFrameRGB','-append')
end

%% Execute extract

extract_nopop()

