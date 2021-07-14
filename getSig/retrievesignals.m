function infoSigPar = retrievesignals(filename)
% retrieve timetrace signals from a transposed file _Trans.dat using ROIs
% 
% input: filename = (string) spectral filename with '.mat', possibly with path
% In case no input is given the script will ask for the filename
%
% Output: infoSigPar = info on Signal retrieval Parameters (struct)
%
% Neuropil removal is attempted by simple substraction of signal of
% surrounding pixels (the background ROIs): 
% corrected sig = raw - 0.7 * background sig
%
% version changes:
%   version 1.3 (2019-3-11): 'a' set to 0.7
% 	version 1.2 (2019-2-26): Changes in plotting and SEAL variables corrected
%	version 1.1            : Changes in plotting
%
% Chris van der Togt
% edited by Leander
% 2019-2-6

% Variables of the neuropil properties
alpha = 0.7; % for 'y = a*x + b'. 0.7 is used in other research groups
bufSize = 4; % n pixels buffer around ROIs before background ROI starts
surround = 11; % take area of n pixels around the buffer as background ROI

infoSigPar = struct(); % info with script parameters
infoSigPar.scriptVersion  = 1.3;
infoSigPar.alpha  = alpha; % alpha of neuropil substraction
infoSigPar.bufSize = bufSize;
infoSigPar.surround = surround;

% Activate much tighter subplots
% [subplot margin top&side],[figure bottomspace,topspace],[leftspace,rightspace]
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.06], [0.1 0.1], [0.1 0.02]);

% Select file if code is being called as script
try
    asfunction = nargin == 1;
catch
    asfunction = false;
end
if ~asfunction
    [fn, pn] = uigetfile('*_SPSIG.mat');
    if fn == 0 % if cancelled
        fprintf('no file selected: exiting\n')
        return
    end
    filename = [pn, fn];
end

% Load the data
load(filename, 'Mask', 'PP', 'BImg');

% memorymap transposed data file
filenameTrans = [filename(1:end-9) '_Trans.dat'];
if ~isfile(filenameTrans) % Resort to manual selection of transposed data
    [fnTrans, pnTrans] = uigetfile('*.dat', 'Where is the transposed data file?');
    if fnTrans == 0
        fprintf('No transposed data found, unable to do signal retrieval\n')
        return
    else
        filenameTrans = [pnTrans fnTrans];
    end
end
[sbxt, dim, freq] = transmemap(filenameTrans); 
fprintf('\nloaded file %s\n',filenameTrans(1:end-4))

% Create background ROIs
[~, ~, buf] = BufferMask(Mask, bufSize);
% Background contours and seperate masks
[backCon, zones, backMask] = BufferMask(buf+Mask, surround); 

% Plotting the contours
figure('Units','normalized', 'Position', [0 0.2 0.475 0.5], 'Name', 'signal retrieval');
subplot(1,2,1); 
imagesc(BImg)
colormap gray, hold on
cmap = winter(length(backCon));
for i = 1:PP.Cnt
    for j = 1:length(backCon(i).C)
        plot(backCon(i).C(j).x, backCon(i).C(j).y,'linewidth',1.5,'color',cmap(i,:))
    end
    plot(PP.Con(i).x, PP.Con(i).y, 'r')
end
title(sprintf('%s',filename(1:end-4)))
drawnow

Mt = Mask'; %the mask needs to be transposed
idx = find(Mt(:)> 0); %get the image indices for all rois

% Extracting soma signals: enough memory for all ROIs
fprintf('extracting neuron ROI signals...\n')
Sigrois = sbxt.Data.y(:,idx);   
Sigrois = double(Sigrois);
sigraw = zeros(dim(1),PP.Cnt);
tic
for i = 1:PP.Cnt
    ichan = find(Mt(:)== i);
    [~, ia] = intersect(idx, ichan); %get indices corresponding to this roi
    sigraw(:,i) = mean(Sigrois(:,ia),2); %and average this selection 
end
fprintf('done: %.0f sec\n', toc)

% Extracting neuropil signals
Mt = backMask'; %the mask needs to be transposed
sigrawBack = zeros(dim(1), PP.Cnt);
fprintf('extracting neuropil ROI signals.....\n')
tic
for i = 1:PP.Cnt
    % The memory is probably not large enough to get all the ROIs in one go
    % thats why the memorymapped data is being extracted in the for loop
    Mt(zones(i).y, zones(i).x) = zones(i).mask'; % add transposed background mask i
    idx = find(Mt(:)==i); % get the image indices for this roi
    % indices = indices + double(dim(2)); %add a line?! OR NO LINE ?!?!!
    sigrawBack(:,i) = mean(double(sbxt.Data.y(:,idx)),2);
    
    if mod(i, 50)==1
        fprintf('ROI %d/ %d done\n', i, PP.Cnt)
    end
end
fprintf('done: %.0f sec\n', toc)


% The neuropil signals are lowered to around 0 so the corrected raw signal
% is not much lower then the original raw signal
sigrawCorrected = sigraw - alpha*(sigrawBack - repmat(prctile(sigrawBack,2), [size(sigrawBack,1), 1]));

xas = (1:length(sigraw))./freq;

% Baseline correction
window = round(5000*freq/15); %adapted for framerate
if length(xas) < window
    window = length(xas)-1; %otherwise you will get inf
end

% baseline estimation due to Pnevmatikakis et al
sig = basecorrect(sigraw, window);
sigBack = basecorrect(sigrawBack,window);
sigCorrected = basecorrect(sigrawCorrected, window);


% Save the retrieved signal
save(filename, 'sig', 'sigraw', 'sigBack', 'sigCorrected', ...
               'freq', 'infoSigPar', '-append')
fprintf('saved fluorescence and dF\\F0\n')


% Plot retrieved df/f
subplot(1,2,2);
RGB = CreateRGB({log1p(sig)', sigrawBack'}, 'rg b');
imagesc(xas./60, 1:size(sig,2), RGB); 
title('df/f signal (yellow=soma, blue=background'); xlabel('time (minutes)')
xlim([0 xas(end)/60]); ylabel('ROI (n)')

