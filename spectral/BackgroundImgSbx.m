function [BImgMax, BImgAverage] = BackgroundImgSbx(varargin)
% Create GFP channel images from long recording
% 
% [BImgMax, BImgAverage] = BackgroundImgSbx(sbxName, spigfileName);
% 
% Optional input: 1. sbx file name (string) (with foldername if necissary)
%                 2. spsig file name to save the images to
%
% In case only one input is given (the sbx filename), the spectral file is
% assumed to be present in the same folder with the usual naming scheme.
% In case no input is given the file names will be requested.
% 
% The background images are created by taking a requested amount of chunks
% from the data, each chunk evenly spaced in time from the next chunk. 
% A chunk averages a small number of frames to reduce noise, 
% and is then max projected onto the  image that is being build up by
% the different chunks.
%
% Leander de Kraker
% 2018-10-17
%

%% You can also execute this as a script

nchunks = 300; % How many pieces of the data to take
nframes = 25; % How many frames are those pieces long

if exist('varargin', 'var') && nargin == 2 % In case both sbx and spsig file name are given
    % When called as a function, make sure the normcorr sbx file is being used
    sbxName = varargin{1};
    spsigName = varargin{2};
    % find if the filename has to be edited, mat and SPSIG has to be removed
    delimeterPos = regexp(sbxName,{'_SPSIG','.mat','sbx'});
    delimeterPos = delimeterPos(~cellfun('isempty',delimeterPos));
    if ~isempty(delimeterPos)
    	sbxName = sbxName(1:delimeterPos{1}-1);
    end
    plotter = false;
    
elseif exist('varargin', 'var') && nargin == 1 % in case only sbx filename is given
    sbxName = varargin{1};
    % find if the filename has to be edited, mat and SPSIG has to be removed
    delimeterPos = regexp(sbxName,{'_SPSIG','.mat', 'sbx'}); 
    delimeterPos = delimeterPos(~cellfun('isempty',delimeterPos));
    if ~isempty(delimeterPos)
    	sbxName = sbxName(1:delimeterPos{1}-1);
    end
    spsigName = [sbxName '_SPSIG'];
    plotter = false;
    
else  % If not called/ executing with input, request file
    [fn , pn] = uigetfile('*.sbx','get the sbx file to read the data from');
    sbxName = [pn fn(1:end-4)];
    if fn == 0 % if nothing was selected quit
        fprintf('No file selected')
        return
    end
    
    % Guess for the SPSIG file name
    spsigName = [sbxName '_SPSIG']; 
    if exist(spsigName,'file') ~= 2 % if not found, request it
        [fn , pn] = uigetfile('*_SPSIG.mat','get the SPSIG file to save the images to');
        spsigName = [pn fn(1:end-4)];
    end
    plotter = true;
    
end

sbxread(sbxName, 0,1);
global info
fprintf('Creating background image for %s\n', sbxName)

% Activate much tighter subplots
% [subplot margin top&side],[figure bottomspace,topspace],[leftspace,rightspace]
subplot = @(m,n,p) subtightplot (m, n, p, [0.075 0.005], [0.05 0.07], [0.07 0.01]);

colors = cmapL('greenFancy', 256);

nf = 250; % number of frames for example images
img = sbxread(sbxName, 0,nf);

nchannels = size(img, 3);
fprintf('%d color channel(s) seem(s) present\n', nchannels)

% Create example frames to check between beginning and end of recording
imgAv1 = mean(squeeze(img(:,:,1,:)),3); 
img = sbxread(sbxName, info.max_idx-nf-1, nf);
imgAv2 = mean(squeeze(img(:,:,1,:)),3);

% Take an average using different times in the recording.
BImgMax = zeros(size(imgAv1));
BImgAverage = zeros(size(imgAv1));
if nchannels >= 2
    BImgMaxC2 = zeros(size(imgAv1));
    BImgAverageC2 = zeros(size(imgAv1));
end

% The pieces of data are taken at equal intervals
chunkStart = floor(linspace(0, info.max_idx-nframes-1, nchunks));

if plotter
    figure('units','normalized','position',[0.25 0.2 0.5 0.4]);
    tic
end

% calculate the images
for i = 1:nchunks
    chunkImage = sbxread(sbxName, chunkStart(i), nframes);
    chunkImage2 = mean(squeeze(chunkImage(:,:,1,:)),3);
    BImgMax = max(BImgMax, chunkImage2);
    BImgAverage = BImgAverage + chunkImage2;

    if mod(i,10)==0 && plotter
        subplot(1,2,1),
        imagesc(BImgMax); colorbar
        title(sprintf('MaxProjection: iteration %d/%d. used %.3f%% of sbx file: time %.2f sec',...
            i, nchunks, (nframes*i./info.max_idx*100), toc))
        subplot(1,2,2),
        imagesc(BImgAverage/i); colorbar
        title(sprintf('average projection: maxval=%.1f',max(BImgAverage(:))))
        colormap(colors)
        pause(0.001)
        tic
    end
    
    if nchannels >= 2 % Create a background image for the 2nd color channel
        chunkImage2 = mean(squeeze(chunkImage(:,:,2,:)),3);
        BImgMaxC2 = max(BImgMaxC2, chunkImage2);
        BImgAverageC2 = BImgAverageC2 + chunkImage2;
    end
end
BImgAverage = BImgAverage/nchunks;

% Fix horizontal alignment, error from incorrect bi-directional scanning
% BImgMaxLines = BImgMax;
% BImgMaxLines(1:2:end, 1:end-1) = BImgMaxLines(1:2:end, 2:end);
% BImgAvLines = BImgAverage;
% BImgAvLines(1:2:end, 1:end-1) = BImgAvLines(1:2:end, 2:end);

figure('units','normalized','position',[0.01 0.1 0.25 0.8]);
h = gobjects(4,1);
h(1) = subplot(2,1,1);
imagesc(log1p(imgAv1))
title(sprintf('mean projection of first %d frames', nf))
h(2) = subplot(2,1,2);
imagesc(log1p(imgAv2));
title(sprintf('mean projection of last %d frames', nf))
colormap(colors)
figtitle(sbxName)

figure('units','normalized','position',[0.26 0.1 0.25 0.4]);
h(3) = subplot(1,1,1);
imagesc(log1p(BImgAverage))
colormap(colors)
title(sprintf('average %d chunks of %d frames = %d frames',nchunks, nframes, nchunks*nframes))

figure('units','normalized','position',[0.26 0.5 0.25 0.4]);
h(4) = subplot(1,1,1);
imagesc(log1p(BImgMax))
colormap(colors)
title(sprintf('max projections of %d chunks of %d averaged frames = %d frames',...
    nchunks, nframes, nchunks*nframes))
% print(gcf, [fn(1:end-18) 'img_averageFluorescence.png'],'-dpng', '-r400')

linkaxes(h, 'xy')



if nchannels >= 2
    figure('units','normalized','position',[0.53 0.1 0.3 0.8]);
    h = gobjects(2,1);
    h(1) = subplot(2,1,1);
    imagesc(log1p(BImgAverageC2))
    title(sprintf('average %d chunks of %d frames = %d frames',nchunks, nframes, nchunks*nframes))
    h(2) = subplot(2,1,2);
    imagesc(log1p(BImgMaxC2))
    title(sprintf('max projections of %d chunks of %d averaged frames = %d frames',...
        nchunks, nframes, nchunks*nframes))
    colormap([colors(:,2) colors(:,3) colors(:,1)]) % Custom shuffled black-greenish-white colormap
    figtitle(sbxName)
    linkaxes(h, 'xy')
end

infoimg.percent = nframes*i./info.max_idx*100;
infoimg.framesPerChunk = nframes;
infoimg.nchunks = nchunks;
infoimg.comment = {'BImgMax= maximum projection';'BImgAverage= average projection';...
    'BImg= low frequency spectral components';...
    'SpatialCorr= signal correlation per pixel with ROI brightest pixel'};
if nchannels >= 2
    infoimg.comment(end+1) = {'BImgMaxC2 = maximum projection of 2nd color channels'};
    infoimg.comment(end+1) = {'BImgAverageC2 = average projection of 2nd color channels'};
end
fprintf('done calculating :)\n')


%% check how the spectral looks

load(spsigName, 'Sax', 'SPic')

if exist('SPic','var') && exist('Sax', 'var')

    imgStack = log(SPic(:,:,2:end));
    Sax(1) = []; %first spectral component is the average power over al components

    vals = sort(unique(imgStack(:))); % remove -inf values
    if vals(1)==-inf
        imgStack(imgStack==-inf) = vals(2);
    end
    imgStack = permute(imgStack,[2 1 3]); % transpose the SPic variable so it's same as BImg

    nSpec = size(imgStack,3);

    % A spectral image for all the different frequencies
    colors = flip(jet(nSpec));
    norma = false; % normalize each frequency?
    figure('units','normalized','position',[0.52 0.5 0.25 0.4]);
    subplot(1,1,1)
    SPicTCell = squeeze(num2cell(imgStack, [1 2]));
    RGBspectral = CreateRGB2(SPicTCell, 1:nSpec, colors, norma);
    
    imagesc(RGBspectral)
    title('The intensity of the spectral density at all different frequencies')
    % Set the correct colorbar for this image
    colormap(colors)
    h = colorbar;
    hTicks = linspace(Sax(1),Sax(end),length(h.Ticks));
    h.Ticks = h.Ticks; % prevent ticks from changing
    h.TickLabels = num2str(hTicks', '%2.2f');
    ylabel(h, 'spectral frequency', 'FontSize',12)
else
    fprintf('Unable to load SPic and Sax (Spectral Axis) from spsigName\n')
    fprintf('Not creating colored spectral image\n')
end

%% The end: save it

save(spsigName, 'BImgMax','BImgAverage','infoimg','-append')
if nchannels >= 2
    save(spsigName, 'BImgMaxC2','BImgAverageC2','-append')
end
fprintf('saved background images to %s\n\n',spsigName)

end

