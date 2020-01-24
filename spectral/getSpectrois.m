function getSpectrois(varargin)
% 
% get rois from cross spectral components
%
% input : SPSIG file : series of spectral images 
%         spar : parameter structure (optional)
% output: Mask: an image with the same size as BImg, with numbers that
%               denote which ROI is present at which pixel
%         
%         PP: struct with ROI information
%         PP.Cnt: The number of ROIs
%         PP.Con contours: for every ROI defined by x and y values
%         PP.A: area of the contour in pixels
%         PP.Rvar : average pixel covariance
%         PP.SpecProfile: the mean values of spectral density for every
%                         frequency in every ROI
%         PP.peakFreq: Which frequency has the highest spectral power value
%         PP.P: Values for every ROI
%         PP.P: row 1: x values of the pixel with highest value in the BImg
%         PP.P: row 2: y values of the pixel with highest value in the BImg
%         PP.P: row 3: the maximum value of BImg in each ROI
%
%         SpatialCorr:
%         spar: Spectral PARameters used to find the ROIs
%
%
global DISPLAY
DISPLAY = 0;
if exist('varargin', 'var') && nargin == 2
    filename = varargin{1};

    spar = varargin{2};
    flds = fieldnames(spar);
    aflds = { 'cutOffHz', 'border', 'areasz', 'roundedness', 'voxel', 'cutoffcorr'};
    if sum(ismember(aflds, flds)) < 6
        disp('Error; number of input values is not valid; please run : spar = Spectroiparm()')
        return;
    end

elseif exist('varargin', 'var') && nargin == 1
    filename = varargin{1};

else
    [fn, pn] = uigetfile('*_SPSIG.mat');
    filename = [pn fn];

    spar = Spectroiparm(); %reads roi segmentation parameters from file
end

%% Load and Process Spectral Images:  load('SPic.mat')
fprintf('\nloading...')
load(filename, 'SPic', 'Sax')
Transfile = [filename(1:end-9) 'Trans.dat'];
[sbxt, ~, freq] = transmemap(Transfile);
fprintf('\nloaded\n')

% obtain average spectral density for each image: decays exponentially
imgStack = log(SPic(:,:,2:end));
Sax(1) = []; %first spectral component is the average power over al components

vals = sort(unique(imgStack(:))); % remove -inf values
if vals(1)==-inf
    imgStack(imgStack==-inf) = vals(2);
end
imgStackT = permute(imgStack,[2 1 3]); % transpose the SPic variable so it's same as BImg

cutOffHz = spar.cutOffHz;
Slow = (Sax <= cutOffHz);
Spect = imgStack(:,:,Slow);
Spect = permute(Spect, [2 1 3]); %transposes the images,
Spect = setminlevel(Spect); %set minimum level to zero


% Activate much tighter subplots
% [subplot margin top&side],[figure bottomspace,topspace],[leftspace,rightspace]
%subplot = @(m,n,p) subtightplot (m, n, p, [0.075 0.005], [0.05 0.07], [0.09 0.01]);

BImg = max(Spect, [], 3);

figure('units','normalized','position',[0.51 0.1 0.25 0.4]);
subplot(1,1,1)
imagesc(BImg), caxis(prctile(BImg(:), [0.01 99.9])), hold on
title(sprintf('spectral components max projection <%.1fHz', cutOffHz))

%% ROI selection based on cross spectral density
% find ROIs with higest local maxima in low frequency component images

PP = [];
PP.Cnt = 0; 
dim = size(Spect); %width, height, z 
Mask = zeros(dim(1:2)); 
SpatialCorr = zeros(dim(1:2));

for i = 1:dim(3)
    tic
    fprintf('iteration %2d/%2d: ', i, dim(3))
    Img = Spect(:,:,i); 
    Img(Mask>0) = 0;
    figure(1), hold off,imagesc(Img), colormap gray, hold on %, caxis([h0 h1])     
    [PP, Mask, SpatialCorr] = roisfromlocalmax(Img, PP, Mask, spar, sbxt, freq, SpatialCorr);   

    str = sprintf('%.2fHz, number of ROIs found: %5d. time elapsed = %.2fminutes\n', Sax(i), PP.Cnt, toc/60);
    fprintf(str)
    
    Con = PP.Con;
    for k = 1:PP.Cnt
        plot(Con(k).x, Con(k).y, 'r')
    end
    title(str)
    pause(0.01)
end

% swap x and y values, otherwise positions are not correct in Displayrois
P = PP.P; % plot(x,y), plots at the position (y, x) in an image
PP.P(1,:) = P(2,:);
PP.P(2,:) = P(1,:);

% Retrieve the average spectral profile of the ROI
[PP.SpecProfile, PP.peakFreq] = SpecProfileCalcFun(imgStackT, Mask, Sax);

for i = 1:PP.Cnt
    %also redefine peak maxima based on projected maxima in BImg
    ROIi = Mask == i; 
    PP.P(3,i) = max(BImg(ROIi)); % maximum
end

save(filename, 'PP', 'Mask', 'BImg', 'spar', 'SpatialCorr', '-append')
fprintf('saved.\n')
