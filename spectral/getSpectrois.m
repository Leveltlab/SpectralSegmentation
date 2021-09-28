function getSpectrois(varargin)
% 
% get rois from cross spectral components
%
% GUI for setting spar (Spectral search PARameters) will be started except
% if getSpectrois is called with a spar variable that has all the correct
%   fields
%
% input : SPSIG file : series of spectral images 
%         spar : parameter structure (optional)
%
% output: Mask (H x W)(2D double): Numbers denote which ROI is present at 
%                                  which pixel
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
%         BImg (H x W) (2D double): maximum projection of the used spectral
%                                   images.
%         
%         SpatialCorr (H x W) (2D double): Signal correlation values for
%                 every ROI, from seedpoint to the other pixels in that ROI
%
%         spar: Spectral PARameters used to find the ROIs
%
%
global DISPLAY
global spar
DISPLAY = false;

if exist('varargin', 'var') && nargin >= 1
    filenameSPSIG = varargin{1};
else
    [fn, pn] = uigetfile('*_SPSIG.mat');
    filenameSPSIG = [pn fn];
end

if exist('varargin', 'var') && nargin >= 2
    spar = varargin{2};
    runSparArm = false;
    flds = fieldnames(spar);
    aflds = {'cutOffHzMax', 'cutOffHzMin', 'border', 'areasz',...
             'roundedness', 'voxel', 'cutOffCorr','useFluorescenceImg'};
    if sum(ismember(aflds, flds)) < 8
        warning('number of input values is not valid')
        runSparArm = true;
    end
else
    runSparArm = true;
end

%% Load and Process Spectral Images:  load('SPic.mat')
fprintf('\nloading...')
load(filenameSPSIG, 'SPic', 'Sax')

sfn = regexp(filenameSPSIG,'SPSIG', 'split');
filenameTrans = [sfn{1} 'DecTrans.dat'];

if ~isfile(filenameTrans) %Using decimated data
    filenameTrans = [sfn{1} 'Trans.dat'];
    if ~isfile(filenameTrans) % Resort to manual selection of transposed data
        [fnTrans, pnTrans] = uigetfile('*.dat', 'Where is the transposed data file (prefering decimated data)?');
        if fnTrans == 0
            fprintf('No transposed data found, unable to do automatic ROI creation\n')
            return
        else
            filenameTrans = [pnTrans fnTrans];
        end
    end
end
[sbxt, ~, freq] = transmemap(filenameTrans);
fprintf('Memory mapped %s\n', filenameTrans)

% obtain average spectral density for each image: decays exponentially
imgStack = log(SPic(:,:,2:end));
Sax(1) = []; %first spectral component is the average power over al components

imgStackT = permute(imgStack,[2 1 3]); % transpose the SPic variable so it's same as BImg
imgStackT = setminlevel(imgStackT); %replaces -infs and subtracts minimum

if runSparArm
%      spar = Spectroiparm(); %reads roi segmentation parameters from file
	h = SpectParArm(imgStackT, Sax); % arm the spectral parameters (spar variable)
    waitfor(h)
end

%%
selectedFreq = (Sax >= spar.cutOffHzMin) & (Sax <= spar.cutOffHzMax);
Spect = imgStackT(:,:,selectedFreq);
SaxUsed = Sax;
if spar.useFluorescenceImg % Add fluorescence images to the spectral images
    load(filenameSPSIG, 'BImgMax', 'BImgAverage')
    if exist('BImgMax', 'var')
        maxSpect = max(Spect(:));
        if max(BImgAverage(:)) > maxSpect
            BImgAverage = BImgAverage ./ max(BImgAverage(:)) .* maxSpect;
        end
        if max(BImgMax(:)) > maxSpect
            BImgMax = BImgMax ./ max(BImgMax(:)) .* maxSpect;
        end
        Spect(:,:,end+1) = BImgMax;
        Spect(:,:,end+1) = BImgAverage;
        SaxUsed(end+1:end+2) = 0;
    else
        warning('Tried to use BImgMax and BImgAverage in ROI search but they are absent')
    end
end

BImg = max(Spect, [], 3);

figure('units','normalized','position',[0.51 0.1 0.25 0.4]);
imagesc(BImg), caxis(prctile(BImg(:), [0.01 99.9])), hold on
title(sprintf('spectral components max projection >%.2fHz & <%.2fHz',...
                spar.cutOffHzMin, spar.cutOffHzMax))

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
    figure(1), hold off, imagesc(Img), colormap gray, hold on   
    [PP, Mask, SpatialCorr, rlog] = roisfromlocalmax(Img, PP, Mask, spar, sbxt, freq, SpatialCorr);   

    str = sprintf('number of ROIs found (%.2fHz): %5d. time elapsed = %.2fminutes\n', SaxUsed(i), PP.Cnt, toc/60);
    fprintf(str)
    
    SummaryGetRois(rlog, spar)
    
    if isfield(PP, 'Con')
        figure(1)
        Con = PP.Con;
        for k = 1:PP.Cnt
            plot(Con(k).x, Con(k).y, 'r')
        end
        title(str,'FontSize',10, 'FontWeight', 'normal')
        pause(0.05)
    end
end

% swap x and y values, otherwise positions are not correct in Displayrois
P = PP.P; % plot(x,y), plots at the position (y, x) in an image
PP.P(1,:) = P(2,:);
PP.P(2,:) = P(1,:);
PP.P(3,:) = []; % remove maximum as it is specific for the frequency in which it was found

% Retrieve the average spectral profile of the ROI
[PP.SpecProfile, PP.peakFreq, PP.peakVal] = SpecProfileCalcFun(imgStackT, Mask, 1:PP.Cnt, Sax);

% Save
specUsed = Spect;
save(filenameSPSIG, 'PP', 'Mask', 'BImg', 'spar', 'SpatialCorr', 'specUsed', '-append')
fprintf('saved.\n')



