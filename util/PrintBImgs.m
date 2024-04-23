function [s, m] = PrintBImgs(filepath, filename, filepathImg, zoom, varargin)
% Quickly save the average & max fluorescence & colored spectral as png
% with their actual resolution
% BImgAverage, BImgMax and BImgRGB are created by FluorescenceImgSbx.m
% 
% 
% Leander de Kraker
% 2022-10-20
%
%%
doPlotForOurEyes = true;
if ~exist('varargin', 'var') || (exist('varargin', 'var') && nargin == 0)
    zoom = 1.3; % MICROSCOPE ZOOM LEVEL, FROM LOGBOOK
    [filename, filepath] = uigetfile('*SPSIG.mat');
    filepathImg = filepath;
    doScalebarPlot = true;
else
    doScalebarPlot = varargin{1};
    if length(varargin)==2
        doPlotForOurEyes = varargin{2};
    end
end

% Load data
load([filepath, filename], 'BImgMax', 'BImgAverage', 'BImg', 'BImgRGB', 'SPic', 'Sax')


delimeterPos = regexp(filename,{'_SPSIG','.mat','sbx'});
delimeterPos = delimeterPos(~cellfun('isempty',delimeterPos));
if ~isempty(delimeterPos)
    filename = filename(1:delimeterPos{1}-1);
end

lineShiftCorrection = 0;

%% Scale bar

dims = size(BImgMax);

% Set these parameters manually!
fov = 1000; % amount of um the microscope images with 1x zoom
scale =  dims(2) / (fov/zoom); % number of pixels / µm images (depends on zoom)
um = 100; % How big you want your scalebar

scaleBar = round(scale*um);

scaleBary = round(dims(1)*0.9); % Set vertical position of scalebar
scaleBary = scaleBary:scaleBary+2; % Set thickness of scalebar
scaleBarx = round(dims(2)*0.9); % Setting the x position of scalebar
scaleBarx = scaleBarx-scaleBar+1:scaleBarx;
scaleBar = {scaleBarx, scaleBary};


%% Image settings
baseName =  [filepathImg, sprintf('img_scalebar%dum', um)];

ncolors = 256;
colors = cmapL('greenFancy', ncolors);
doPerc = true; % limits are given as percentiles of image values (not actual amount)
doLineShift = false;
pos = GetCenteredFigPos(dims([2 1])); % figure that does not get saved


%% BImgMax (Maximum fluorescence projection) 
limsPerc = [0.1, 99.99];
saveName = [baseName, '_maxFluorescence_', filename, '.png'];
img = MidtoneBalance(BImgMax, 0.3);
[rgb, s, m] = PrintBImg(img, limsPerc, colors,  doPlotForOurEyes, saveName, scaleBar, doPerc, doLineShift, pos);


%% BImgAverage (average fluorescence projection)
limsPerc = [0.1, 99.99];
saveName = [baseName, '_averageFluorescence_', filename, '.png'];
rgb = PrintBImg(BImgAverage, limsPerc, colors,  doPlotForOurEyes, saveName, scaleBar, doPerc, doLineShift, pos);


%% Spectral image
limsPerc = [0.1, 99.99];
saveName = [baseName, '_spectral_', filename, '.png'];
rgb = BImg;
if sum(rgb(:)==-inf)>0
    vals = unique(rgb(:));
    rgb(rgb==-inf) = vals(2);
end
rgb = PrintBImg(rgb, limsPerc, gray(256),  doPlotForOurEyes, saveName, scaleBar, doPerc, doLineShift, pos);


%% colored spectral image

% imgRGB = BImgRGB;
freqToUse = [0, 0.25]; % frequency range in Hz

figure('Units', 'Pixels', 'InnerPosition', pos)
subplot('Position', [0 0 1 1])
[imgRGB, colorsSpectral, colorbarVals] = SpectralColorImg('data', {SPic, Sax}, freqToUse, true);
imgRGB = MidtoneBalance(imgRGB, 0.4);
imgRGB(scaleBary, scaleBarx, :) = 1;
h = gca;
h.Children.CData = imgRGB;

% save image
imwrite(imgRGB, [baseName, '_coloredSpectral', filename, '.png'])

%% Color spectral image colorbar

if doScalebarPlot
    im = Colorbar2D(colorsSpectral, MidtoneBalance(gray(100), 0.6));
    
    figure('Position', [680 638 120 340])
    imagesc(1:size(im,2), colorbarVals, im)
    h = gca;
    h.YDir = 'normal';
    h.YTick = colorbarVals([1 end]); % prevent ticks from changing
    h.YTickLabel = num2str(colorbarVals([1 end]), '%.2f');
    h.XTick = [];
    xlabel('spectral power');
    ylabel('spectral frequency')
    h.FontSize = 11;
    h.XLabel.Position = [0, 0.0129, 1];
    
    SaveImg({'svg', 'png'},  [baseName, '_colorSpectralLegend'])
end










