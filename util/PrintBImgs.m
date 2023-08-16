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
if ~exist('varargin', 'var')
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


%% Image settings
baseName =  [filepathImg, sprintf('img_scalebar%dum', um)];

ncolors = 256;
colors = cmapL('greenFancy', ncolors);

pos = GetCenteredFigPos(dims([2 1])); % figure that does not get saved


%% BImgMax (Maximum fluorescence projection) 

rgb = BImgMax;
rgb = MidtoneBalance(rgb, 0.3);


brightnessCutOff = 99.99; % Values above this percentile will be max color value. 
                          % 2p data brightness can be highly skewed 
                          % this causes most pixels to be dark. This is not
                          % nice to look at. So select lower than 100% to
                          % cut off the brightness and thereby making the
                          % image brighter/ more contrast
darknessCutOff = 0.1; % percentile of values below which the values will be the darkest color.
lims = prctile(rgb, [darknessCutOff, brightnessCutOff]); % color axis limits

[rgb, s,  m]      = ShiftLinesImg(rgb);
[rgb, s(2), m(2)] = ShiftLinesImg(rgb);
rgb = SetLimits(rgb, lims, colors);

% Apply scalebar
rgb(scaleBary, scaleBarx, :) = 1;

% For our eyes only
if doPlotForOurEyes
    figure('Units', 'Pixels', 'InnerPosition', pos)
    subplot('Position', [0 0 1 1])
    image(rgb);
end

% save image
imwrite(rgb, [baseName, '_maxFluorescence2_', filename, '.png'])


%% BImgAverage (average fluorescence projection)

rgb = BImgAverage;

brightnessCutOff = 99.999;
darknessCutOff = 0.1; % percentile of values below which the values will be the darkest color.
lims = prctile(rgb(:), [darknessCutOff, brightnessCutOff]); % color axis limits

[rgb, s,  m]      = ShiftLinesImg(rgb);
[rgb, s(2), m(2)] = ShiftLinesImg(rgb);
rgb = SetLimits(rgb, lims, colors);
rgb(scaleBary, scaleBarx, :) = 1;


% For our eyes only
if doPlotForOurEyes
    figure('Units', 'Pixels', 'InnerPosition', pos)
    subplot('Position', [0 0 1 1])
    image(rgb);
end

% save image
imwrite(rgb, [baseName, '_averageFluorescence_', filename, '.png'])

%% Spectral image

rgb = BImg;

brightnessCutOff = 99.99;
darknessCutOff = 0.1; % percentile of values below which the values will be the darkest color.
if sum(rgb(:)==-inf)>0
    vals = unique(rgb(:));
    rgb(rgb==-inf) = vals(2);
end
lims = prctile(rgb(:), [darknessCutOff, brightnessCutOff]); % color axis limits

colorsGray = gray(256);
rgb = SetLimits(rgb, lims, colorsGray);
rgb(scaleBary, scaleBarx, :) = 1;

% For our eyes only
if doPlotForOurEyes
    figure('Units', 'Pixels', 'InnerPosition', pos)
    subplot('Position', [0 0 1 1])
    image(rgb);
end

% save image
imwrite(rgb, [baseName, '_spectral_', filename, '.png'])


%% colored spectral image

% % imgRGB = BImgRGB;
% freqToUse = [0, 0.25]; % frequency range in Hz
% 
% figure('Units', 'Pixels', 'InnerPosition', pos)
% subplot('Position', [0 0 1 1])
% [imgRGB, colorsSpectral, colorbarVals] = SpectralColorImg('data', {SPic, Sax}, freqToUse, true);
% imgRGB = MidtoneBalance(imgRGB, 0.4);
% imgRGB(scaleBary, scaleBarx, :) = 1;
% h = gca;
% h.Children.CData = imgRGB;
% 
% % save image
% imwrite(imgRGB, [baseName, '_coloredSpectral', filename, '.png'])

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










