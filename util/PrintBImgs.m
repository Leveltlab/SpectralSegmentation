% Quickly save the average & max fluorescence & colored spectral as png
% with their actual resolution
% BImgAverage, BImgMax and BImgRGB are created by FluorescenceImgSbx.m
% 
% 
% Leander de Kraker
% 2022-10-20
%

% Load data
[filename, filepath] = uigetfile('*SPSIG.mat');
load([filepath, filename], 'BImgMax', 'BImgAverage', 'BImgRGB', 'SPic', 'Sax')

delimeterPos = regexp(filename,{'_SPSIG','.mat','sbx'});
delimeterPos = delimeterPos(~cellfun('isempty',delimeterPos));
if ~isempty(delimeterPos)
    filename = filename(1:delimeterPos{1}-1);
end


%% Scale bar

dims = size(BImgMax);

% Set these parameters manually!
fov = 1000; % amount of um the microscope images with 1x zoom
zoom = 1.6; % MICROSCOPE ZOOM LEVEL, FROM LOGBOOK
scale =  dims(2) / (fov/zoom); % number of pixels / µm images (depends on zoom)
um = 100; % How big you want your scalebar

scaleBar = round(scale*um);

scaleBary = round(dims(1)*0.9); % Set vertical position of scalebar
scaleBary = scaleBary:scaleBary+2; % Set thickness of scalebar
scaleBarx = round(dims(2)*0.9); % Setting the x position of scalebar
scaleBarx = scaleBarx-scaleBar+1:scaleBarx;


%% Image settings
baseName =  [filepath, filename, sprintf('_img_scalebar%dum', um)];

ncolors = 256;
colors = cmapL('greenFancy', ncolors);

pos = GetCenteredFigPos(dims([2 1])); % figure that does not get saved


%% BImgMax (Maximum fluorescence projection) 
brightnessCutOff = 99.99; % Values above this percentile will be max color value. 
                          % 2p data brightness can be highly skewed 
                          % this causes most pixels to be dark. This is not
                          % nice to look at. So select lower than 100% to
                          % cut off the brightness and thereby making the
                          % image brighter/ more contrast
darknessCutOff = 0.1; % percentile of values below which the values will be the darkest color.
lims = prctile(BImgMax(:), [darknessCutOff, brightnessCutOff]); % color axis limits


rgb = SetLimits(BImgMax, lims, colors);

% Apply scalebar
rgb(scaleBary, scaleBarx, :) = 1;

% For our eyes only
figure('Units', 'Pixels', 'InnerPosition', pos)
subplot('Position', [0 0 1 1])
image(rgb);

% save image
imwrite(rgb, [baseName, '_maxFluorescence.png'])


%% BImgAverage (average fluorescence projection)

brightnessCutOff = 99.999;
darknessCutOff = 0.1; % percentile of values below which the values will be the darkest color.
lims = prctile(BImgAverage(:), [darknessCutOff, brightnessCutOff]); % color axis limits

rgb = SetLimits(BImgAverage, lims, colors);
rgb(scaleBary, scaleBarx, :) = 1;


% For our eyes only
figure('Units', 'Pixels', 'InnerPosition', pos)
subplot('Position', [0 0 1 1])
image(rgb);

% save image
imwrite(rgb, [baseName, '_averageFluorescence.png'])


%% colored spectral image

imgRGB = BImgRGB;
freqToUse = [0, 0.25]; % frequency range in Hz

figure('Units', 'Pixels', 'InnerPosition', pos)
subplot('Position', [0 0 1 1])
[imgRGB, colorsSpectral, colorbarVals] = SpectralColorImg('data', {SPic, Sax}, freqToUse, true);
imgRGB = MidtoneBalance(imgRGB, 0.4);
imgRGB(scaleBary, scaleBarx, :) = 1;
h = gca;
h.Children.CData = imgRGB;

% save image
imwrite(imgRGB, [baseName, '_coloredSpectral.png'] )

%% Color spectral image colorbar

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








