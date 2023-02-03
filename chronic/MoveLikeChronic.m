function imgs = MoveLikeChronic(varargin)
% Opens images that are in individual SPSIG files using a chronic file 
% and registers them like the chronic matching moved each recording
% 
% Leander de Kraker
% 2022-7-11
% 

%% Load chronic file
if exist('varargin', 'var') && nargin==2
    chronicPath = varargin{1};
    chronicName = varargin{2};
else
    [chronicName, chronicPath] = uigetfile('*chronic.mat','get Chronic file');
end

load([chronicPath, chronicName], 'transformed',...
            'filenames','filepaths','filenamesShort', 'nfiles')

%% Load and chronically register requested images

% Which variables to load from each file?
toLoad = {'BImg', 'BImgAverage', 'BImgMax', 'BImgR', 'BImgG'};
% toLoad = {'BImg', 'BImgAverage', 'BImgMax'};
useimgsMean = true; % Use imgsMean from files which contain that variable?

nToLoad = length(toLoad);
imgs = struct();

for i = 1:nToLoad
    % Load the requested images
    imgsi = cell(nfiles, 1);
    for f = 1:nfiles
        notdone = true;
        if useimgsMean && ~isempty(who('-file', [filepaths{f}, filenames{f}], 'imgsMean'))
            % Use imgsMean if present because it has data from multiple files
            if i==1; fprintf('Trying to use imgsMean struct for file %d\n', f);end
            imgsi{f} = load([filepaths{f} filenames{f}], 'imgsMean');
            if isfield(imgsi{f}.imgsMean, toLoad{i})
                imgsi{f} = imgsi{f}.imgsMean.(toLoad{i});
                notdone = false; % loading is done
            end
        end
        if notdone
            imgsi{f} = load([filepaths{f} filenames{f}], toLoad{i});
            imgsi{f} = imgsi{f}.(toLoad{i});
        end
        imgsi{f} = imgsi{f}-min(imgsi{f}(:));
    end
    % Register them like the chronic matchign did
    imgsi = MoveLikeChronicFun(imgsi, transformed);
    imgs.(toLoad{i}) = imgsi;
end

%% Save the registered images to the chronic file

save([chronicPath, chronicName], 'imgs', '-append')
fprintf('\nsaved imgs to %s\n\n', chronicName)

clearvars imgsi i f notdone

%% Plot them overlayed on each other


colorsRGB = jet(nfiles);
norma = false; % normalize each image in colored overlay?
whiteBalance = true; % white balance resulting image? (makes colorbar untrue)

% Activate tighter subplots
% [subplot margin top&side],[figure bottomspace,topspace],[leftspace,rightspace]
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.1], [0.08 0.2]);

hax = gobjects(nToLoad, 1);
for i = 1:nToLoad
    imgsi  = imgs.(toLoad{i});
    
    figure;
    hax(i) = subplot(1,1,1);
    imagesc(CreateRGB2(imgsi, colorsRGB, norma, whiteBalance))
    colormap(colorsRGB)
    hcbar = colorbar;
    hcbar.YTick = ((1:nfiles)-0.5)/nfiles;
    hcbar.YTickLabel = filenamesShort;
    figtitle(sprintf('%s overlayed',toLoad{i}))
    
end
linkaxes(hax, 'xy')

%% Plot them next to each other


colors = cmapL('inferno', 256);

% Activate much tighter subplots
% [subplot margin top&side],[figure bottomspace,topspace],[leftspace,rightspace]
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.1], [0.08 0.06]);
% zoomtoD1 = [350 500]; % Vertical zoom for image
% zoomtoD2 = [300 600]; % Horizontal  zoom
zoomtoD1 = [200 500]; % Vertical zoom for image
zoomtoD2 = [200 500]; % Horizontal  zoom

% Set figure around center of screen
pos = GetCenteredFigPos([1200, 700]);

% subplots to make
nsubV = round(sqrt(nfiles));
nsubH = ceil(sqrt(nfiles));

for i = 1:nToLoad
    imgsi  = imgs.(toLoad{i});
    
    figure('Position', pos)
    h = gobjects(nfiles, 1);
    for f = 1:nfiles
        h(f) = subplot(nsubV, nsubH, f);
        imagesc(imgsi{f})
        title(filenamesShort{f})
        xlim(zoomtoD2)
        ylim(zoomtoD1)
        if f<(nfiles-nsubH+1)
            h(f).XTickLabel = [];
        end
        if mod(f,nsubH)~=1
            h(f).YTickLabel = [];
        end
    end
    colormap(colors)
    figtitle(toLoad{i})
    linkaxes(h, 'xy')
end


