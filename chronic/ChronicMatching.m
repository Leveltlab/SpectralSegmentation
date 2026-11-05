% Register spectral images of many days, and link their ROIs.
% Registration uses 2D cross-correlation via fft, with a chosen reference
% 
% Run the script step by step
% 
% Leander de Kraker
% 
% 2018-2-2
% History:
% 2020-2-24 
% 2020-5-22 Memory and speed improvement in registration.
% 2021-5-17 Big update to chronic viewer: Editing matches is now practical
% 2021-11-5 Cleaned up code. Shortening filenames & cutting edges in functions
% 2022-3-1  Rotation uses nearest neighbour interpolation to prevent noise
%               reduction which causes higher correlation with more rotation
% 2022-4-12 Improved figure spam.
% 2022-9-24 Changed some plotting, bugfixed contour registration
% 2025-8-27 Seperated the linking and saving to integrate with ImageReg2P
% 
% 

clc
%% Select SPSIG files
% load multiple files, ordered on recording time. Press cancel when done

filenames = {};
filepaths = {};
filedates = NaT(1);

selecting = true; % true as long as files are being selected
i = 0;
while selecting
    i = i + 1;
    if i == 1
        str = sprintf('load spectrally analyzed file %d', i);
    else
        str = sprintf('load spectrally analyzed file %d. previous: %s. Press cancel when done'...
            , i, datestr(filedates(i-1)));
    end
    [filenames{i}, filepaths{i}] = uigetfile('*SPSIG.mat', str);
    
    if filenames{i} == 0 % Cancel is pressed probably: stop with selecting
        filenames(i) = [];
        filepaths(i) = [];
        selecting = false;
    else % Get the date
        filedates(i) = FindDate(filepaths{i}, filenames{i});
    end
end
filenames = filenames';
filepaths = filepaths';
filedates = filedates';
nfiles = length(filenames); % The number of files that have been selected

%% Load the data of the selected files
tic

filenamesShort = ShortenFileNames(filenames, 1, filedates);

clearvars data
fprintf('\nloading in %d files:\n', nfiles)
for i = 1:nfiles
    fprintf('\nloading file %d...',i)
    data(i) = load([filepaths{i} filenames{i}], 'Mask', 'PP', 'BImg');
    fprintf('\nloaded "%s"\n', filenames{i})
end
toc

clearvars selecting str pos i questionStr answerStr strdate

%% Backing up data variable WHEN GOING TO RESIZE SOME RECORDINGS
% This of course only resizes the images in workspace. actual recordings will not be edited
data2 = data;
%% Resize some recordings? ONLY IF RESIZE IS NECESSARY! Manual parameters!
% This step is only necessary when recordings are made with different zoom

data = data2; % Take the back-up data (so you can retry easily from original sizes)

rec = [1 2 3 4 5]; % Which recordings to resize

xc = 313; % x Position of correctly registered cells
x1 = 351; % x position of neuron that is shifted in recording that will be resized 
x2 = 348;%  x position same cell as other in other recording
yc = 276; y1 = 422; y2 = 407;
scalefac    = 1 + (x1-x2)/(xc-x2); % x stretch
scalefac(2) = 1 + (y1-y2)/(yc-y2); % y stretch

% scalefac = [1.0725 1.1125];
fprintf('\n\nRECORDINGS ARE GETTING RESIZED!\n\n')
for i = 1:length(rec)
    outSize = round(scalefac .* [size(data(rec(i)).BImg,2), size(data(rec(i)).BImg,1)]);
    data(rec(i)).BImg = imresize(data(rec(i)).BImg, 'OutputSize', outSize([2 1]),'method','bilinear');
    data(rec(i)).Mask = imresize(data(rec(i)).Mask, 'OutputSize', outSize([2 1]),'method','nearest');
    for j = 1:data(rec(i)).PP.Cnt
        data(rec(i)).PP.Con(j).x = data(rec(i)).PP.Con(j).x .* scalefac(1);
        data(rec(i)).PP.Con(j).y = data(rec(i)).PP.Con(j).y .* scalefac(2);
    end
    data(rec(i)).PP.P(1,:) = data(rec(i)).PP.P(1,:) .* scalefac(1);
    data(rec(i)).PP.P(2,:) = data(rec(i)).PP.P(2,:) .* scalefac(2);
    fprintf('recording %d: %s resized!\n', rec(i), filenamesShort{rec(i)})
end


%% Get important variables from data and prepare for registration

% Preallocate variables
BImgs = cell(nfiles, 1); % plural BImg = BImgs
Masks= cell(nfiles, 1);  % plural Mask = Masks
mini = zeros(1,nfiles);
maxi = zeros(1,nfiles);

for i = 1:nfiles
    if isfield(data, 'BImg')
        img = data(i).BImg;
    elseif isfield(data, 'BImgMax')
        img = data(i).BImgMax;
    else
        img = data(i).BImgAverage;
    end
    
    % Normalize and remove -infs
    vals = unique(img);
    mini(i) = vals(round(end/20)); % 5th percentile
    maxi(i) = vals(end);
    img(img < mini(i)) = mini(i);
    img = (img - mini(i)) ./ range([mini(i), maxi(i)]);
    
    % Save mask and background image
    BImgs{i} = img;
    Masks{i} = data(i).Mask;
end
try
    PPs = [data.PP];
catch
    % probably crashing because of dissimilar structs because specprofile field
    % So let's leave those out
    PPs = struct('P', [],  'A', [], 'Cnt', [], 'Con', []);
    for i = 1:nfiles
        PPs(i).P = data(i).PP.P;
        PPs(i).A = data(i).PP.A;
        PPs(i).Cnt = data(i).PP.Cnt;
        PPs(i).Con = data(i).PP.Con;
    end
end

toplot = 1:nfiles; % sessions to plot

% Create Rainbowy color channel allocations for nfiles, with bias to blue,
colors = flipud(cmapL([0 0 1; 0 1 1; 0 1 0; 1 0.7 0; 1 0 0; 0.7 0 1], nfiles));
if size(colors,1)==3
    colors = [1 0 0; 0 1 0; 0.3 0.3 1];
elseif size(colors,1)==2
    colors = [1 0 0; 0 1 0.5];
end

% Plot image of unregistered overlayed recordings
figure('Units','Normalized', 'Position', [0.1766 0.1407 0.4203 0.7074]) 
subplot('Position', [0 0 1 1])
rgb = CreateRGB2(BImgs(toplot), colors(toplot,:));
imagesc(rgb)
for i = 1:length(toplot)
    text(20, 20+i*17, datestr(filedates(toplot(i))),'color',colors(toplot(i),:),...
        'fontweight','bold','fontsize',12)
end
hold on
for r = 1:nfiles
    for i = 1:PPs(r).Cnt
        plot(PPs(r).Con(i).x, PPs(r).Con(i).y, 'color', [colors(r,:) 0.2])
    end
    plot(PPs(r).P(1, :), PPs(r).P(2, :), '.', 'color', colors(r,:))
end
figtitle('non-registered images', 'Color', 'w')


% Prepare for registration
transformed = struct();
transformed.vert = [];
transformed.hori = [];
transformed.rotation = [];
transformed.cutting = [];
transformed.dims = [];
iter = 0;
corrScoreMean = [];

% Calculate similarity between the images
corrScore = BImgOverlapScore(BImgs);
corrScoreMean(:,end+1) = [sum(corrScore)./(nfiles-1)]'; %#ok<NBRAK>

clearvars vals mini maxi i img RGB

%% Register images with manual GUI.

%% Register images automatically (optional: run multiple times)
referenceNum = 2; % which background image to take as reference
lockRecs =  []; % Lock certain recordings (to the reference), do not move them

buf = 40; % Buffer to remove around the to-register image

moveRecs = 1:nfiles;
moveRecs(lockRecs) = [];


if ~exist('iter', 'var')
    warning('Had to recreate iter and corrScore variable')
    iter = 0;
    corrScore = BImgOverlapScore(BImgs);
    corrScoreMean(:,1) = [sum(corrScore)./(nfiles-1)]'; %#ok<NBRAK>
end
if ~exist('toplot','var'); toplot=1:nfiles;end
if ~exist('colors','var'); colors=flipud(cmapL([0 0 1; 0 1 1; 0 1 0; 1 0.7 0; 1 0 0; 0.7 0 1], nfiles));end

iter = iter + 1;

if referenceNum>nfiles
    error('Select a lower reference recording number! referenceNum exceeds number of files')
end

tic
fprintf('\nstarting registration. with recording %d as reference\n', referenceNum)
% Determine which minimum size image since recordings are different sizes
dims = cellfun(@size, BImgs,'UniformOutput',false);
dims = cat(1,dims{:});
transformed(iter).dims = dims;

h = min(dims(:,1))-buf;
w = min(dims(:,2))-buf;

moveVert = zeros(1, nfiles);
moveHori = zeros(1, nfiles);

% The base session to register others with
BImgRef = BImgs{referenceNum};
BImgRef = BImgRef - mean(BImgRef(:));

% Do cross correlations with the base image
for i = 1:nfiles
    if ismember(i, moveRecs)
        section = BImgs{i}(buf:h, buf:w); % take smaller part of image

        % normalize images
        section = section - mean(section(:));

        % Get x and y offset based on 2D fft cross corelation
        correl = xcorr2_fft(BImgRef,section);
        [ssr,snd] = max(correl(:));
        [ij, ji] = ind2sub(size(correl),snd);
        moveVert(i) = -(size(BImgRef,1) - ij - buf);
        moveHori(i) = -(size(BImgRef,2) - ji - buf);
    else
        moveVert(i) = 0;
        moveHori(i) = 0;
    end
end

% Calculate how much bigger the corrected image needs to be
moveVertMin = min(moveVert);
if moveVertMin>0, moveVertMin=1; end
moveHorMin = min(moveHori);
if moveHorMin>0, moveHorMin=1; end

% Calculate the offset of each image (final offset)
regiVert = zeros(nfiles, 2);
regiHori = zeros(nfiles, 2);
for i = 1:nfiles
    regiVert(i, 1) = moveVert(i) + abs(moveVertMin)+1;
    regiVert(i, 2) = regiVert(i, 1) + dims(i, 1)-1;
    regiHori(i, 1) = moveHori(i) + abs(moveHorMin)+1;
    regiHori(i, 2) = regiHori(i, 1) + dims(i, 2)-1;
end

% Apply offsets
for i = 1:nfiles
    imgtemp = BImgs{i};
    BImgs{i} = zeros(max(regiVert(:,2)), max(regiHori(:,2)));
    BImgs{i}(regiVert(i,1):regiVert(i,2), regiHori(i,1):regiHori(i,2)) = imgtemp;
    imgtemp = Masks{i};
    Masks{i} = zeros(max(regiVert(:,2)), max(regiHori(:,2)));
    Masks{i}(regiVert(i,1):regiVert(i,2), regiHori(i,1):regiHori(i,2)) = imgtemp;
    
    for j = 1:PPs(i).Cnt % Update coordinates of the contours
        PPs(i).Con(j).x = PPs(i).Con(j).x + regiHori(i,1)-1;
        PPs(i).Con(j).y = PPs(i).Con(j).y + regiVert(i,1)-1;
    end
    PPs(i).P(1,:) = PPs(i).P(1,:) + regiHori(i,1)-1; % And contour 'centres'
    PPs(i).P(2,:) = PPs(i).P(2,:) + regiVert(i,1)-1;
    
    fprintf('rec %d shift:%3d px horizontal,%3d px vertical\n',...
            i, regiHori(i,1)-1, regiVert(i,1)-1)
end

transformed(iter).vert = regiVert;
transformed(iter).hori = regiHori;

% Cut empty edges away
[BImgs, Masks, PPs, cutting] = CutEmptyEdges(BImgs, Masks, PPs, nfiles);

transformed(iter).cutting = cutting;

% Calculate correlation between the images
corrScore = BImgOverlapScore(BImgs);
corrScoreMean(:,end+1) = [sum(corrScore)./(nfiles-1)]'; %#ok<NBRAK>


% Automatic Rotation registration % % %
fprintf('rotating registration...\n')
rotations = -1.5:0.05:1.5; % rotations to apply: counterclockwise to clockwise degrees
buf = 10; % buffer image edges because rotation makes part fall off

rotCorrel  = zeros(nfiles, length(rotations));
BImgRef = BImgs{referenceNum}(buf:end-buf,buf:end-buf);

for i = 1:nfiles
    for j = 1:length(rotations)
        BImgRot = imrotate(BImgs{i}, rotations(j),'nearest','crop');
        BImgRot = BImgRot(buf:end-buf, buf:end-buf);
        
        rotCorrel(i,j)= corr2(BImgRot, BImgRef);
    end
end
rotCorrel = smoothG(rotCorrel', 2)';
toc

% Best rotations
[rotBest, rotIdx] = max(rotCorrel,[],2);
rotAng = rotations(rotIdx);
% difference in similarity between 0 rotation and best rotation
rotDiff = rotBest - rotCorrel(:,rotations==0);

% Apply rotations__________________________________________________________
rotAng(abs(rotAng)<0.02) = 0;
rotAng(lockRecs) = 0; % Do not rotate the recordings that shouldn't move

for i = find(abs(rotAng)>0.02) % If rotation angle < 0.02, don't rotate
    
    fprintf('rotating recording %d: %.2f degrees: %.3f better\n', i, rotAng(i), rotDiff(i))
    
    % Rotate image and Mask
    BImgs{i} = imrotate(BImgs{i}, rotAng(i),'bilinear','crop');
    Masks{i} = imrotate(Masks{i},rotAng(i),'nearest','crop');
    
    % Rotate contours coordinates
    PPs(i) = RotatePPCoordinates(PPs(i), -rotAng(i), size(BImgs{1}));
end

transformed(iter).rotation =  rotAng';

% Calculate similarity between the images
corrScore = BImgOverlapScore(BImgs);
corrScoreMean(:,end+1) = [sum(corrScore)./(nfiles-1)]'; %#ok<NBRAK>



% Show the result of registration  % % % % % % % % % % % % % % 
figure('Units', 'normalized', 'Position',[0.15 0.1 0.65 0.75])
% The registered images
subplot('Position', [0.04 0.05 0.65 0.95])
rgb = CreateRGB2(BImgs(toplot), colors);
imagesc(rgb);
hold on
for i = 1:PPs(referenceNum).Cnt
    plot(PPs(referenceNum).Con(i).x, PPs(referenceNum).Con(i).y, 'color', colors(referenceNum,:))
end
plot(PPs(referenceNum).P(1, :), PPs(referenceNum).P(2, :), '.', 'color', colors(referenceNum,:))
figtitle(sprintf('Registration iteration %d', iter))
for i = 1:length(toplot)
    text(20, 20+i*17, filenamesShort{toplot(i)},'color',colors(toplot(i),:),...
        'fontweight','bold','fontsize',12)
end

% Plot the change of correlation of the images over multiple
% registration steps
h = subplot('Position', [0.8 0.7 0.175 0.25]);
% h.PositionConstraint = 'outerposition';
step = (1:size(corrScoreMean,2))-1;
for i = 1:length(toplot)
    plot(step, corrScoreMean(toplot(i),:)', '-o', 'color', colors(toplot(i),:)); hold on
end
h.XTick = 1:size(corrScoreMean,2);
xLabels = cell(size(corrScoreMean,2),1);
for j = 1:size(corrScoreMean,2)
    if mod(j,2)==1; xLabels{j} = sprintf('trans. %d', floor(j/2)+1);
    else;           xLabels{j} = sprintf('rot. %d', floor(j/2)); end
end
h.XTickLabel = xLabels;
h.XTickLabelRotation = 45;
ylabel('correlation coefficient')
title('Average image correlation')
legend(filenamesShort, 'Location', 'southeast')

% Plot the correlation between all seperate background images
h = subplot('Position', [0.8 0.45 0.175 0.175]);
imagesc(corrScore);
caxis([0 1]); colorbar
h.XTick = 1:nfiles;
h.YTick = 1:nfiles;
h.XTickLabel = filenamesShort;
h.YTickLabel = filenamesShort;
h.XTickLabelRotation = 45;
title('Image correlation coefficients')

% Plot the correlation of different rotations
subplot('Position', [0.8 0.1 0.175 0.175])
% h.PositionConstraint = 'outerposition';
for i = 1:length(toplot)
    plot(rotations, rotCorrel(toplot(i),:),'.-','color',colors(toplot(i),:))
    hold on
end
plot(rotAng,rotBest,'xk')
title(sprintf('Image correlation with rec %d', referenceNum))
ylabel('Correlation coefficient r')
xlabel('Rotation (degrees)'), xlim([min(rotations),max(rotations)])

% Check for missing ROIs
for i = 1:nfiles
    missingROIs = find(diff(unique(Masks{i}))>1);
    if ~isempty(missingROIs)
        warning('%d ROIs has fallen off the FOV because of rotation! Ask Leander for better software', length(missingROIs))
    end
end

%% Re-draw all ROI contours.
% ROI contours are hard to register same as the images, since the contours
% can rotate on subpixel levels etc. One solution is to just redraw all
% ROI contours.

for i = 1:nfiles
    PPs(i) = PPfromMask(Masks{i}, [], PPs(i));
end

%% clear some registration variables
clearvars ssr step h j ji ij i regiHori regiVert RGB rotbest rotAng rotation
clearvars section w h xLabels snd rotIdx moveVert moveHori moveHorMin moveVertMin
clearvars imgtemp lockRecs moveRecs dims cutting buf referenceNum
clearvars rotCorrel rotBest rotDiff correl BImgRot BImgRef rotations

%% 
