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
RGB = CreateRGB2(BImgs(toplot), colors(toplot,:));
imagesc(RGB)
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
referenceNum = 3; % which background image to take as reference
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
RGB = CreateRGB2(BImgs(toplot), colors);
imagesc(RGB);
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


%% Checking how ROIs overlap. % % % % % % % %  MATCHING ROIS

% Linking ROIs with each other, by checking if there is a certain percentage 
% overlap between them
thres = 0.6; % OVERLAP THRESHOLD.    range (0.5, 1) = 50% overlap - 100%

if thres > 1
    warning('thres should be equal or smaller then 1')
end
if thres <= 0
    warning('trheshold has to be above 0, and lower then 1 (100 %)')
elseif thres < 0.3
    warning('low threshold may introduce errors in the linked ROIs')
end


% Checking overlap between the ROIs of all recording pairs
inRoi = CalcRoiOverlap(Masks);

% Linking overlapping ROIs between all recording pairs above threshold,
[linked, linked1] = LinkOverlappingRois(inRoi, thres);

% Squishing linked ROIs into one link matrix
linkMat = SquishLinkedRois(linked);

fprintf('\ndone matching ROIs\n')



%% MATCH STATISTIC SECTION % % % % % % % % % % % % % % % % % % % % % % % %

if ~exist('toplot','var'); toplot=1:nfiles;end
if ~exist('colors','var'); colors=flipud(cmapL([0 0 1; 0 1 1; 0 1 0; 1 0.7 0; 1 0 0; 0.7 0 1], nfiles));end

% Clean the linkMat of single ROI matches (those appear after editing in ChronicViewer)
nLinks = sum(linkMat~=0,2);
linkMat(nLinks==1,:)=[];
nLinks = sum(linkMat~=0,2);

nMatches = size(linkMat, 1);

% The 'OVERLAP SCORE' is based on the pixel overlap between the linked ROIs. 
% This also takes into account the ROIs which are not officially linked with
% each other, but were linked together because of mutual links.
score = zeros(nMatches,1);
for i = 1:nMatches
    score(i) = OverlapScore(inRoi, linkMat, i);
end

% Plot scores of the ROIs
h = figure; hold on
plot(nLinks+(rand(length(nLinks),1)./5-0.1), score, 'o')
% plot(nLinks+(score2-0.5), score2, 'o')
h.Children.XTick = 2:nfiles;
ylim([0 1])
title(sprintf('overlap scores for chronic matched ROIs. thr=%.2f',thres))
xlabel('ROI found in n recordings')
ylabel('score (% pixel overlap average)')
% Average score per number of links in the match
scoreMean = zeros(1,nfiles-1);
for i = 2:nfiles
    scoreMean(i-1) = mean(score(nLinks==i));
end
plot(2:nfiles, scoreMean,'ro-')
line([1.5 nfiles+0.5],[thres, thres],'color',[0 0 0 0.4],'linestyle','--')
legend({'match score','average score','threshold'},...
    'location','southeast')

%__________________________________________________________________________
linkMatAllRois = linkMat; % Extend linkMatrix with ROIs that didn't find any matches
for i = 1:nfiles
    % keep roi numbers which are not members of this linkmat column
    allRoisi = 1:PPs(i).Cnt;
    singlesi = allRoisi(~ismember(allRoisi, linkMat(:,i)));
    linkMatAllRois(end+1:end+length(singlesi), i) = singlesi;
end
% extend scores with 0: there is no overlap for single neurons
scoreAllRois = score;
scoreAllRois(end+1:length(linkMatAllRois)) = 0; 

% calculate number of links and plot histogram
nLinks = sum(linkMatAllRois~=0, 2); % number of roi links in each row

figure('Position', [105 755 560 222])
h = histogram(nLinks);
% h.BinCounts = h.BinCounts.*(1:nfiles);
title(sprintf('number of roi links in each row thres = %.2f, actual neuron count',... 
    thres))
ylabel('n neurons')
xlabel('found back')

% ________________________________________________________________________
% TABLE
% Find the number- and % of neurons in each recording that were found back
nFoundBackMat = zeros(nfiles+1);
nameCol = cell(1,nfiles);
for i = 1:nfiles
    present = nLinks(linkMatAllRois(:,i)>0);
    % Count how many neurons with 1:nfiles matches recording i has
    nFoundBackMat(i,1:nfiles) = histc(present, 1:nfiles);
    nameCol{i} = sprintf('%d links',i);
end
nameCol{1} = 'single';
nFoundBackMat(nfiles+1,:) = sum(nFoundBackMat, 1);
nFoundBackMat(:,nfiles+1) = sum(nFoundBackMat, 2);

% Percentage wise confusion matrix
n = double(repmat([PPs.Cnt]', [1, nfiles]));
n = [n; sum(n, 1)];
n = [n, repmat(n(nfiles+1,1),[nfiles+1,1])];
nFoundBackPercMat = round(nFoundBackMat ./ n .* 100);

figure('name','number of ROIs matched per recording')
h = uitable('Data',nFoundBackMat,'Units', 'Normalized', 'Position', [0, 0.5, 1, 0.4]);
h.ColumnName = [nameCol 'summed ROIs'];
h.RowName = [filenamesShort; {'summed ROIs'}];
annotation('textbox', [0, 0.95, 1, 0], 'string', 'Specific number of ROIs which are in a match with that amount of linked ROIs')

h = uitable('Data',round(nFoundBackPercMat,1), 'Units','Normalized', 'Position', [0, 0, 1, 0.4]);
h.ColumnName = [nameCol, {'mean'}];
h.RowName = [filenamesShort; {'mean'}];
annotation('textbox', [0, 0.45, 1, 0], 'string', 'Percentage of ROIs which are in a match with that amount of linked ROIs')


%__________________________________________________________________________
% TABLE figure 2. confusionlike matrix
% Find how much ROIs recording i linked with recording j
confusionFoundMat = zeros(nfiles+1);
for i = 1:nfiles
    for j = 1:nfiles
        confusionFoundMat(i,j) = sum(linkMatAllRois(:,i)>0 & linkMatAllRois(:,j)>0);
    end
end
confusionFoundMat(nfiles+1,1:end-1) = round(mean(confusionFoundMat(1:end-1,1:end-1),1));
confusionFoundMat(1:end-1,nfiles+1) = round(mean(confusionFoundMat(1:end-1,1:end-1),2));

figure
h = uitable('Data', confusionFoundMat, 'Units', 'Normalized', 'Position', [0, 0.5, 1, 0.4]);
h.ColumnName = [filenamesShort; {'mean'}];
h.RowName = [filenamesShort; {'mean'}];
annotation('textbox', [0, 0.95, 1, 0], 'string', 'Between which recordings ROIs are linked')

n = double(repmat([PPs.Cnt]', [1, nfiles]));
confusionFoundMatPerc = round(confusionFoundMat(1:nfiles, 1:nfiles) ./ n * 100);

h = uitable('Data', confusionFoundMatPerc, 'Units', 'Normalized', 'Position', [0, 0, 1, 0.4]);
h.ColumnName = filenamesShort;
h.RowName = filenamesShort;
annotation('textbox', [0, 0.45, 1, 0], 'string', 'Between which recordings ROIs are linked as percentage of recording row')

% 
%__________________________________________________________________________
% Use number of ROIs connected to a ROI to create new image
nLinksMask = cell(1, nfiles);
for i = 2:nfiles
    rois = linkMatAllRois(nLinks==i,:); % rois with i linked rois
    nLinksMask{i} = zeros(size(Masks{i}));
    for j = 1:nfiles
        roij = rois(:,j);
        roij(roij==0) = [];
        idx = ismember(Masks{j}, roij);
    	nLinksMask{i}(idx) = nLinksMask{i}(idx) + 1;
    end
end

RGB = CreateRGB2(nLinksMask(toplot), colors);
figure
imagesc(RGB)  
text(20, 450, filenames{1}(1:14),'color',[1 1 1])
for i = 2:nfiles
    text(20, 20+i*10, sprintf('%d links', i), 'color',colors(i,:), 'fontweight','bold')
    text(20, 440+i*10, filenames{i}(1:14), 'color',[1 1 1])
end
title(sprintf('number of links=%d, threshold = %.2f',size(linkMat,1),thres))

% Sort the matrix and scores so the most links, with the highest scores are
% on the top of the matrices
nLinks = sum(linkMat > 0, 2);
[~, sorted]=  sortrows([nLinks, score], 'descend');
linkMat = linkMat(sorted,:);
score = score(sorted);
nLinks = nLinks(sorted);

clearvars keep idx i j n matchedMasks h h1 RGB rois spaces sorted i nameCol nmasks roij
clearvars singlesi allRoisi

%% Chronic viewer checker UI

ChronicViewer(BImgs, Masks, filenamesShort, nLinksMask, linkMat, PPs, score, inRoi)
% If you saved changes in linkMat with the ChronicViewer, rerun previous
% section to recalculate statistics!!!!

%% Save chronic info
% Important chronic information gets saved
spaces = regexp(filenames{1}, '_');
filename = [filenames{1}(1:spaces(1)), filenames{1}(spaces(3)+1:spaces(4)), 'chronic.mat'];
[filename, pathname, ~] = uiputfile('', 'save the chronic information file', filename);

if exist('scalefac', 'var')
    transformed.scaleFac = scalefac;
    transformed.scaleRecs = rec;
end

if exist(strjoin({filename, '.mat'},''))~=0
    warning('file already exists, do you want to overwr.shit too late')
end

fprintf('saving\n')
cd(pathname)
save(filename, 'nfiles', 'BImgs', 'Masks', 'PPs', 'filenames', 'filenamesShort', 'filepaths', 'filedates',...
    'transformed', 'thres', 'inRoi', 'linked', 'linkMat', 'score', 'nLinks', 'nLinksMask',...
    'nFoundBackMat','nFoundBackPercMat', 'confusionFoundMat')
fprintf('saved %s\n', filename)


