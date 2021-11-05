% Register spectral images of many days, and link their ROIs.
% Registration uses 2D cross-correlation via fft, with a chosen reference
% 
% Run the script step by step
% 
% Leander de Kraker
% 
% 2018-2-2
% Last edit: 2020-2-24

clc
%% Load SPSIG data
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
    else
        pos = regexp(filepaths{i}, '\');
        
        % Try to extract the date from after the first _ in the filename
        try 
            % folder with the data should have recording date as name
            strdate = filepaths{i}(pos(end-1)+1:pos(end)-1);
            strdate = strrep(strrep(strdate, '-', ''), '_', '');
            filedates(i) = datetime(strdate, 'inputformat','yyyyMMdd');
            
        catch ME % Extraction failed
            if strcmp(ME.identifier, 'MATLAB:datetime:ParseErr') % Ask user for input
                questionStr = sprintf('No date detected for file %d, \nfile:%s\nfolder: %s\n',...
                               i, filenames{i},filepaths{i}); 
                questionStr = [questionStr, 'Please type its date here as yyyyMMdd:'];
                
                answerStr = inputdlg(questionStr);
                
                try % try to decipher user input as a datetime
                    filedates(i) = datetime(answerStr, 'inputformat','yyyyMMdd');
                    filedates(i)
                catch METOO % another error: No date for this recording.
                    warning('User input not recognized as date')
                    throw(METOO)
                end
            else % extraction failed because of unkown reason
                throw(ME)
            end
        end
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
    data(i) = load([filepaths{i} filenames{i}], 'SPic', 'Mask', 'PP', 'BImg');
    fprintf('\nloaded "%s"\n', filenames{i})
    
end
toc

clearvars selecting str pos i questionStr answerStr strdate

%% Backing up data variable WHEN GOING TO RESIZE SOME RECORDINGS
data2 = data;
%% Resize some recordings? ONLY IF RESIZE IS NECESSARY! Manual parameters!
% This step is only necessary when recordings are made with different zoom

data = data2;

rec = [1,2,3,4,5]; % Which recordings to resize

scalefac = [0 0]; % x stretch y stretch
xc = 313; % x Position of correctly registered cells
x1 = 351; % x position of neuron that is shifted in recording that will be resized 
x2 = 348;%  x position same cell as other in other recording
scalefac(1) = 1 + (x1-x2)/(xc-x2);
% % yc = 276; y1 = 422; y2 = 405;
% xc = 313; % x Position of correctly registered cells
% x1 = 353; % x position of neuron that is shifted in recording that will be resized 
% x2 = 348;% 
yc = 276; y1 = 422; y2 = 407;
scalefac(2) = 1 + (y1-y2)/(yc-y2);
% scalefac = [0.9125, 0.88];
% scalefac = [1.0075,1.0075];

% scalefac = [1.0725 1.1125]; % x stretch y stretch
fprintf('\n\nRECORDINGS ARE GETTING RESIZED!\n\n')
for i = 1:length(rec)
    outSize = round(scalefac .* [size(data(rec(i)).SPic,1), size(data(rec(i)).SPic,2)]);
    data(rec(i)).SPic = imresize(data(rec(i)).SPic, 'OutputSize', outSize,'method','bilinear');
    data(rec(i)).Mask = imresize(data(rec(i)).Mask, 'OutputSize', outSize([2 1]),'method','nearest');
    for j = 1:data(rec(i)).PP.Cnt
        data(rec(i)).PP.Con(j).x = data(rec(i)).PP.Con(j).x .* scalefac(1);
        data(rec(i)).PP.Con(j).y = data(rec(i)).PP.Con(j).y .* scalefac(2);
    end
    data(rec(i)).PP.P(1,:) = data(rec(i)).PP.P(1,:) .* scalefac(1);
    data(rec(i)).PP.P(2,:) = data(rec(i)).PP.P(2,:) .* scalefac(2);
    fprintf('recording %d: %s resized!\n', rec(i), filenamesShort{rec(i)})
end


%% Create cell arrays from the different datasets' SPic variables & Masks

% Preallocate variables
BImg = cell(1, nfiles);
Masks= cell(1, nfiles);
mini = zeros(1,nfiles);
maxi = zeros(1,nfiles);

for i = 1:nfiles
    %obtain highest contrasted spectral density
%     img = max(real(log(data(i).SPic)),[],3);
%     img = data(i).BImgA;
    img = data(i).BImg';

    % Normalize and remove -infs
    vals = unique(img);
    mini(i) = vals(2); % second smallest value (not -inf)
    maxi(i) = vals(end);
    img(img < mini(i)) = mini(i);
    img = (img - mini(i)) ./ range([mini(i), maxi(i)]);
    
    % Save mask and background image
    BImg{i} = img;
    Masks{i} = data(i).Mask';
end

clearvars mini maxi Img idx i vals channel Maskextra

% PLOTTING 
% Activate much tighter subplots
% [subplot margin top&side],[figure bottomspace,topspace],[leftspace,rightspace]
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.01], [0 0], [0 0]);

toplot = 1:nfiles; % sessions to plot

% Create Rainbowy color channel allocations for nfiles, with bias to blue,
colors = flipud(cmapL([0 0 1; 0 1 1; 0 1 0; 1 0.7 0; 1 0 0; 0.7 0 1], nfiles));
if size(colors,1)==3
    colors = [1 0 0; 0 1 0; 0.3 0.3 1];
elseif size(colors,1)==2
    colors = [1 0 0; 0 1 0.5];
end
RGB = permute(CreateRGB2(BImg(toplot), colors(toplot,:)), [2 1 3]);

figure; subplot(1,1,1)
imagesc(RGB)
for i = 1:length(toplot)
    text(20, 20+i*17, datestr(filedates(toplot(i))),'color',colors(toplot(i),:),...
        'fontweight','bold','fontsize',12)
end
title('non-registered masks')


%% On first run

BImg2 = BImg;
Masks2 = Masks;
try
    PP = [data.PP];
catch
    % probably crashing because of dissimilar structs because specprofile field
    % So let's leave those out
    PP = struct('P', [],  'A', [], 'Cnt', [], 'Con', []);
    for i = 1:nfiles
        PP(i).P = data(i).PP.P;
        PP(i).A = data(i).PP.A;
        PP(i).Cnt = data(i).PP.Cnt;
        PP(i).Con = data(i).PP.Con;
    end
end
transformed = struct();
transformed.xoff = {};
transformed.yoff = {};
transformed.rotation = {};
transformed.recordings = {};
corrScoreMean = [];

% Calculate similarity between the images
corrScore = BImgOverlapScore(BImg);
corrScoreMean(:,end+1) = [sum(corrScore)./(nfiles-1)]'; %#ok<NBRAK>

%% Register images with manual GUI.

BImg = BImg2;
Masks = Masks2;
imgs.BImg = BImg;
imgs.Masks = Masks;
interpol = {'bilinear', 'nearest'}; % How to interpolate each image for rotation?
ManualRegistrationv3(imgs, PP, filenamesShort, interpol)

%% Register images with manual ctrl point WARPING GUI.

BImg = BImg2;
Masks = Masks2;
filenamesShort = ShortenFileNames(filenames, 1, filedates);
imgs.BImg = BImg;
imgs.Masks = Masks;
interpol = {'bilinear', 'nearest'}; % How to interpolate each image for rotation?
ManualRegistration(imgs, PP, filenamesShort, interpol)



%% Register images automatically
referenceNum = 2; % which background image to take as reference
lockRecs = 1:3; % Lock certain recordings, do not move them
buf = 40; % Buffer to remove around the to-register image


moveRecs = 1:nfiles;
moveRecs(lockRecs) = [];

if referenceNum>nfiles
    error('Select a lower reference recording number! referenceNum exceeds number of files')
end

tic
BImg = BImg2;
Masks = Masks2;

fprintf('\nstarting registration. with recording %d as reference\n', referenceNum)
% Determine which minimum size image since recordings are different sizes
maxdim = [0 0];
dims = (cellfun(@size, BImg2,'UniformOutput',false));
dims = cat(1,dims{:});
maxdim(1) = max(dims(:,1));
maxdim(2) = max(dims(:,2));


w = min(dims(:,1))-buf;
h = min(dims(:,2))-buf;

x = zeros(1, nfiles); % x offset for images
y = zeros(1, nfiles);

% The base session to register others with
BImgRef = BImg2{referenceNum};
BImgRef = BImgRef - mean(BImgRef(:));

% Do cross correlations with the base image
for i = 1:nfiles
    if ismember(i, moveRecs)
        section = BImg2{i}(buf:w, buf:h); % take smaller part of image

        % normalize images
        section = section - mean(section(:));

        % Get x and y offset based on 2D fft cross corelation
        correl = xcorr2_fft(BImgRef,section);
        [ssr,snd] = max(correl(:));
        [ij, ji] = ind2sub(size(correl),snd);
        x(i) = -(size(BImgRef,1) - ij - buf);
        y(i) = -(size(BImgRef,2) - ji - buf);
    else
        x(i) = 0;
        y(i) = 0;
    end
end

% Calculate how much bigger the corrected image needs to be
xmin = min(x);
if xmin>0, xmin=1; end
ymin = min(y);
if ymin>0, ymin=1; end


% Calculate the offset of each image (final offset)
xoff = zeros(nfiles, 2);
yoff = zeros(nfiles, 2);
for i = 1:nfiles
    xoff(i, 1) = x(i) + abs(xmin)+1;
    xoff(i, 2) = xoff(i, 1) + dims(i, 1)-1;
    yoff(i, 1) = y(i) + abs(ymin)+1;
    yoff(i, 2) = yoff(i, 1) + dims(i, 2)-1;
end

% Apply offsets
for i = 1:nfiles
    BImg2{i} = zeros(max(xoff(:,2)), max(yoff(:,2)));
    BImg2{i}(xoff(i,1):xoff(i,2), yoff(i,1):yoff(i,2)) = BImg{i};
    Masks2{i} = zeros(max(xoff(:,2)), max(yoff(:,2)));
    Masks2{i}(xoff(i,1):xoff(i,2), yoff(i,1):yoff(i,2)) = Masks{i};
    
    for j = 1:PP(i).Cnt % Update coordinates of the contours
        PP(i).Con(j).x = PP(i).Con(j).x + xoff(i,1)-1;
        PP(i).Con(j).y = PP(i).Con(j).y + yoff(i,1)-1;
    end
    PP(i).P(1,:) = PP(i).P(1,:) + xoff(i,1)-1; % And contour 'centres'
    PP(i).P(2,:) = PP(i).P(2,:) + yoff(i,1)-1;
    
    fprintf('rec %d shift:%3d px horizontal,%3d px vertical\n',...
            i, xoff(i,1)-1, yoff(i,1)-1)
end

transformed.xoff{end+1} = xoff;
transformed.yoff{end+1} = yoff;

% Cut empty edges away
[BImg2, Masks2, PP] = CutEmptyEdges(BImg2, Masks2, PP, nfiles);


% Calculate correlation between the images
corrScore = BImgOverlapScore(BImg2);
corrScoreMean(:,end+1) = [sum(corrScore)./(nfiles-1)]'; %#ok<NBRAK>


% Automatic Rotation registration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('rotating registration...\n')

rotation = -1.5:0.05:1.5; % rotations to apply: counterclockwise to clockwise degrees
rotCorrel  = zeros(nfiles, length(rotation));
buf = 10; % buffer image edges because rotation makes part fall off

BImgRef = BImg2{referenceNum}(buf:end-buf,buf:end-buf);

for i = 1:nfiles
    for j = 1:length(rotation)
        BImgRot = imrotate(BImg2{i}, rotation(j),'bicubic','crop');
        BImgRot = BImgRot(buf:end-buf, buf:end-buf);
        
        rotCorrel(i,j)= corr2(BImgRot, BImgRef);
%       imagesc(permute(CreateRGB({BImgRef, BImgRot},'r g'),[2 1 3]));
%       title(sprintf('session %d. rotation %d: %.1fdeg. simil=%.3f',i,j,rotation(j),simil(i,j)))
    end
end
toc

% Best rotations
[rotBest, rotIdx] = max(rotCorrel,[],2);
rotAng = rotation(rotIdx);
% difference in similarity between 0 rotation and best rotation
rotDiff = rotBest - rotCorrel(:,rotation==0);

figure
% subplot(2,1,1)
for i = 1:length(toplot)
    plot(rotation, rotCorrel(toplot(i),:),'.-','color',colors(toplot(i),:))
    hold on
end
plot(rotAng,rotBest,'xk')
title(sprintf('image correlation with recording %d after different rotations buffer %d',...
                referenceNum, buf))
ylabel('correlation coefficient r')
xlabel('rotation (degrees)'), xlim([min(rotation),max(rotation)])
legend(filenamesShort)

% Apply rotations__________________________________________________________
% If difference between best and 0 is smaller than 0.01, don't bother
% If the difference is that small, it's probably because of the artifact
thres = 0.0005;
rotAng(rotDiff<thres) = 0;
rotAng(abs(rotAng)<0.03) = 0;
rotAng(lockRecs) = []; % Do not rotate the recordings that shouldn't move

for i = find(abs(rotAng)>0.03) % If rotation angle < 0.03, don't rotate
    
    fprintf('rotating recording %d: %.2f degrees: %.3f better\n', i, rotAng(i), rotDiff(i))
    
    % Rotate image and Mask
    BImg2{i} = imrotate(BImg2{i}, rotAng(i),'bilinear','crop');
    Masks2{i} = imrotate(Masks2{i},rotAng(i),'nearest','crop');
    
    % Rotate contours coordinates
    PP(i) = RotatePPCoordinates(PP(i), rotAng(i), size(BImg2{1}));
end

% Show the rotated images
RGB = CreateRGB2(BImg2(toplot), colors);
figure; subplot(1,1,1)
imagesc(permute(RGB, [2,1,3]));
title('registered background rotation corrected')
for i = 1:length(toplot)
    text(20, 20+i*17, filenamesShort{toplot(i)},'color',colors(toplot(i),:),...
        'fontweight','bold','fontsize',12)
end

transformed.rotation{end+1} =  rotAng;
transformed.recordings{end+1} = 1:nfiles;

% Calculate similarity between the images
corrScore = BImgOverlapScore(BImg2);
corrScoreMean(:,end+1) = [sum(corrScore)./(nfiles-1)]'; %#ok<NBRAK>


% Plot how well the background images overlap % % % % % % % % % % % % % % %
clearvars subplot

% Plot the change of correlation of the images over multiple
% registration steps
figure('units', 'normalized', 'OuterPosition', [0.4 0.3 0.3 0.5])
h = subplot(1,1,1);
step = (1:size(corrScoreMean,2))-1;
for i = 1:length(toplot)
    plot(step, corrScoreMean(toplot(i),:)', '-o', 'color', colors(toplot(i),:)); hold on
end
xlabel('transformation i')
ylabel('correlation coefficient')
legend(filenamesShort,'location', 'southeastoutside')
h.XTick = 1:size(corrScoreMean,2);
xLabels = cell(size(corrScoreMean,2),1);
for j = 1:size(corrScoreMean,2)
    if mod(j,2)==1; xLabels{j} = sprintf('trans. %d', floor(j/2)+1);
    else;           xLabels{j} = sprintf('rot. %d', floor(j/2)); end
end
h.XTickLabel = xLabels;
h.XTickLabelRotation = 45;

% The correlation between all seperate background images
figure('units', 'normalized', 'OuterPosition', [0 0.3 0.4 0.6]);
h = subplot(1,1,1);
imagesc(corrScore);
caxis([0 1]); colorbar
h.XTick = 1:nfiles;
h.YTick = 1:nfiles;
h.XTickLabel = filenamesShort;
h.YTickLabel = filenamesShort;
h.XTickLabelRotation = 45;
title('correlation coefficients')


%% Checking how ROIs overlap. % % % % % % % %  MATCHING ROIS

% Linking ROIs with each other, by checking if there is a certain percentage 
% overlap between them
thres = 0.525; % OVERLAP THRESHOLD range (0.5, 1) = 50% overlap - 100%

if thres > 1
    warning('thres should be equal or smaller then 1')
end
if thres <= 0
    warning('trheshold has to be above 0, and lower then 1 (100 %)')
elseif thres < 0.3
    warning('low threshold may introduce errors in the linked ROIs')
end

% Checking overlap between the ROIs
% inRoi:  column1: linked ROI of other mask
%         column3: number of pixels overlap
%         column4: percentage overlap of this ROI with other ROI
nmasks = 1:nfiles;
inRoi = {1,length(nmasks)};
fprintf('Calculating every ROI overlap from every recording to all recordings...\n')
for m = nmasks % all masks
    
    otherMasks = nmasks;
    % for each mask, save info for each ROI
    inRoi{m} = cell(PP(m).Cnt, length(otherMasks)); 
    for r = 1:PP(m).Cnt % for the ROIs in this mask
        
        roipos = Masks2{m}==r;
        for c = otherMasks % Check with all other Masks
            % Save which ROIs are in ROI n, how many pixels overlap and fraction overlap
            overlap = Masks2{c}(roipos);
            [vals, ~, idx] = unique(overlap); % which ROIs in this ROI?
            inRoi{m}{r,c} = [vals, accumarray(idx, 1)]; % number of pixels overlap
            inRoi{m}{r,c}(:,3) = inRoi{m}{r,c}(:,2)./sum(inRoi{m}{r,c}(:,2)); % frac overlap
            inRoi{m}{r,c}(inRoi{m}{r,c}==0,:) = []; % delete overlap with non-ROIs
        end
    end
end


% linked: column1: ROI of mask i, 
%         column2: linked ROI of other mask
%         column3: number of pixels overlap
%         column4: percentage overlap of this ROI with other ROI
fprintf('linking ROIs overlapping above threshold...\n')
linked = cell(nfiles);
for m = nmasks % m = current mask
    otherMasks = nmasks(nmasks~=m); % Skip own overlap
    for c = otherMasks % c = compared to
        for r = 1:PP(m).Cnt % r = each ROI
            link = inRoi{m}{r,c}(:,3)>thres;
            if any(link)
                thisLink = [ones(sum(link),1).*r, inRoi{m}{r,c}(link,:)];
                linked{m, c} = [linked{m, c}; thisLink];
            end
        end
    end
end

clearvars i j k m n c r roipos otherMasks thisLink link iteration vals idx

%__________________________________________________________________________
% Creating consistent links by looking at how much % the ROIs overlap,
% viewing from both the masks, instead of 1 towards the other

linked2 = cell(nfiles);
fprintf('Checking whether ROI overlap is high enough between the two recordings...\n')
% linked2: column1: ROI of recording m, 
%          column2: linked ROI of other recording
%          column3: number of pixels overlap
%          column4: percentage overlap of this ROI with other ROI
%          column5: percentage overlap for both the ROIs
for m = nmasks
    otherMasks = nmasks(nmasks~=m);
    for c = otherMasks
        for i = 1:size(linked{m,c},1)
            r = linked{m,c}(i,1); %  roi of mask
            rc = linked{m,c}(i,2); % r was linked with roi from other mask
            overlap1 = linked{m,c}(i,4);
            % Overlap looking from other recording
            found = inRoi{c}{rc, m}(inRoi{c}{rc, m}(:,1)==r, :); 
            
            overlap2 = found(3);
            overlap = (overlap1 + overlap2) ./ 2;
            if overlap > thres % true if mutual overlap is large enough!
                linked2{m,c} = [linked2{m,c}; [r, rc, overlap1, overlap2, overlap]];
            end
        end
    end
end
% Now combine the linked ROIs from recording 1-> recording 2 and 2->1 etc
% to one array
for c = nmasks
    otherMasks = c+1:nmasks(end);
    for m = otherMasks
        % Switch columns to match the order from other recording comparison
        if ~isempty(linked2{m,c})
            linked2{m,c} = [linked2{m,c}(:,2),linked2{m,c}(:,1),linked2{m,c}(:,4),linked2{m,c}(:,3),linked2{m,c}(:,5)];
            linked2{c,m} = [linked2{c,m}; linked2{m,c}];
            linked2{m,c} = [];
            [~, sorted]  = sort(linked2{c,m}(:,1));
            linked2{c,m} = linked2{c,m}(sorted,:);
            linked2{c,m} = unique(linked2{c,m},'rows');
        end
    end
end
clearvars found check overlap overlap1 overlap2 otherMasks i c r m

% all links in a matrix
fprintf('Compressing link data into match matrix...\n')
linkMat = cell(nfiles-1, 1);
linkMat2 = cell(nfiles-1, 1);
for m = 1:nfiles-1
    otherMasks = m+1:nmasks(end);
    linkMat{m} = zeros(1,nfiles);
    for c = otherMasks
        if ~isempty(linked2{m,c})
            own = linked2{m,c}(:,1);
            compared = linked2{m,c}(:,2);
            overlap = linked2{m,c}(:,5);
            linkMat{m}(end+1:end+length(own), [m, c]) = [own, compared];
        end
    end
    linkMat{m}(1,:) = [];
    % Sort to the neuron of current m
    [~, idx] = sort(linkMat{m}(:,m));
    linkMat{m} = linkMat{m}(idx,:);
    
    % Compress the matrix m so the neuron of mask m only appears once
    [rois, limits] = unique(linkMat{m}(:,m));
    limits(end+1) = size(linkMat{m},1)+1;
    linkMat2{m} = zeros(length(rois), nfiles);
    for j = 1:length(rois) % for each roi. put others in one row
        linkMat2{m}(j,:) = max(linkMat{m}(limits(j):limits(j+1)-1,:),[],1);
    end
end
linkMat2 = cell2mat(linkMat2);

% Sort & squish the matrix
for i = 2:nfiles
    [~,sorted] = sort(linkMat2(:,i));
    linkMat2 = linkMat2(sorted,:);
    [rois, limits] = unique(linkMat2(:,i));
    if length(limits) > 1
        limits(end+1) = size(linkMat2,1)+1;
        dif = diff(limits); % number of elements with ROI j
        % If one of the recordings matched another ROI keep it. So even if
        % the matched ROI is not found via some of the recordings, it can still
        % get linked with those recordings.
        for j = 2:length(rois)
            if dif(j) > 1
                rows = limits(j):(limits(j)+dif(j)-1);
                linkMat2(limits(j),:) = max(linkMat2(rows,:));
            end
        end
        % Now we made sure we don't lose any links, squish the matrix
        linkMat2 = linkMat2([1:limits(2), limits(3:end-1)'],:);
    end
end

% clearvars i j c otherMasks own compared idx rois limits

% Creating Mask with only linked ROIs

matchedMasks = Masks2;
for i = 1:nfiles
    idx = ismember(matchedMasks{i},linkMat2(:,i));
    matchedMasks{i}(~idx) = 0; 
end
% 
% matchedMasksm = CreateRGB(matchedMasks, toplot, colors, 'binary');
% matchedMasksm = permute(matchedMasksm, [2,1,3]);
% figure
% imagesc(matchedMasksm)
% title(sprintf('linked masks, links=%d, threshold = %.2f',size(linkMat2,1),thres2))
% 
clearvars i j r c m dif vals idx thisLink limits rois own compared otherMasks
fprintf('\ndone matching ROIs\n')
%
% figure
% imagesc(permute(CreateRGB(matchedMasks([1 2 3]), 'r g b', 'binary'), [2 1 3]))
% title('recordings 1 2 3')
% figure
% imagesc(permute(CreateRGB(matchedMasks([4 5 6]), 'r g b', 'binary'), [2 1 3]))
% title('recordings 4 5 6')



%% __________________________________________________________________________
% The 'OVERLAP SCORE' is based on the pixel overlap between the linked ROIs. 
% This also takes into account the ROIs which are not officially linked with
% each other, but were linked together because of mutual links.
score = zeros(length(linkMat2),1);
for i = 1:length(linkMat2)
    score(i) = OverlapScore(inRoi, linkMat2, i);
end

% Plot scores of the ROIs
h = figure; hold on
nLinks = sum(linkMat2~=0,2);
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
linkMat3 = linkMat2; % Extend linkMatrix with ROIs that didn't find any matches
for i = 1:nfiles
    % keep roi numbers which are not members of this linkmat column
    allRoisi = 1:PP(i).Cnt;
    singlesi = allRoisi(~ismember(allRoisi, linkMat2(:,i)));
    linkMat3(end+1:end+length(singlesi), i) = singlesi;
end
% extend scores with 0: there is no overlap for single neurons
score2 = score;
score2(end+1:length(linkMat3)) = 0; 

% calculate number of links and plot histogram
nLinks = sum(linkMat3~=0, 2); % number of roi links in each row

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
    present = nLinks(linkMat3(:,i)>0);
    % Count how many neurons with 1:nfiles matches recording i has
    nFoundBackMat(i,1:nfiles) = histc(present, 1:nfiles);
    nameCol{i} = sprintf('%d links',i);
end
nameCol{1} = 'single';
nFoundBackMat(nfiles+1,:) = sum(nFoundBackMat, 1);
nFoundBackMat(:,nfiles+1) = sum(nFoundBackMat, 2);

% Percentage wise confusion matrix
n = repmat([PP.Cnt]', [1, nfiles]);
n = [n; sum(n, 1)];
n = [n, repmat(n(nfiles+1,1),[nfiles+1,1])];
nFoundBackPercMat = round(nFoundBackMat ./ n .* 100);

figure('name','number of ROIs matched per recording')
h = uitable('Data',nFoundBackMat,'Units', 'Normalized', 'Position', [0, 0.6, 1, 0.3]);
h.ColumnName = [nameCol 'summed ROIs'];
h.RowName = [filenamesShort; {'summed ROIs'}];

h = uitable('Data',round(nFoundBackPercMat,1), 'Units','Normalized', 'Position', [0, 0.3, 1, 0.3]);
h.ColumnName = [nameCol, {'mean'}];
h.RowName = [filenamesShort; {'mean'}];


%
%__________________________________________________________________________
% TABLE 2
% Find how much ROIs recording i linked with recording j
confusionFoundMat = zeros(nfiles+1);
for i = 1:nfiles
    for j = 1:nfiles
        confusionFoundMat(i,j) = sum(linkMat3(:,i)>0 & linkMat3(:,j)>0);
    end
end
confusionFoundMat(nfiles+1,1:end-1) = round(mean(confusionFoundMat(1:end-1,1:end-1),1));
confusionFoundMat(1:end-1,nfiles+1) = round(mean(confusionFoundMat(1:end-1,1:end-1),2));

h = uitable('Data', confusionFoundMat, 'Units', 'Normalized', 'Position', [0, 0, 1, 0.3]);
h.ColumnName = [filenamesShort; {'mean'}];
h.RowName = [filenamesShort; {'mean'}];

% 
%__________________________________________________________________________
% Use number of ROIs connected to a ROI to create new image
nLinksMask = cell(1, nfiles);
for i = 2:nfiles
    rois = linkMat3(nLinks==i,:); % rois with i linked rois
    nLinksMask{i} = zeros(size(Masks2{i}));
    for j = 1:nfiles
        roij = rois(:,j);
        roij(roij==0) = [];
        idx = ismember(Masks2{j}, roij);
    	nLinksMask{i}(idx) = nLinksMask{i}(idx) + 1;
    end
end

RGB = CreateRGB2(nLinksMask(toplot), colors);
figure
imagesc(permute(RGB, [2 1 3]))  
text(20, 450, filenames{1}(1:14),'color',[1 1 1])
for i = 2:nfiles
    text(20, 20+i*10, sprintf('%d links', i), 'color',colors(i,:), 'fontweight','bold')
    text(20, 440+i*10, filenames{i}(1:14), 'color',[1 1 1])
end
title(sprintf('number of links=%d, threshold = %.2f',size(linkMat2,1),thres))

% Sort the matrix and scores so the most links, with the highest scores are
% on the top of the matrices
nLinks = sum(linkMat2 > 0, 2);
[~, sorted]=  sortrows([nLinks, score], 'descend');
linkMat2 = linkMat2(sorted,:);
score = score(sorted);
nLinks = nLinks(sorted);

% %% histograms of percentage overlap
% overlaps = {inRoi{1}{:,2}, inRoi{2}{:,1},inRoi{1}{:,3}, inRoi{2}{:,3}, inRoi{3}{:,1},inRoi{3}{:,2}}';
% overlaps = cell2mat(overlaps);
% overlaps = overlaps(:,3);
% [b, a] = hist(overlaps, 25);
% h1 = bar(a,b);
% h1.FaceColor = [0 0 1];
% h1.FaceAlpha = 0.8;

clearvars keep idx i j n matchedMasks h h1 RGB rois spaces

%% Chronic viewer checker UI
filenamesShort = ShortenFileNames(filenames, 1, filedates);

ChronicViewer(BImg2, Masks2, filenamesShort, nLinksMask, linkMat2, PP, score, inRoi)
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
save(filename, 'nfiles', 'BImg2', 'Masks2', 'PP', 'filenames', 'filepaths', 'filedates',...
    'transformed', 'thres', 'linked2', 'linkMat2', 'score', 'nLinks', 'nLinksMask', 'inRoi',...
    'nFoundBackMat','nFoundBackPercMat', 'confusionFoundMat')
fprintf('saved %s\n', filename)


