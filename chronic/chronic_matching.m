% Register spectral images of many days, and link their ROIs.
% Registration uses 2D cross-correlation with a chosen reference dataset
% 
% Leander de Kraker
% 
% 2018-2-2
% Last edit: 2019-9-12


%% Load SPSIG data
% load multiple files, ordered on recording time. Press cancel when done

filenames = {};
filepaths = {};
filedate = NaT(1);

selecting = true; % true as long as files are being selected
i = 0;
while selecting
    i = i + 1;
    if i == 1
        str = sprintf('load spectrally analyzed file %d', i);
    else
        str = sprintf('load spectrally analyzed file %d. previous: %s. Press cancel when done'...
            , i, datestr(filedate(i-1)));
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
            filedate(i) = datetime(filepaths{i}(pos(end-1)+1:pos(end)-1), 'inputformat','yyyyMMdd');
            
        catch ME % Extraction failed
            if strcmp(ME.identifier, 'MATLAB:datetime:ParseErr') % Ask user for input
                questionStr = sprintf('No date detected for file %d, \nfile:%s\nfolder: %s\n',...
                               i, filenames{i},filepaths{i}); 
                questionStr = [questionStr, 'Please type its date here as yyyyMMdd:'];
                
                answerStr = inputdlg(questionStr);
                
                try % try to decipher user input as a datetime
                    filedate(i) = datetime(answerStr, 'inputformat','yyyyMMdd');
                    filedate(i)
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
filedate = filedate';
nfiles = length(filenames); % The number of files that have been selected

tic
clearvars data
fprintf('\nloading in %d files:\n', nfiles)
for i = 1:nfiles
    fprintf('\nloading file %d...',i)
    data(i) = load([filepaths{i} filenames{i}], 'SPic', 'Sax', 'freq', 'Mask', 'PP', 'BImg');
%     data(i) = load([filepaths{i} filenames{i}], 'SPic', 'Sax', 'freq', 'Mask', 'PP', 'BImg', 'BImgAverage');
    fprintf('\nloaded "%s"\n', filenames{i})
    
end
toc

clearvars selecting str pos i questionStr answerStr

%% Backing up data variable
data2 = data;
%% Resize some recordings?
% This step is only necessary when recordings are made with different zoom
% levels.

data = data2;

rec = [2,3,4];

% scalefac = [0 0]; % x stretch y stretch
% xc = 259; % x Position of correctly registered cells
% x1 = 58; % x position of neuron that is shifted in recording that will be resized 
% x2 = 81;%  x position same cell as other in other recording
% scalefac(1) = 1 + (x1-x2)/(xc-x2);
% yc = 151;
% y1 = 470;
% y2 = 442;
% scalefac(2) = 1 + (y1-y2)/(yc-y2);
scalefac = [0.9125, 0.88];
% scalefac = [1.0075,1.0075];

% scalefac = [1.0725 1.1125]; % x stretch y stretch
fprintf('\n\nRECORDINGS ARE GETTING RESIZED!\n\n')
for i = 1:length(rec)
    outSize = round(scalefac .* [size(data(rec(i)).SPic,1), size(data(rec(i)).SPic,2)]);
    data(rec(i)).SPic = imresize(data(rec(i)).SPic, 'OutputSize', outSize,'method','nearest');
%     data(rec(i)).SPic(data(rec(i)).SPic<0) = 0; % interpolation causes values below zero which 
    data(rec(i)).Mask = imresize(data(rec(i)).Mask, 'OutputSize', outSize([2 1]),'method','nearest');
    for j = 1:data(rec(i)).PP.Cnt
        data(rec(i)).PP.Con(j).x = data(rec(i)).PP.Con(j).x .* scalefac(1);
        data(rec(i)).PP.Con(j).y = data(rec(i)).PP.Con(j).y .* scalefac(2);
    end
    data(rec(i)).PP.P(1,:) = data(rec(i)).PP.P(1,:) .* scalefac(1);
    data(rec(i)).PP.P(2,:) = data(rec(i)).PP.P(2,:) .* scalefac(2);
    fprintf('recording %d: %s resized!\n', rec(i), datestr(filedate{rec(i)}))
end


%% Create cell arrays from the different datasets' SPic variables & Masks

% Preallocate variables
BImg = cell(1, nfiles);
Masks= cell(1, nfiles);
mini = zeros(1,nfiles);
maxi = zeros(1,nfiles);

for i = 1:nfiles
    %obtain highest contrasted spectral density
    Img = max(real(log(data(i).SPic)),[],3);
%     Img = data(i).BImgA;
%     Img = real(log(data(i).SPic(:,:,1))); % Average fluoresence picture

    % Normalize and remove -infs
    vals = unique(Img);
    mini(i) = vals(2); % second smallest value (not -inf)
    maxi(i) = vals(end);
    Img(Img < mini(i)) = mini(i);
    Img = (Img - mini(i)) ./ range([mini(i), maxi(i)]);
    
    % Save mask and background image
    BImg{i} = Img;
    Masks{i} = data(i).Mask';
end

clearvars mini maxi Img idx i vals channel Maskextra

%% Plot

% Activate much tighter subplots
% [subplot margin top&side],[figure bottomspace,topspace],[leftspace,rightspace]
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.01], [0 0], [0 0]);

toplot = 1:nfiles; % sessions to plot
% Create Rainbowy color channel allocations for nfiles, with bias to blue,
% green, red unmixed allocation: 'r', 'g' and 'b'
% colors = strsplit('b b bg g g gr r r rb'); 
% chosenColors = round(linspace(1,length(colors),nfiles));
% colors = colors(chosenColors);
% % Because rb is the last color it will always get chosen, but I don't want
% % that, I would prefer red if it isn't present already
% if ~ismember(colors, 'r') 
%     colors{end} = 'r';
% end
% % toplot = [1 2 3 4];
% % colors = 'r g b gb';
% 
% figure; subplot(1,1,1)
% [RGB, cvals] = CreateRGB(BImg, colors);
% imagesc(permute(RGB,[2 1 3]))
% for i = toplot
%     text(20, 20+i*17, datestr(filedate(i)),'color',cvals(i,:),...
%         'fontweight','bold','fontsize',12)
% end
% 
% figure; subplot(1,1,1)
% imagesc(permute(CreateRGB(Masks, colors, 'binary'),[2 1 3]))
% title('non-registered masks')

colors = flipud(cmapL([0 0 1; 0 1 1; 0 1 0; 1 0.7 0; 1 0 0; 0.7 0 1], length(toplot)));
% colors = flipud(cmapL([0 0 1; 0 1 1; 0 1 0; 1 0.7 0; 1 0 0], length(toplot)));
if length(toplot)==3
    colors = [1 0 0; 0 1 0; 0.3 0.3 1];
elseif length(toplot)==2
    colors = [1 0 0; 0 1 0.5];
end
RGB = CreateRGB2(BImg, toplot, colors);

figure; subplot(1,1,1)
imagesc(RGB)
for i = toplot
    text(20, 20+i*17, datestr(filedate(i)),'color',colors(i,:),...
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


%% Register images
referenceNum = 3; % which background image to take as reference

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


% buffersize = 100;
buf = 16;

w = dims(1,1)-buf;
h = dims(1,2)-buf;

x = zeros(1, nfiles); % x offset for images
y = zeros(1, nfiles);

% The base session to register others with
base = BImg2{referenceNum};
base = base - mean(base(:));

% Do cross correlations with the base image
for i = 1:nfiles
    section = BImg2{i}(buf:w, buf:h); % take smaller part of image
    
    % normalize images
    section = section - mean(section(:));
    
    % Get x and y offset based on 2D cross corelation
    correl = xcorr2(base,section);
    [ssr,snd] = max(correl(:));
    [ij, ji] = ind2sub(size(correl),snd);
    x(i) = -(size(base,1) - ij - buf);
    y(i) = -(size(base,2) - ji - buf);
end

% Calculate how much bigger the corrected image needs to be
xmin = min(x);
if xmin>0, xmin=1; end
ymin = min(y);
if ymin>0, ymin=1; end


% Calculate the offset of each image (final offset)
xoff = zeros(nfiles, 2);
yoff = zeros(nfiles, 2);
% xoff(1, 1) = abs(xmin)+1;
% xoff(1, 2) = dims(1, 1)+xoff(1,1)-1;
% yoff(1, 1) = abs(ymin)+1;
% yoff(1, 2) = dims(1, 2)+yoff(1,1)-1;
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

clearvars h w i ij ji Freq correl x y ymax ymin xmin xmax ssr snd extra

% Calculate similarity between the images
corrScore = BImgOverlapScore(BImg2);
corrScoreMean(:,end+1) = [sum(corrScore)./(nfiles-1)]'; %#ok<NBRAK>


% Automatic Rotation registration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('rotating registration...\n')

rotation = -1.75:0.05:1.75; % rotations to apply: counterclockwise to clockwise degrees
simil  = zeros(nfiles, length(rotation));
% correl = zeros(nfiles, length(rotation));
buf = 10; % buffer image edges because rotation makes part fall off

BImgRef = BImg2{referenceNum}(buf:end-buf,buf:end-buf);

for i = 1:nfiles
    for j = 1:length(rotation)
        BImgRot = imrotate(BImg2{i}, rotation(j),'bicubic','crop');
        BImgRot = BImgRot(buf:end-buf, buf:end-buf);
        
        simil(i,j)= corr2(BImgRot, BImgRef);
%         imagesc(permute(CreateRGB({BImgRef, BImgRot},'r g'),[2 1 3]));
%         subplot(1,1,1)        
%         imagesc(simmap)
%         title(sprintf('session %d. rotation %d: %.1fdeg. simil=%.3f',...
%                         i,j,rotation(j),simil(i,j)))
%         pause(0.1)
    end
end
toc

% Best rotations
[rotBest, rotIdx] = max(simil,[],2);
rotAng = rotation(rotIdx);
% difference in similarity between 0 rotation and best rotation
rotDiff = rotBest - simil(:,rotation==0);

figure
% subplot(2,1,1)
% plot(rotation, simil', '.-')
for i = 1:size(simil(:,1))
    plot(rotation, simil(i,:),'.-','color',colors(i,:))
    hold on
end
plot(rotAng,rotBest,'xk')
title(sprintf('image similarity with recording 1 after different rotations buffer %d',buf))
ylabel('correlation coefficient r')
xlabel('rotation (degrees)'), xlim([min(rotation),max(rotation)])
legend(datestr(filedate,'yyyy-mm-dd'))

% Apply rotations__________________________________________________________
% If difference between best and 0 is smaller than 0.01, don't bother
% If the difference is that small, it's probably because of the artifact
thres = 0.0005;
rotAng(rotDiff<thres) = 0;
rotAng(abs(rotAng)<0.03) = 0;
dims = size(BImg2{1});
shift = dims./2; % Amount to shift contour coordinates

for i = find(abs(rotAng)>0.03) % If rotation angle < 0.03, don't rotate
    
    fprintf('rotating recording %d: %.2f degrees: %.3f better\n', i, rotAng(i), rotDiff(i))
    
    % Rotate image and Mask
    BImg2{i} = imrotate(BImg2{i}, rotAng(i),'bilinear','crop');
    Masks2{i} = imrotate(Masks2{i},rotAng(i),'nearest','crop');
    
    % Rotate contours coordinates
    rotRad = deg2rad(rotAng(i));
    for j = 1:PP(i).Cnt
        % shift the coordinate points so the centre coordinate is 0
        x = PP(i).Con(j).x - shift(1); % Adjust the contour coordinates
        y = PP(i).Con(j).y - shift(2);
        % Rotate the shifted coordinates (around the new centre [0, 0])
        x =  x  * cos(rotRad) - y  * sin(rotRad);
        y =  x  * sin(rotRad) + y  * cos(rotRad);
        % Shift them back where they should be
        PP(i).Con(j).x = x + shift(1);
        PP(i).Con(j).y = y + shift(2);
    end
    % Also adjust the contour 'centres'
    xp = PP(i).P(1,:) - shift(1);
    yp = PP(i).P(2,:) - shift(2);
    xp = xp * cos(rotRad) - yp * sin(rotRad);
    yp = xp * sin(rotRad) + yp * cos(rotRad);
    PP(i).P(1,:) = xp + shift(1);
    PP(i).P(2,:) = yp + shift(2);
end

% Show the rotated images
RGB = CreateRGB2(BImg2, toplot, colors);
figure; subplot(1,1,1)
imagesc(permute(RGB, [2,1,3]));
title('registered background rotation corrected')
for i = 1:nfiles
    text(20, 20+i*17, datestr(filedate(i), 'yyyy-mm-dd'),'color',colors(i,:),...
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
for i = 1:nfiles
    plot(step, corrScoreMean(i,:)', '-o', 'color', colors(i,:))
    hold on
end
xlabel('transformation i')
ylabel('correlation coefficient')
legend(datestr(filedate, 'yyyy-mm-dd'),'location', 'southeastoutside')
h.XTick = 1:size(corrScoreMean,2);
xLabels = cell(size(corrScoreMean,2),1);
for j = 1:size(corrScoreMean,2)
    if mod(j,2)==1
        xLabels{j} = sprintf('trans. %d', floor(j/2)+1);
    else
        xLabels{j} = sprintf('rot. %d', floor(j/2));
    end
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
h.XTickLabel = datestr(filedate, 'yyyy-mm-dd');
h.YTickLabel = datestr(filedate, 'yyyy-mm-dd');
h.XTickLabelRotation = 45;
title('correlation coefficients')

%% Cut empty edges away from the BImgs and Mask % % % % % % % % % % % % % % 
%% backup the background images before removing the empty edges
BImgBackup = BImg2;
MasksBackup = Masks2;
PPBackup = PP;

%% Use the backup to revert the background imagesjj
BImg2 = BImgBackup;
Masks2 = MasksBackup;
PP = PPBackup;
% Determine which minimum size image since recordings are different sizes
maxdim = [0 0];
dims = (cellfun(@size, BImg2,'UniformOutput',false));
dims = cat(1,dims{:});
maxdim(1) = max(dims(:,1));
maxdim(2) = max(dims(:,2));

%% Cut empty edges away from the BImgs and Mask % EXPERIMENTAL
% This code needs to detect the edges in order to not crash.
BImgMat = cat(3, BImg2{:});
BImgMat = max(BImgMat,[],3);
MaskMat = cat(3, Masks2{:});
MaskMat = mean(MaskMat>0,3);

edges1 = find(diff(all(BImgMat==0))~=0);
edges2 = find(diff(all(BImgMat'==0))~=0);

figure
imagesc(CreateRGB({BImgMat,MaskMat},'g br'));
title('do not "cut the edges" if the lines do not denote the edges of the image')
hold on
plot([edges1(1) edges1(1)],[1 dims(1)],'y')
plot([edges1(2) edges1(2)],[1 dims(1)],'y')
plot([1 dims(2)], [edges2(1) edges2(1)],'w')
plot([1 dims(2)], [edges2(2) edges2(2)],'w')

%% cut the edges
% This code needs to detect the edges in order to not crash.
for i = 1:nfiles
    BImg2{i}(edges2(2):end,:) = [];
    Masks2{i}(edges2(2):end,:) = [];
    BImg2{i}(:,edges1(2):end) = [];
    Masks2{i}(:,edges1(2):end) = [];
    BImg2{i}(1:edges2(1)-1,:) = []; % remove the top and left side
    Masks2{i}(1:edges2(1)-1,:) = [];
    BImg2{i}(:,1:edges1(1)-1) = []; % remove the bottom and right side
    Masks2{i}(:,1:edges1(1)-1) = [];
    for j = 1:PP(i).Cnt
        PP(i).Con(j).y = PP(i).Con(j).y - edges1(1);
        PP(i).Con(j).x = PP(i).Con(j).x - edges2(1);
    end
    PP(i).P(2,:) = PP(i).P(2,:) - edges1(1);
    PP(i).P(2,:) = PP(i).P(3,:) - edges2(1);
end

% Determine which minimum size image since recordings are different sizes
maxdim = [0 0];
dims = (cellfun(@size, BImg2,'UniformOutput',false));
dims = cat(1,dims{:});
maxdim(1) = max(dims(:,1));
maxdim(2) = max(dims(:,2));

%% Checking how ROIs overlap. % % % % % % % %  MATCHING ROIS

% Linking ROIs with each other, by just checking if there is a certain % 
% overlap between them
thres = 0.51; % overlap needed

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
clearvars found check overlap overlap1 overlap2 toSave otherMasks i c r m

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
    limits(end+1) = length(linkMat2)+1;
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
clearvars i j r c m idx limits rois own compared otherMasks
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

figure
h = histogram(nLinks);
% h.BinCounts = h.BinCounts.*(1:nfiles);
title(sprintf('number of roi links in each row thres = %.2f, actual neuron count',... 
    thres))
ylabel('n neurons')
xlabel('found back')

% ________________________________________________________________________
% TABLE
% Find the number- and % of neurons in each recording that were found back
confusion = zeros(nfiles+1);
nameCol = cell(1,nfiles+1);
nameRow = cell(1,nfiles+1);
for i = 1:nfiles
    present = nLinks(linkMat3(:,i)>0);
    % Count how many neurons with 1:nfiles matches recording i has
    confusion(i,1:nfiles) = histc(present, 1:nfiles);
    spaces = regexp(filenames{i}, '_');
    nameRow{i} = [filenames{i}(1:spaces(1)-1), ' ', filenames{i}(spaces(1)+1:spaces(2)-1)];
    nameCol{i} = sprintf('%d links',i);
end
nameCol{1} = 'single';
nameCol{end} = 'summed ROIs';
nameRow{end} = 'summed ROIs';
confusion(nfiles+1,:) = sum(confusion, 1);
confusion(:,nfiles+1) = sum(confusion, 2);

% Percentage wise confusion matrix
n = repmat([PP.Cnt]', [1, nfiles]);
n = [n; sum(n, 1)];
n = [n, repmat(n(nfiles+1,1),[nfiles+1,1])];
confusionPerc = round(confusion ./ n .* 100);

figure('name','number of ROIs matched per recording')
h = uitable('Data',confusion,'Units', 'Normalized', 'Position',[0, 0.5, 1, 0.5]);
h.ColumnName = nameCol;
h.RowName = nameRow;

h = uitable('Data',round(confusionPerc,1), 'Units','Normalized', 'Position',[0, 0, 1, 0.5]);
h.ColumnName = nameCol;
h.RowName = nameRow;


%
%__________________________________________________________________________
% TABLE 2
% Find how much ROIs recording i linked with recording j
confusionMore = zeros(nfiles+1);
for i = 1:nfiles
    for j = 1:nfiles
        confusionMore(i,j) = sum(linkMat3(:,i)>0 & linkMat3(:,j)>0);
    end
end
confusionMore(nfiles+1,1:end-1) = round(mean(confusionMore(1:end-1,1:end-1),1));
confusionMore(1:end-1,nfiles+1) = round(mean(confusionMore(1:end-1,1:end-1),2));

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

RGB = CreateRGB2(nLinksMask, toplot, colors);
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

clearvars keep idx i matchedMasks

% %% histograms of % overlap
% 
% overlaps = {inRoi{1}{:,2}, inRoi{2}{:,1},inRoi{1}{:,3}, inRoi{2}{:,3}, inRoi{3}{:,1},inRoi{3}{:,2}}';
% overlaps = cell2mat(overlaps);
% overlaps = overlaps(:,3);
% 
% [b, a] = hist(overlaps, 25);
% h1 = bar(a,b);
% h1.FaceColor = [0 0 1];
% h1.FaceAlpha = 0.8;


%% Chronic viewer checker UI
ChronicViewer(BImg2, Masks2, filenames, nLinksMask, linkMat2, PP, score, inRoi)


%% Save chronic info
% Important chronic information gets saved
spaces = regexp(filenames{1}, '_');
filename = [filenames{1}(1:spaces(1)), filenames{1}(spaces(3)+1:spaces(4)), 'chronic.mat'];
[filename, pathname, ~] = uiputfile('', 'save the chronic information file', filename);

filenames = filenames;
filepaths = filepaths;

if exist('scalefac')
    transformed.scaleFac = scalefac;
    transformed.scaleRecs = rec;
end

if exist(strjoin({filename, '.mat'},''))~=0
    warning('file already exists, do you want to overwr.shit too late')
end

fprintf('saving\n')
cd(pathname)
save(filename, 'nfiles', 'BImg2', 'Masks2', 'PP', 'filenames', 'filepaths', 'transformed',...
    'thres', 'linked2', 'linkMat2', 'score', 'nLinksMask', 'inRoi', 'filedate',...
    'confusion','confusionMore')
fprintf('saved %s\n', filename)


%% Load chronic info

[filename] = uigetfile('');
load(filename)

if exist('filedate', 'var')
    % Convert filedate variable to proper datetime array
    if iscell(filedate)
        filedate = [filedate{:}];
        fprintf('converted filedates from cell array to datetime array\n')
    end
else
    fprintf('o dear no file dates were extracted from the titles\n')
end

nLinks = sum(linkMat2~=0, 2); % number of roi links in each row

clearvars keep idx i j roij idx rois RGB matchedMasks


%% Register other images in the files like the spectral BImgs
% Other images could be present in SPSIG files. This will register
% requested variables, using the same offsets and rotations as done with
% the spectral BImgs
%
%

whichone = 1;

interImgs = cell(nfiles,1);
hasImg = zeros(nfiles,1);
for i = 1:nfiles
    if exist([filepaths{i} filenames{i}],'file')
        presentVars = who('-file',[filepaths{i} filenames{i}]);
        if ismember('BImgA',presentVars) && whichone == 1
            hasImg(i) = 1;
            load([filepaths{i} filenames{i}], 'BImgA')
            tempImg = BImgA;
            tempImg = (tempImg - min(tempImg(:))) / (max(tempImg(:)-min(tempImg(:))));
            lims = prctile(tempImg, [1 99.5]);
            tempImg(tempImg<lims(1)) = lims(1);
            tempImg(tempImg>lims(2)) = lims(2);
    %         interImgs{i} = tempImg';
            interImgs{i} = log1p(tempImg)';
    %         if size(data(i).SPic,1)~=size(interImgs{i},1) || size(data(i).SPic,2)~=size(interImgs{i},2)
    %             fprintf('\nSize of BImgA not the same as expected from SPic for rec %d', i)
    %             size(data(i).SPic)
    %             size(interImgs{i})
    %         end
            clearvars BImgA
            
        elseif ismember('BImgMax',presentVars) && whichone == 2
                hasImg(i) = 1;
                load([filepaths{i} filenames{i}], 'BImgMax')
                tempImg = BImgMax;
                tempImg = (tempImg - min(tempImg(:))) / (max(tempImg(:)-min(tempImg(:))));
                lims = prctile(tempImg, [1 99.5]);
                tempImg(tempImg<lims(1)) = lims(1);
                tempImg(tempImg>lims(2)) = lims(2);
        %         interImgs{i} = tempImg';
                interImgs{i} = log1p(tempImg)';
        %         if size(data(i).SPic,1)~=size(interImgs{i},1) || size(data(i).SPic,2)~=size(interImgs{i},2)
        %             fprintf('\nSize of BImgA not the same as expected from SPic for rec %d', i)
        %             size(data(i).SPic)
        %             size(interImgs{i})
        %         end
                clearvars BImgA

        elseif ismember('BImgAverage',presentVars) && whichone == 3
            hasImg(i) = 1;
            load([filepaths{i} filenames{i}], 'BImgAverage')
            tempImg = BImgAverage;
            tempImg = (tempImg - min(tempImg(:))) / (max(tempImg(:)-min(tempImg(:))));
            lims = prctile(tempImg, [1 99.5]);
            tempImg(tempImg<lims(1)) = lims(1);
            tempImg(tempImg>lims(2)) = lims(2);
            interImgs{i} = tempImg';
    %         interImgs{i} = log1p(tempImg)';
    %         if size(data(i).SPic,1)~=size(interImgs{i},1) || size(data(i).SPic,2)~=size(interImgs{i},2)
    %             fprintf('\nSize of BImgA not the same as expected from SPic for rec %d', i)
    %             size(data(i).SPic)
    %             size(interImgs{i})
    %         end
            clearvars BImgAverage
        else
            fprintf('No Img found for recording %d\n', i)
        end
    else
        fprintf('file %d not found\n', i)
    end
end

for i = 1:nfiles
    if hasImg(i)>0
        for j = 1:length(transformed.xoff)
            interImgsTemp = interImgs;
            xoff = transformed.xoff{j};
            yoff = transformed.yoff{j};
            rotation = transformed.rotation{j};
            
            % Apply the translation and rotation
            interImgs{i} = zeros(max(xoff(:,2)), max(yoff(:,2)));
            interImgs{i}(xoff(i,1):xoff(i,2), yoff(i,1):yoff(i,2)) = interImgsTemp{i};
            if rotation(i)>0.03
                interImgs{i} = imrotate(interImgs{i}, rotation(i),'bilinear','crop');
            end
        end
    end
end

% Show the registered images
toplot = find(hasImg>0);
colors = flipud(cmapL([0 0 1; 0 1 1; 0 1 0; 1 0.7 0; 1 0 0; 0.7 0 1], length(toplot)));

% toplot(1:2:end) = [];
RGB = CreateRGB2(interImgs, toplot, colors(toplot,:));

resizing = 2;
RGB = imresize(RGB, 2);

RGB = permute(RGB, [2,1,3]);
figure; subplot(1,1,1)
imagesc(RGB);
title('registered background rotation corrected')
for i = toplot'
    text(20, (20+i*17)*resizing, datestr(filedate(i)),'color',colors(i,:),...
        'fontweight','bold','fontsize',12)
end


%% Match red cells
% The Res file should have a field 'red' in 'info.rois'
%
%

interROIs = cell(1,nfiles); % Interesting ROIs/ stained interneuron ROIs
hasID = false(1,nfiles);
for i = 1:nfiles
    load([filepaths{i} filenames{i}(1:end-4), '_Res.mat'],'info')
    if isfield(info.rois, 'red')
        hasID(i) = true;
        interROIs{i} = find([info.rois(:).red]);
    else
        fprintf('No interest ROIs for recording %d\n', i)
    end
end
hasID = find(hasID); % So which recordings had redROIs identified

linkMatinterROI = false(size(linkMat2));
for i = 1:length(hasID)
    reci = hasID(i);
    linkMatinterROI(:,reci) = ismember(linkMat2(:,reci),interROIs{reci});
end
confidence = sum(linkMatinterROI,2)./nLinks;
idx = find(confidence);
values = confidence(idx);


ChronicViewer(BImg2, Masks2, filenames, nLinksMask, linkMat2, PP, [score confidence], inRoi)


