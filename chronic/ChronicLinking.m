% The second part of the chronic matching process. Using the registered
% masks to calculate ROI overlap and matches
% 
% Leander de Kraker
% 2020-2-24
% 2025-9-24: Added ability to import registration_data
% 
%% Use output of the ImageReg2P app
if exist('registration_data', 'var')
    filepaths = cellstr(registration_data.Path);
    filenames = cellstr(registration_data.Name);
    filenamesShort = cellstr(registration_data.IDName);
    filedates = cat(1, registration_data.Date{:});
    colors = registration_data.imgColor;
    Masks = registration_data.Mask_warped;
    transformed.method = 'affine2D';
    transformed.transforms = registration_data.tForm;
    BImgs = registration_data.Image_warped;
    nfiles = length(filenames);
    
    PPs = struct('P', [],  'A', [], 'Cnt', [], 'Con', []);
    Masks = cell(nfiles, 1);
    for i = 1:nfiles
        load([filepaths{i}, filenames{i}], 'PP', 'Mask');
        PPs(i).A = PP.A;
        PPs(i).Cnt = PP.Cnt;
        PPs(i).Con = PP.Con;
        PPs(i).P = PP.P;
        Masks{i} = Mask;
    end
    Masks = warpImagesfromTforms(transformed.transforms, Masks);
    for i = 1:nfiles
        PPs(i) = PPfromMask(Masks{i}, [], PPs(i));
    end
end

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
nrois = [PPs.Cnt];

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
linkMatAllRois = ExtendLinkMat(linkMat, nrois);
% scoreAllRois   = ExtendScore(score, size(linkMatAllRois, 1));

% calculate number of links and plot histogram
nLinksAllRois = sum(linkMatAllRois~=0, 2); % number of roi links in each row

figure('Units', 'Normalized', 'Position', [0.05 0.7 0.3 0.2])
h = histogram(nLinksAllRois);
hax = gca;
set(hax, 'XTick', 1:nfiles)
title(sprintf('Number of ROI links in each row. Minimum overlap thres = %.2f',... 
    thres))
ylabel('n neurons')
xlabel('found back')


% ________________________________________________________________________
% TABLES
% Find the number- and % of neurons in each recording that were found back
[nFoundBack, nFoundBackPerc, nameCol] = CalcFoundBack(linkMatAllRois, nLinksAllRois, nrois);

figure('name','number of ROIs matched per recording')
h = uitable('Data',nFoundBack,'Units', 'Normalized', 'Position', [0, 0.5, 1, 0.4]);
h.ColumnName = [nameCol 'summed ROIs'];
h.RowName = [filenamesShort; {'summed ROIs'}];
annotation('textbox', [0, 0.95, 1, 0], 'string', 'Specific number of ROIs which are in a match with that amount of linked ROIs')

h = uitable('Data',round(nFoundBackPerc,1), 'Units','Normalized', 'Position', [0, 0, 1, 0.4]);
h.ColumnName = [nameCol, {'mean'}];
h.RowName = [filenamesShort; {'mean'}];
annotation('textbox', [0, 0.45, 1, 0], 'string', 'Percentage of ROIs which are in a match with that amount of linked ROIs')


%__________________________________________________________________________
% TABLE figure 2. confusionlike matrix
[confusionFoundMat, confusionFoundMatPerc] = CalcConfusionFoundMat(linkMatAllRois, nrois);

figure
h = uitable('Data', confusionFoundMat, 'Units', 'Normalized', 'Position', [0, 0.5, 1, 0.4]);
h.ColumnName = [filenamesShort; {'mean'}];
h.RowName = [filenamesShort; {'mean'}];
annotation('textbox', [0, 0.95, 1, 0], 'string', 'Between which recordings ROIs are linked')

h = uitable('Data', confusionFoundMatPerc, 'Units', 'Normalized', 'Position', [0, 0, 1, 0.4]);
h.ColumnName = filenamesShort;
h.RowName = filenamesShort;
annotation('textbox', [0, 0.45, 1, 0], 'string', 'Between which recordings ROIs are linked as percentage of recording row')

% 
%__________________________________________________________________________
% Use number of ROIs connected to a ROI to create new image
nLinksMask = CreateNLinksMask(Masks, linkMatAllRois, nLinksAllRois);

colors2 = hot(nfiles);
colors2 = [0.5 0.5 0.5; colors2(1:end-1,:)];

figure
rgb = CreateRGB2(nLinksMask, colors2);
imagesc(rgb)
for i = 1:nfiles
    text(20, 20+i*17, sprintf('%d links', i), 'color',colors2(i,:), 'fontweight','bold')
    text(size(rgb,2)-10, size(rgb,1)-(nfiles+1)*15+i*15, filenamesShort{i}, 'color', [1 1 1], 'HorizontalAlignment','right')
end
title(sprintf('number of links=%d, threshold = %.2f', size(linkMat,1), thres))


%__________________________________________________________________________
% Sort the matrix and scores so the most links, with the highest scores are
% on the top of the matrices
nLinks = sum(linkMat > 0, 2);
[~, sorted]=  sortrows([nLinks, score], 'descend');
linkMat = linkMat(sorted,:);
score = score(sorted);
nLinks = nLinks(sorted);

clearvars i j h h1  rois sorted i nameCol

%% Chronic viewer checker UI

ChronicViewer(BImgs, Masks, filenamesShort, linkMat, PPs, score, inRoi)
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
    'nFoundBack','nFoundBackPerc', 'confusionFoundMat')
fprintf('saved %s\n', filename)


