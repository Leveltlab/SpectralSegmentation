% Determine if ROIs seem to be missed for valid reason using all possible
% images of a recording
% 
% Step 1. Do chronic matching
% Step 2. Run MoveLikeChronic to save imgs with all possible image types of
% the recordings
%  
% 
% Leander de Kraker
% 2023-10-18
%

[chronicName, chronicPath] = uigetfile('*chronic.mat','get Chronic file');

%% Load the chronic file
load([chronicPath, chronicName], 'filepaths', 'filenames', 'filedates', 'nfiles',...
                                 'linkMat', 'transformed', 'nLinks', 'imgs', 'PPs', 'Masks')
imgTypes = fieldnames(imgs);
if ~exist('imgs','var');fprintf('Run MoveLikeChronic\n');end

%% For which individiual background types to check if ROI is visible?

imgCheck = {'BImg', 'BImgAverage', 'BImgMax'};
nimgCheck = length(imgCheck);
dims = size(imgs.(imgTypes{1}){1});

imgsRGB = cell(nfiles, 1);
for i = 1:nfiles
    imgsRGB{i} = zeros(dims(1), dims(2), 3);
    for j = 1:nimgCheck
        img = imgs.(imgCheck{j}){i};
        vals = unique(img(:));
        valCutoff = prctile(vals(:), 2); % darken values
        img(img<valCutoff) = valCutoff;
        valCutOff = prctile(vals(:), 99.9); % brighten values
        img(img>valCutOff) = valCutOff;
        
        img = rescale(img);
        img = MidtoneBalance(img, 0.3);
        imgsRGB{i}(:,:,j) = img;
    end
end


%% 

from = 3; % from which recordings to select the existing ROIs
to = 1;   % Which recording to note existence of (missing) ROIs
doMatched = false; % Also note for matched?

existing = nan(size(linkMat, 1), size(linkMat, 2), nimgCheck);

clc
% close all
figure('WindowStyle', 'docked')
pos = [0.05 0.05 0.9 0.9];
buffer = [0.1 0 0.1 0.1];
titles = cellstr(datestr(filedates));
cutbuffer = 35;
distr = 'rect';

if doMatched % show for all of the from recording
    roisToDo = find(linkMat(:,from)>0)';
else % Only show ROIs which are not linked in the other recordings
    roisToDo = find(linkMat(:,from)>0 & ~all(linkMat(:,to),2))';
end


[hAxes, hImgs] = PlotManyImgs(imgsRGB, pos, distr, buffer, titles, []);
colormap(gray)

colors = repmat([0 0 1 0.5], [nfiles, 1]);
for i = 1:nfiles
    axes(hAxes(i))
    hold on
    PlotCon(PPs(i), 'r')
    colorsi = colors;
    colorsi(i, :) = [1 0 0 1];
    hlines = PlotChronicConnectedROIs(linkMat, PPs, colorsi);
end

hDot = plot(1,1,'x');
hDot2 = plot(1,1,'x');
for i = roisToDo([29:end])
    
    idx = find(linkMat(i, :));
    roi = linkMat(i, idx);
    titles = cell(nfiles, 1);
    xlims = [];
    ylims = [];
    x = []; % for average ROI seedpoint coordinate
    y = []; 
    delete(hDot)
    hDot = gobjects(length(idx), 1);
    for j = 1:length(idx)
        xlims = cat(2, xlims, PPs(idx(j)).Con(roi(j)).x);
        ylims = cat(2, ylims, PPs(idx(j)).Con(roi(j)).y);
        roiID = linkMat(i,idx(j));
        xj = PPs(idx(j)).P(1,roi(j));
        yj = PPs(idx(j)).P(2,roi(j));
        x = [x xj]; 
        y = [y yj];
        hDot(j) = plot(xj, yj, 'bx', 'Parent', hAxes(idx(j)));
    end
    
    cutx = [min(xlims)-cutbuffer, max(xlims)+cutbuffer];
    cuty = [min(ylims)-cutbuffer, max(ylims)+cutbuffer];
    cutx(cutx<1) = 1; cutx(cutx>dims(2)) = dims(2);
    cuty(cuty<1) = 1; cuty(cuty>dims(1)) = dims(1);
    
    x = mean(x);
    y = mean(y);
    
    set(hAxes, 'XLim', cutx, 'YLim', cuty)
    linkaxes(hAxes, 'xy')

    delete(hDot2)
    hDot2 = gobjects(nfiles, 1);
    for j = 1:nfiles
        hDot2(j) = plot(x, y, 'bo', 'Parent', hAxes(j));
    end
%     set(hAxes, 'XTick', [], 'YTick', [])
    
    % Going through images to check if ROI is visually present
    roiStr = sprintf('ROI (%d/%d) ', find(roisToDo==i), length(roisToDo));
    fprintf(roiStr)
    for j = 1:nimgCheck
        for f = 1:nfiles
            hImgs(f).CData = imgsRGB{f}(:,:,j);
        end
        title(hAxes(from), [roiStr, imgCheck{j}])
        fprintf('%s', imgCheck{j})
        for f = to %[to from]
            title(hAxes(f), 'Neuron present? (Y/N)')
            answer = 'x';
            while ~contains('yn', answer)
                answer = lower(input(sprintf(' %d: ', f), 's'));
            end
            if contains(answer,'y')
                existing(i,to,j) = 1;
            elseif contains(answer, 'n')
                existing(i,to,j) = 0;
            else
                fprintf('This should not be possible\n')
                existing(i,to,j) = -1;
            end
        end
    end
end

%% Make summary

% existingSmall: for each checked ROI (1st dimension), was it labeled as not visible (0 or 1) 
%      in each of the used image types (2nd dimension)?
% y: how many ROIs were labeled as either not visible (first colom), or
% visible (second colom). Third colom is present if ROI was labeled as -1

existingSmall = existing;
existingSmall(:, all(all(isnan(existing),3)),:) = [];
existingSmall(all(all(isnan(existing),3),2),:,:) = [];
existingSmall = squeeze(existingSmall);
y = zeros(nimgCheck, 3);
meanings = [0 1 -1];
for j = 1:nimgCheck
    for i = 1:length(meanings)
        y(j, i) = sum(existing(:,to(1),j)==meanings(i));
    end
end
if sum(y(:,3))>0
    fprintf('some ROIs are labeled as -1 which should not happen. Try to edit these entries\n')
else
    y(:,3) = [];
end

%% Save

save([chronicName(1:end-4) '_missingVisual'], 'existing', 'from', 'imgCheck', 'existingSmall', 'y')


