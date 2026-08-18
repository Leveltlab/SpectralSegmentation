function [unifiedMasks, unifiedMask, roiSource] = mergeChronicROIs(Masks, PPs, linkMat, transforms, filepaths, filenames, pixelAspectRatio, options)
% 
% MERGECHRONICROIS combine the per-recording ROI masks of a chronic
% recording into a single mask + PP struct
% 
% load('myfile_chronic.mat', 'Masks', 'PPs', 'linkMat')
% option.minAreaPixels = 30;
% [unifiedMask, unifiedPP, roiSource] = mergeChronicROIs(Masks, PPs, linkMat, options);
% 
% Matched ROIs (rows of linkMat, each with >= 2 nonzero entries) are
% fused into one averaged shape:
%
% ROIs that were never matched to anything else are added unchanged -
% ROIs whose evidence falls below minEvidence are dropped entirely
%
% Claude Sonnet 5 Medium & Leander de Kraker
% 2026-7-22
%

arguments
    Masks (:, 1) cell
    PPs (1, :) struct
    linkMat (:, :) double
    transforms affinetform2d
    filepaths
    filenames
    pixelAspectRatio (1, 1) double
    options.smoothSigma double = 1.5    % px, blur applied before averaging matched shapes
    options.shapeThreshold double = 0.5 % fraction of matched recordings a pixel must appear in to survive
    options.minAreaPixels double = 20    % drop stray fragments smaller than this after thresholding
    options.minEvidence double = 0      % final ROIs with less total evidence than this are discarded
    options.showEvidencePlot logical = true % print stats + plot the evidence distribution before filtering
    options.maskVarName = 'Mask';
end

nfiles = length(Masks);
if length(PPs) ~= nfiles
    error("PPs has %d entries but %d Masks were given", length(PPs), nfiles)
end
if size(linkMat, 2) ~= nfiles
    error("linkMat has %d columns but %d recordings (Masks) were given", size(linkMat, 2), nfiles)
end
imgSize = size(Masks{1});

%% 1. Assign every (recording, ROI) pair a "final" merged ROI ID
nMatchedRows = size(linkMat, 1);
assign = arrayfun(@(i) zeros(PPs(i).Cnt, 1), 1:nfiles, 'UniformOutput', false); % assign{rec}(roi) = finalID

for row = 1:nMatchedRows
    for rec = 1:nfiles
        roi = linkMat(row, rec);
        if roi > 0
            assign{rec}(roi) = row;
        end
    end
end

nextID = nMatchedRows;
for rec = 1:nfiles
    unassigned = find(assign{rec} == 0)';
    for roi = unassigned
        nextID = nextID + 1;
        assign{rec}(roi) = nextID;
    end
end
nFinal = nextID;

nLinks = sum(linkMat>0, 2);
nLinks(end+1:nFinal) = 1;

%% 2. Collect, for every final ID, which (recording, roi) pairs feed it
roiSource = repmat(struct('recordings', [], 'rois', []), nFinal, 1);
for rec = 1:nfiles
    for roi = 1:PPs(rec).Cnt
        id = assign{rec}(roi);
        roiSource(id).recordings(end+1) = rec;
        roiSource(id).rois(end+1) = roi;
    end
end

%% 3. Score each final ROI's evidence: more/better-supporting
%     recordings = more trustworthy, regardless of shape
hasRvar = isfield(PPs, 'Rvar');
evidence = zeros(nFinal, 1);
for id = 1:nFinal
    recs = roiSource(id).recordings;
    rois = roiSource(id).rois;
    if hasRvar
        evidence(id) = sum(arrayfun(@(c) PPs(recs(c)).Rvar(rois(c)), 1:numel(recs)));
    else
        evidence(id) = numel(recs); % no quality score available - fall back to recording count
    end
end

%% 4. Show the evidence distribution before anything gets dropped,
%     so minEvidence can be picked by looking at real numbers instead
%     of guessing blind
if options.showEvidencePlot
    isSingleton = nLinks==1;
    if hasRvar
        evidenceLabel = 'evidence (sum of Rvar across contributing recordings)';
        evidenceKind  = 'sum of Rvar';
        evidenceNorm = evidence;
    else
        evidenceLabel = 'evidence (recording count)';
        evidenceKind  = 'recording count';
        evidenceNorm = evidence;
    end
    
    fprintf('mergeChronicROIs: evidence stats (%s)\n', evidenceKind);
    fprintf('  all       n=%4d  median=%.3f  [%.3f, %.3f]\n', ...
        nFinal, median(evidenceNorm), min(evidenceNorm), max(evidenceNorm));
    fprintf('  singleton n=%4d  median=%.3f  [%.3f, %.3f]\n', ...
        nnz(isSingleton), median(evidenceNorm(isSingleton)), min(evidenceNorm(isSingleton)), max(evidenceNorm(isSingleton)));
    fprintf('  merged    n=%4d  median=%.3f  [%.3f, %.3f]\n', ...
        nnz(~isSingleton), median(evidenceNorm(~isSingleton)), min(evidenceNorm(~isSingleton)), max(evidenceNorm(~isSingleton)));
    
    figure('Name', 'mergeChronicROIs: evidence distribution');
    hold on
    bins = 0:0.05:max(evidenceNorm);
    histogram(evidenceNorm(isSingleton), bins, 'FaceColor', [0.85 0.35 0.1], 'DisplayName', 'singleton (1 recording)');
    histogram(evidenceNorm(~isSingleton), bins, 'FaceColor', [0.1 0.45 0.8], 'DisplayName', 'merged (>=2 recordings)');
    if options.minEvidence > 0
        xline(options.minEvidence, 'k--', 'minEvidence', 'LineWidth', 1.5, 'DisplayName', 'current cutoff');
    end
    xlabel(evidenceLabel)
    ylabel('number of final ROIs')
    legend('Location', 'best')
    title('Pick minEvidence by where singleton/merged separate')
    hold off
end

%% 5. Drop low-evidence ROIs entirely and renumber what's left
keep = find(evidence >= options.minEvidence);
nDropped = nFinal - numel(keep);
roiSource = roiSource(keep);
evidence = evidence(keep);
nFinal = numel(keep);

%% 6. Build each surviving ROI's shape, then paint it into one mask,
%     letting the higher-evidence shape win any pixel conflicts
unifiedMask = zeros(imgSize);
confidence  = zeros(imgSize);
nMerged = 0;

for id = 1:nFinal
    recs = roiSource(id).recordings;
    rois = roiSource(id).rois;
    k = numel(recs);
    
    if k == 1
        shape = Masks{recs(1)} == rois(1);
        conf  = double(shape) .* evidence(id);
    else % Blur masks together to get average shape
        nMerged = nMerged + 1;
        summed = zeros(imgSize);
        for c = 1:k
            bw = Masks{recs(c)} == rois(c);
            summed = summed + imgaussfilt(double(bw), options.smoothSigma);
        end
        avgProb = summed ./ k;
        shape = avgProb >= options.shapeThreshold;
        shape = bwareaopen(shape, options.minAreaPixels);
        conf = avgProb .* evidence(id); % spatial agreement scaled by how well-evidenced this ROI is overall
    end
    
    px = shape & (conf > confidence); % To see if this ROI has the highest confidence on each pixel where it could be
    unifiedMask(px) = id;
    confidence(px) = conf(px); % Update the found confidences for ROI presence
end

fprintf('mergeChronicROIs: %d final ROIs (%d merged from >=2 recordings, %d kept as-is, %d dropped below minEvidence)\n', ...
    nFinal, nMerged, nFinal - nMerged, nDropped);


%% unwarp1. get original size, and rebuild the exact output view it was originally warped into
originalSizes = cell(nfiles, 1);
outViews    = cell(nfiles, 1);
for i = 1:nfiles
    m = matfile(fullfile(filepaths{i}, filenames{i})); % metadata only, doesn't load the array
    originalSizes{i} = size(m, options.maskVarName);
    outViews{i} = affineOutputView(originalSizes{i}, transforms(i), "BoundsStyle", "FollowOutput");
end

%% unwarp2. Reconstruct the shared canvas (finalView) the unified mask lives on
if ~isfield(transforms, 'finalView')
    finalView = transforms.finalView;
else % old version didn't save finalView
    finalView = getLargestImgRef(vertcat(outViews{:}));
end

%% unwarp3. Invert each recording's transform, warping the unified mask from
%     the shared canvas back onto that recording's own native grid
unifiedMasks = cell(nfiles, 1);
identity = affinetform2d();
ownView = cell(nfiles, 1);
for i = 1:nfiles
    ownView{i} = affineOutputView(originalSizes{i}, identity, "BoundsStyle", "SameAsInput");
    unifiedMasks{i} = imwarp(unifiedMask, finalView, invert(transforms(i)), ...
        "OutputView", ownView{i}, "interp", 'nearest');
end

%% Save new SPSIG files for the unifiedMask
filenamesBase = strrep(filenames, '_SPSIG.mat', '');
for i = 1:nfiles
    filenamesUnified = [filenamesBase{i}, '_unified_SPSIG.mat'];
    unifiedFilei.PP = PPFromMask(unifiedMasks{i});
    unifiedFilei.Mask = unifiedMasks{i};
    unifiedFilei.transforms.ownView = ownView{i};
    save([filepaths{i}, filenamesUnified{i}], '-struct', 'unifiedFilei', '-append')
end
% Calculate and update the PP
doSave = true;
for i = 1:nfiles
    unifiedPPi = PPModernize(filepaths{i}, filenamesUnified{i}, doSave);
    if i==1
        unifiedPPs = unifiedPPi;
    else
        unifiedPPs(i) = unifiedPPi;
    end
end

%% 
figure
colors = jet(nfiles);
imagesc(CreateRGB2(Masks, colors))
hold on
for i = 1:nfiles
    PlotCon(PPs(i), colors(i,:));
end
title('original masks')
figure
imagesc(CreateRGB2(unifiedMasks, colors))
hold on
title('unified masks')
for i = 1:nfiles
    PlotCon(unifiedPPs(i), colors(i,:))
end
