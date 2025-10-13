function inRoi = CalcRoiOverlap(Masks)
% Checking overlap between every ROI of a recording to each other recording
% 
% output:
% inRoi ([1 x nfiles] cell with inside[nrois_i x nfiles] cell):  
%         column1: linked ROI of other mask
%         column2: number of pixels overlap
%         column3: percentage overlap of this ROI with other ROI
% 

nfiles = length(Masks);
nmasks = 1:nfiles;
inRoi = cell(1,nfiles);

for m = nmasks % all masks
    % for each mask, save info for each ROI
    otherMasks = nmasks;
    nrois = max(Masks{m}(:));
    inRoi{m} = cell(nrois, length(otherMasks));
    
    for r = 1:nrois % for the ROIs in this mask
        roipos = Masks{m}==r;
        
        for c = otherMasks % Check with all other Masks
            % Save which ROIs are in ROI n, how many pixels overlap and fraction overlap
            overlap = Masks{c}(roipos);
            [vals, ~, idx] = unique(overlap); % which ROIs in this ROI?
            inRoi{m}{r,c} = [vals, accumarray(idx, 1)]; % number of pixels overlap
            inRoi{m}{r,c}(:,3) = inRoi{m}{r,c}(:,2)./sum(inRoi{m}{r,c}(:,2)); % frac overlap
            inRoi{m}{r,c}(inRoi{m}{r,c}==0,:) = []; % delete overlap with non-ROIs
        end
    end
end