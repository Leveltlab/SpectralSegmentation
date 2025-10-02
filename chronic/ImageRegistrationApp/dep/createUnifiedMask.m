function [linkMat, unifiedMasks] = createUnifiedMask(filepaths, transforms)
    
    function masks = fetchMasks(filepaths)
        fetchStr = 'Mask';
        
        nFiles = length(filepaths);    
        masks = cell(size(filepaths));
        for i = 1:nFiles
            loadStruct = load(filepaths(i), fetchStr);
            masks{i} = loadStruct.(fetchStr);
        end
    end
    
    function unifiedMasks = linkMat2Masks(masks, linkmat)
        [n_Rois, nVars] = size(linkMat);
        if length(masks) ~= nVars
            error("inconsistent number of files between Masks array and linkmat")
        end
        
        doubled_mat = linkmat + n_Rois;
        doubled_mat(linkmat == 0) = 0;
        
        for i=1:nVars
            cur_mask = masks{i};
            doubled_mask = cur_mask + n_Rois;
            doubled_mask(cur_mask == 0) = 0;
            masks{i} = doubled_mask;
        end
        
        unifiedMasks = cell(size(masks));
        for j=1:nVars
            cur_mask = masks{j};
            for i=1:n_Rois
                cur_idx = doubled_mat(i,j);
                if cur_idx == 0
                    continue
                end
                cur_mask(cur_mask==cur_idx) = i;
            end
            unifiedMasks{j} = cur_mask;
        end
    end
    
    raw_masks = fetchMasks(filepaths);
    warped_masks = warpImagesfromTforms(transforms, raw_masks);
    
    thres = 0.6; % OVERLAP THRESHOLD.    range (0.5, 1) = 50% overlap - 100%
    
    % Checking overlap between the ROIs of all recording pairs
    inRoi = CalcRoiOverlap(warped_masks);
    
    % Linking overlapping ROIs between all recording pairs above threshold,
    [linked, linked1] = LinkOverlappingRois(inRoi, thres);
    
    % Squishing linked ROIs into one link matrix
    linkMat = SquishLinkedRois(linked);
    
    [n_Rois, nVars] = size(linkMat);
    
    is_valid = true; % flag to track validity
    for j = 1:nVars  % loop over each mask (column)
        col_vals = linkMat(:, j);
        col_vals(col_vals == 0) = []; % ignore zeros
        if numel(unique(col_vals)) ~= numel(col_vals)
            msg = printf('Duplicate ROI ID found in column %d\n', j);
            is_valid = false;
            error(msg);
        end
    end

    unifiedMasks = linkMat2Masks(raw_masks, linkMat);
end