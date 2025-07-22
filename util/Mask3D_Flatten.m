function [Mask, maskImg, missingROIs, presentROIs] = Mask3D_Flatten(mask3d, threshold, doPlot)
% Flatten a 3D mask. After you have the mask you can create ROI information
% variable PP by running PPfromMask
% 
% Input: 
%   - mask3d (3D double [x, y, nROIs]): Each ROI has a seperate image, 
%                              enabling weighing of pixels and ROI overlap
%   - threshold (scalar double): Weights lower than this get exluded from 
%                               ROI mask.
%   - doPlot (boolean): true=plot results (default), false=don't plot
% 
% Output: 
%   - Mask (2D double, [x, y]): Each pixel is assigned to a ROIid.
%   - maskImg (2D double, [x, y]): weights of the ROIs image.
%   - missingROIs ([1D double]): which ROI ids are missing,
%       because it is possible that some ROIs were completely overlapping
%       some other ROIs, making it impossible to have them in the 2D mask.
%   - presentROIs ([1D double]): Which ROI IDs are present in the new mask.
% 
% If you have missing ROIs by mask flattening. Make sure to delete info
% related to the ROIs in other places if you continue with this mask: 
% % sig = [oldnrois x time] -> [newnrois x time]
% sig(missingROIs, :) = [];
% 
% 
% Leander de Kraker
% 2021-6-9
% 2025-7-16: changed into a function
% 
arguments
    mask3d
    threshold (1,1) double = 0.03
    doPlot (1,1) logical = true
end

d1 = size(mask3d, 1);
d2 = size(mask3d, 2);
nrois3D = size(mask3d, 3);


%% Flatten the mask

mask3d_firstadded = cat(3, ones(d1, d2).*threshold, mask3d);
[maskImg, Mask] = max(mask3d_firstadded, [], 3);
maskImg = maskImg-threshold;
Mask = Mask - 1;

nroisNew = length(unique(Mask(:)))-1;
nroisMax = max(unique(Mask(:)));
presentROIs = unique(Mask(:));
presentROIs(presentROIs==0) = [];
missingROIs = 1:nrois3D;
missingROIs(presentROIs) = [];
fprintf('%d ROIs missing by flattening mask\n', length(missingROIs))
if ~isempty(missingROIs)
    while nroisNew ~= nroisMax
        missingRoi = find(diff(unique([Mask(:); 0]))>1);
        Mask(Mask>missingRoi(1)) = Mask(Mask>missingRoi(1))-1;
        nroisMax = max(Mask(:));
    end
end


%% Plot
if doPlot
    figure
    imagesc(mean(mask3d, 3))
    title('original mean projection of mask')

    figure
    imagesc(Mask)
    title(sprintf('threshold for mask at %.5f. Original nr ROIs: %d, new %d', threshold, nrois3D, nroisNew))
    cmap = [0 0 0; jet(nroisNew)];
    colormap(cmap); 
    h = colorbar;
    h.Ticks = unique(round(h.Ticks));
    ylabel(h, 'ROI ID')
    
    figure
    imagesc(maskImg)
    title('Representative image (weights of the masks): maskImg')
    h = colorbar;
    ylabel(h, 'ROI weight value')
end


