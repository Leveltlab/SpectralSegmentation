% Flatten a 3D mask
% Input: mask3d [3D double] (x, y, nROIs). Each ROI has a seperate image, 
%                              enabling weighing of pixels and ROI overlap
% 
% output: Mask [2D double]. (x, y) Each pixel is assigned to a ROIid
% 
% 
% 
% Leander de Kraker
% 2021-6-9
% 

% mask3d = compSpatial;
% mask3d = compSpatialSc;
d1 = size(mask3d, 1);
d2 = size(mask3d, 2);
nrois = size(mask3d, 3);


%% Flatten the mask

threshold = 0.03; % Weights lower than this get exluded from ROI mask.

mask3d_firstadded = cat(3, ones(d1, d2).*threshold, mask3d);
[BImg, Mask] = max(mask3d_firstadded, [], 3);
Mask = Mask - 1;

nroisNew = length(unique(Mask(:)))-1;
nroisMax = max(Mask(:));
missingRoisByMaskFlattening = find(diff(unique([Mask(:); 0]))>1);
fprintf('%d ROIs are now gone because flattening\n', length(missingRoisByMaskFlattening));
while nroisNew ~= nroisMax
    missingRoi = find(diff(unique([Mask(:); 0]))>1);
    Mask(Mask>missingRoi(1)) = Mask(Mask>missingRoi(1))-1;
    
    nroisMax = max(Mask(:));
end

%% If you have signals and ROIs were deleted! Update sig.
% sig = [oldnrois x time] -> [newnrois x time]
sig(missingRoisByMaskFlattening, :) = [];

%% Plot
figure
imagesc(Mask)
title(sprintf('threshold for mask at %.5f. Original nr ROIs: %d, new %d', threshold, nrois, nroisNew))
cmap = [0 0 0; jet(nroisNew)];
colormap(cmap)

figure
imagesc(BImg)
title('Representative image (weights of the masks): BImg')

%% Save mask to file

[filename, filepath] = uigetfile('*.mat','save the mask to a mat file');
save([filepath, filename], 'Mask','BImg','missingRoisByMaskFlattening', '-append')
fprintf('\nsaved Mask to %s\n\n', filename)

%% Get PP and stuff from Mask

PPfromMask


