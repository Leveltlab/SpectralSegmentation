function mask3d = Mask2D_To3D(mask2d)
% Convert a 2D ROI Mask to 3D
% 
% Input:
%   mask2d ([h x w] double): Mask of ROIs, for pixel with 0 indicating no 
%       ROI present at pixel, and a integer indicating the ROI id present
% 
% Output:
%   mask3d ([h x w x nrois] double): Mask of ROIs, with each ROI being put
%       in one of the slices of the 3D matrix
% 
% Leander de Kraker
% 2022-4-29
% 

dims = size(mask2d);
nrois = max(mask2d(:));
mask3d = zeros(dims(1), dims(2), nrois);
for i = 1:nrois
    mask3d(:,:,i) = mask2d==i;
end