function Mask = SortROIsMask(Mask, idx)
% Re-Sort the ROIs in Mask and PP
% [Mask, PP] = RemoveROIs(Mask, PP, badROIs)
% 
% Input:
%   Mask - a 2D ROI Mask.
%   idx ([n x 1] double) - numbers indicating how to sort the ROIs
% 
% Output:
%   Mask - Updated Mask.
% 
% Leander de Kraker
% 2023-9-6
% 


maskNew = zeros(size(Mask));
for i = 1:length(idx)
    maskNew(Mask==i) = idx(i); 
end


