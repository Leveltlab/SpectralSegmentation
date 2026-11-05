function [Mask, PP] = SortROIs(Mask, PP, idx)
% Resort the ROIs according to idx, for both Mask and PP
% [Mask, PP] = SortROIs(Mask, PP, idx)
% 
% Input:
%   Mask - a 2D ROI Mask.
%   PP - The standard ROI information struct from SpecSeg
%   idx ([n x 1] double) - new order of the ROIs: idx(i) is the OLD ROI
%                          number that becomes the NEW ROI number i.
% 
% Output:
%   Mask - Resorted Mask.
%   PP - Resorted ROI information struct.
% 
% Leander de Kraker
% 2026-9-9
% See also: SortROIsMask, SortROIsPP, RemoveROIs
% 

Mask = SortROIsMask(Mask, idx);
PP = SortROIsPP(PP, idx);
