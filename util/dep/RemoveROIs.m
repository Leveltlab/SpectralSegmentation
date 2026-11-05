function [Mask, PP] = RemoveROIs(Mask, PP, badROIs)
% Remove the ROIs indicated by badROIs 
% [Mask, PP] = RemoveROIs(Mask, PP, badROIs)
% 
% Input:
%   Mask - a 2D ROI Mask.
%   PP - The standard ROI information struct from SpecSeg
%   badROIs ([n x 1] double) - numbers indicating which ROIs should be
%                              deleted
% 
% Output:
%   Mask - Updated Mask.
%   PP - Updated ROI information struct
% 
% Leander de Kraker
% 2022-5-31
% 

Mask = RemoveROIsMask(Mask, badROIs);
PP = RemoveROIsPP(PP, badROIs);