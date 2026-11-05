function PP = RemoveROIsPP(PP, badROIs)
% Remove the ROIs indicated by badROIs 
% PP = RemoveROIsPP(PP, badROIs)
% 
% Input:
%   PP - The standard ROI information struct from SpecSeg
%   badROIs ([n x 1] double) - numbers indicating which ROIs should be
%                              deleted
% 
% Output:
%   PP - Updated ROI information struct
% 
% Leander de Kraker
% 2026-9-9
% See also: RemoveROIsMask, RemoveROIs, SortROIs
% 

PP.Con(badROIs) = [];
PP.A(badROIs) = [];
PP.P(:,badROIs) = [];
if isfield(PP, 'SpecProfile')
    PP.SpecProfile(:,badROIs) = [];
end
if isfield(PP, 'peakFreq')
    PP.peakFreq(badROIs) = [];
end
if isfield(PP, 'peakVal')
    PP.peakVal(badROIs)  = [];
end
if isfield(PP, 'Roundedness')
    PP.Roundedness(badROIs) = [];
end
if isfield(PP, 'Rvar')
    PP.Rvar(badROIs) = [];
end
if isfield(PP, 'creationMethod')
    PP.creationMethod(badROIs) = [];
end
PP.Cnt = size(PP.P, 2);
