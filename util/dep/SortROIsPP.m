function PP = SortROIsPP(PP, idx)
% Re-Sort the ROIs in Mask and PP
% [Mask, PP] = RemoveROIs(Mask, PP, badROIs)
% 
% Input:
%   PP - The standard ROI information struct from SpecSeg
%   idx ([n x 1] double) - numbers indicating how to sort the ROIs
% 
% Output:
%   PP - Updated ROI information struct
% 
% Leander de Kraker
% 2023-9-6
% 


PP.Con = PP.Con(idx);
PP.A = PP.A(idx);
PP.P = PP.P(:,idx);
if isfield(PP, 'SpecProfile')
    PP.SpecProfile = PP.SpecProfile(idx);
end
if isfield(PP, 'peakFreq')
    PP.peakFreq = PP.peakFreq(idx);
end
if isfield(PP, 'peakVal')
    PP.peakVal  = P.peakVal(idx);
end
if isfield(PP, 'Roundedness')
    PP.Roundedness = PP.Roundedness(idx);
end
if isfield(PP, 'Rvar')
    PP.Rvar = PP.Rvar(idx);
end
if isfield(PP, 'creationMethod')
    PP.creationMethod = PP.creationMethod(idx);
end
PP.Cnt = size(PP.P, 2);


