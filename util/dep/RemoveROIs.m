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


Mask(ismember(Mask, badROIs)) = 0; % Set all badROIs to 0
vals = unique(Mask(:)); % which ROIs are left in the Mask
vals = vals(2:end); % remove the 0       
% Update the ROI ID numbers
for i = 1:length(vals)
    if vals(i) ~= i
        Mask(Mask == vals(i)) = i;
    end
end
