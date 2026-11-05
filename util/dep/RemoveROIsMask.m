function Mask = RemoveROIsMask(Mask, badROIs)
% Remove the ROIs indicated by badROIs 
% Mask = RemoveROIs(Mask, badROIs)
% 
% Input:
%   Mask - a 2D ROI Mask.
%   badROIs ([n x 1] double) - numbers indicating which ROIs should be
%                              deleted
% 
% Output:
%   Mask - Updated Mask.
% 
% Leander de Kraker
% 2026-9-9
% 

Mask(ismember(Mask, badROIs)) = 0; % Set all badROIs to 0
vals = unique(Mask(:)); % which ROIs are left in the Mask
vals = vals(2:end); % remove the 0       
% Update the ROI ID numbers
for i = 1:length(vals)
    if vals(i) ~= i
        Mask(Mask == vals(i)) = i;
    end
end
