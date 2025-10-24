function [idx, vals] = ChronicVal2i(linkMat, nrois, linkVal, i)
% map values for ROI matches back to the original recordings
% 
% Input: 
%   - linkMat ([nmatches x nfiles] double): matrix which stores ROI ID
%               which are linked in a row
%   - nrois ([nfiles x 1] double): how many ROIs were present in each file
%   Optional:
%   - linkVal ([nmatches x 1] double): values belonging to every match
%   - i ([ni x 1] double): which recordings to calculate output for.
% 
% output:
%  - idx ([ni x 1] cell array with doubles [nrois x 1] inside): For each 
%                          ROI, in which location of the linkMat it appears
%  - vals ([ni x 1] cell array with doubles [nrois x 1] inside) For each
%                           ROI, what value it had for the given values of
%                           the match
% 
% Leander de Kraker
% 2025-10-16
% 
arguments
    linkMat double
    nrois double
    linkVal double = nan(size(linkMat, 1), 1)
    i double = 1:size(linkMat, 2)
end
if size(linkVal, 2)>1 && size(linkVal, 1)==1
    linkVal = linkVal';
end
linkVal = [NaN; linkVal]; % Add a NaN for when the ROI was mot matched match1 is now in place2
ni = length(i);
nmatches = size(linkMat,1);

idx = cell(ni, 1);
vals = cell(ni, 1);
for ii = i
    idxii = linkMat(:,i(ii)) + 1;
    idxii = accumarray(idxii, (1:nmatches)', [nrois(i(ii))+1 1]);
    idxii(1) = []; % Get rid of the one value which is related to the ROI not being present
    idx{ii} = idxii;
    idxii = idxii+1; % Turn into indexes by turning 0 into 1. ROI 1 is not value 2
    vals{ii} = linkVal(idxii);
end