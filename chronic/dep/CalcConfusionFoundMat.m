function [confusionFoundMat, confusionFoundMatPerc] = CalcConfusionFoundMat(linkMatAllRois, nrois)
% Calculate confusion-like matrix what number of ROIs recording i linked 
% with recording j
% 
% input:
%   - linkMatAllRois
%   - nLinksAllRois
%   - nrois ([1 x nfiles] double). total number of ROIs for each recording
% 
% output:
%   - confusionFoundMat ([nfiles x nfiles] double). number of ROIs linked
%       between recording i and j. Values symmetrical along diagonal.
%   - confusionFoundMatPerc ([nfiles x nfiles] double). percentage of ROIs
%       linked with recording i and j (recording total ROI row wise). 
%       So values are not symmetrical along diagonal.
% 
% Leander de Kraker
% 2025-4-2

nfiles = size(linkMatAllRois, 2);
confusionFoundMat = zeros(nfiles);
for i = 1:nfiles
    for j = 1:nfiles
        confusionFoundMat(i,j) = sum(linkMatAllRois(:,i)>0 & linkMatAllRois(:,j)>0);
    end
end

n = double(repmat(nrois', [1, nfiles]));
confusionFoundMatPerc = round(confusionFoundMat ./ n * 100);
