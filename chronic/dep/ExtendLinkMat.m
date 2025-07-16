function [linkMatAllRois] = ExtendLinkMat(linkMat, nrois)
% Extend linkMatrix and score with the ROIs that didn't find any matches
% 
% Leander de Kraker
% 2025-4-2

linkMatAllRois = linkMat;
nfiles = size(linkMat, 2);
for i = 1:nfiles
    % keep roi numbers which are not members of this linkmat column
    allRoisi = 1:nrois(i);
    singlesi = allRoisi(~ismember(allRoisi, linkMat(:,i)));
    linkMatAllRois(end+1:end+length(singlesi), i) = singlesi;
end