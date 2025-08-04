function [linkMat, nLinks] = FullLinkMat(linkMat, nrois)
% Extend the linkMat to include all single ROIs as well
% 
% Leander de Kraker
% 2025-7-25
% 
 
for i = 1:length(nrois)
    % keep roi numbers which are not members of this linkmat column
    allRoisi = 1:nrois(i);
    singlesi = allRoisi(~ismember(allRoisi, linkMat(:,i)));
    linkMat(end+1:end+length(singlesi), i) = singlesi;
end
nLinks = sum(linkMat>0, 2);
