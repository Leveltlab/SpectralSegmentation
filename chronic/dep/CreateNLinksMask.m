function nLinksMask = CreateNLinksMask(Masks, linkMatAllRois, nLinksAllRois)
% Create images of which ROIs are part of a match with a specific number of
% links/ROIs
% 
% Leander de Kraker
% 2025-4-2

nfiles = length(Masks);
nLinksMask = cell(nfiles, 1);
for i = 1:nfiles
    rois = linkMatAllRois(nLinksAllRois==i,:); % rois with i linked rois
    nLinksMask{i} = zeros(size(Masks{i}));
    for j = 1:nfiles
        roij = rois(:,j);
        roij(roij==0) = [];
        idx = ismember(Masks{j}, roij);
    	nLinksMask{i}(idx) = nLinksMask{i}(idx) + 1;
    end
end