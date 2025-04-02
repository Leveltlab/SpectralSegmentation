function [nFoundBack, nFoundBackPerc, nameCol] = CalcFoundBack(linkMatAllRois, nLinksAllRois, nrois)
% Calculate how many matches have a specific number of links
% 
% input: 
%   - linkMatAllRois
%   - nLinksAllRois
%   - nrois ([1 x nfiles] double). total number of ROIs for each recording
% 
% output: 
%   - nFoundBack ([nfiles+1 x nfiles+1] double): the number of ROIs of
%          each recording (rows) with a specific number of links (columns). 
%          Right margin contains number of ROIs in each recording.
%          Bottom margin contains specific number of ROIs involved in a
%          match with that number of links. (so same as nLinks histogram
%          but multiplied by number of respective links).
%   - nFoundBackPerc ([nfiles+1 x nfiles+1] double): percentage of ROIs
%       of the rowth recording.
%       Right margin contains percentage of ROIs in this recording compared
%       to total number of ROIs in the dataset.
%       Bottom margin contains 
%       
% Leander de Kraker
% 2025-4-2
% 

nfiles = size(linkMatAllRois, 2);
nFoundBack = zeros(nfiles+1);
nameCol = cell(1,nfiles);
for i = 1:nfiles
    present = nLinksAllRois(linkMatAllRois(:,i)>0);
    % Count how many neurons with 1:nfiles matches recording i has
    nFoundBack(i,1:nfiles) = histcounts(present, 1:(nfiles+1));
    nameCol{i} = sprintf('%d ROIs',i);
end
nameCol{1} = 'single';
nFoundBack(nfiles+1,:) = sum(nFoundBack, 1);
nFoundBack(:,nfiles+1) = sum(nFoundBack, 2);

% Percentage wise confusion matrix
n = double(repmat(nrois', [1, nfiles]));
n(end+1,:) = sum(nrois);
n(:,end+1) = sum(nrois);
nFoundBackPerc = round(nFoundBack ./ n .* 100);