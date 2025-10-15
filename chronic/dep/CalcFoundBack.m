function [nFoundBack, nFoundBackPerc, nameCol] = CalcFoundBack(linkMatAllRois, nLinksAllRois, nrois, filenamesShort, doPlot)
% Calculate how many matches have a specific number of links
% 
% input: 
%   - linkMatAllRois ([nmatches x nfiles] double)
%   - nLinksAllRois ([nmatches x 1] double)
%   - nrois ([1 x nfiles] double). total number of ROIs for each recording
%   Optional:
%   - filenamesShort (cell array with strings)
%   - doPlot (boolean)
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
%   - nameCol ([nfiles, 1] cell with strings)
% 
% Leander de Kraker
% 2025-4-2
% 
arguments
    linkMatAllRois
    nLinksAllRois
    nrois
    filenamesShort cell = {}
    doPlot logical = false
end
if size(nrois,1)>1
    nrois = nrois';
end


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

if doPlot
    figure('name','number of ROIs matched per recording')
    h = uitable('Data',nFoundBack,'Units', 'Normalized', 'Position', [0, 0.5, 1, 0.4]);
    h.ColumnName = [nameCol 'summed ROIs'];
    h.RowName = [filenamesShort; {'summed ROIs'}];
    annotation('textbox', [0, 0.95, 1, 0], 'string', 'Specific number of ROIs which are in a match with that amount of linked ROIs')
    
    h = uitable('Data',round(nFoundBackPerc,1), 'Units','Normalized', 'Position', [0, 0, 1, 0.4]);
    h.ColumnName = [nameCol, {'mean'}];
    h.RowName = [filenamesShort; {'mean'}];
    annotation('textbox', [0, 0.45, 1, 0], 'string', 'Percentage of ROIs which are in a match with that amount of linked ROIs')
end