function valid = SparCheckValid(spar)
% Checks if spar (spectral parameters/ search parameters) for roi detection
% seem valid
% 
% Leander de Kraker
% 2026-3-16


flds = fieldnames(spar);
aflds = { 'cutOffHzMin', 'cutOffHzMax', 'border', 'areasz',...
    'roundedness', 'voxel', 'cutOffCorr', 'useFluorescenceImg',...
    'doPlot'};
if sum(ismember(aflds, flds)) == 9
    valid = true;
end