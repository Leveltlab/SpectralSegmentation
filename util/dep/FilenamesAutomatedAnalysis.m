function fileAllnames = FilenamesAutomatedAnalysis(file, do, nSlices)

fileAll

if nSlices>1
    sliceName = strcat(repmat('Split_', [nSlices, 1]), string([1:nSlices]'));
else
    sliceName = '';
end

if do.LineShift && do.LineShiftValue~=0
    
else
    
end
