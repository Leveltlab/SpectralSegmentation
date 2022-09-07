function Z = XYgetZ(indices, sbxt, stacksz)
%return array of traces, given XY coordinates
% Chris van der Togt, 2017, 
% Netherlands Institute for Neuroscience 

if ~isobject(sbxt)
    return
end

lngth = length(indices);
indices = indices(:)';
indices = (indices-1) .* stacksz; %each index is a startpoint for the next z pixel trace

ix = repmat((1:stacksz)', 1, lngth) + repmat(indices, stacksz, 1); 
Z = sbxt.Data(ix(:));