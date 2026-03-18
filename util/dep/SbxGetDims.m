function [dwidth, dheight] = SbxGetDims(info)
% Get the width and height of an sbx file 
% 
if isfield(info, 'simon') || (isfield(info, 'scanbox_version') &&  info.scanbox_version == 2.5 )
    dwidth = info.Shape(info.Perm(2));
    dheight = info.Shape(info.Perm(1));
else
    dwidth = info.sz(2);
    dheight = info.sz(1);
end
