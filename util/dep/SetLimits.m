function img = SetLimits(img, lims, varargin)
% Set the range limits of the img, changing the values to between 0 and 1
% And if a colormap is inputted, apply that colormap to the image
% 
% img = SetLimits(img, lims, colors)
% 
% input:
%   - img (2D double).
%   - lims ([1 x 2] double, or empty, to just apply the colormap to img)
%   Optional:
%   - colors ([n x 3] double), colormap
% 
% output:
%   - img (2D double)
% 
% Leander de Kraker
% 2022-10-20
% 

if isempty(lims)
    lims = [min(img(:)), max(img(:))];
end
img = (img-lims(1))./range(lims);
img(img<0) = 0;
img(img>1) = 1;
if exist('varargin', 'var') && nargin ==3
    colors = varargin{1};
    ncolors = size(colors, 1);
    dims = size(img);
    img = ceil(img*(ncolors-1))+1;
    img = colors(img(:),:);
    img = reshape(img, [dims, 3]);
end