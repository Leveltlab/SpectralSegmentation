function img = SetLimits(img, lims, varargin)
% Set the range limits of the img, 
% And if a colormap is inputted, apply that colormap to the image
% 
% img = SetLimits(img, lims, varargin)
% 
% input:
%   - img (2D double).
%   - lims ([1 x 2] double)
%   Optional:
%   - colors ([n x 3] double), colormap
% 
% output:
%   - img (2D double)
% 
% Leander de Kraker
% 2022-10-20
% 

dims = size(img);

img = (img-lims(1))./range(lims);
img(img<0) = 0;
img(img>1) = 1;
if exist('varargin', 'var') && nargin ==3
    colors = varargin{1};
    ncolors = size(colors, 1);
    img = ceil(img*(ncolors-1))+1;
    img = colors(img(:),:);
    img = reshape(img, [dims, 3]);
end