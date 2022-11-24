function [im, h] = Colorbar2D(colors1, colors2, varargin)
% Create a colorbar image based on two colormaps which are combined
% [im, h] = Colorbar2D(colors1, colors2, varargin)
% 
% 
% input:
%   colors1 ([n x 3] double): rgb colormap
%   colors2 ([m x 3] double): different rgb colormap
% 
% ouput:
%   im ([n x m x 3] double): 2D colorbar image
% 
% Leander de Kraker
% 2022-10-20
% 

%%
ncolors1 = size(colors1, 1);
ncolors2 = size(colors2, 1);

im1 = repmat(permute(colors1, [1 3 2]), [1, ncolors2, 1]);
im2 = repmat(permute(colors2, [3 1 2]), [ncolors1, 1, 1]);

% Different mixing options..
im = im1.*im2;
% im = mean(cat(4, im1, im2), 4);
% im = max(cat(4, im1, im2), [], 4);

