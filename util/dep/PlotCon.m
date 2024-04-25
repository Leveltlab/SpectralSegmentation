function [h, cons] = PlotCon(PP, colors, varargin)
% [h, cons] = PlotCon(PP, colors, xshift, yshift);
% PlotCon(PP, colors)
% 
% Quickly plot the ROI contours
% Is accelerated by concatenating the different ROIs and seperating them by
% a nan so no connecting line is drawn between the ROIs. But it is plotted
% as one graphical object for MATLAB.
% 
% Input:
% - PP (struct): with fields Cnt and Con.x and Con.y
% - colors ([1 x 3 double] or [1 x 4 double] or [scalar color]): e.g.
%           [1 0 0]        or [1 0 0 0.2]    or 'r'
% Optional:
% - xshift: shift to apply to data
% - yshift: shift to apply to data
% 
% 
% Ouput:
%  - h ([scalar] graphics handle)
%  - cons (struct): with fields x and y
% 
% Leander de Kraker
% 2023-4-3
% 

if nargin>2
    xshift = varargin{1};
    yshift = varargin{2};
else
    xshift = 0;
    yshift = 0;
end

cons = struct();
cons.x = [];
cons.y = [];
for r = 1:PP.Cnt
    cons.x = cat(2, cons.x, [PP.Con(r).x, NaN]);
    cons.y = cat(2, cons.y, [PP.Con(r).y, NaN]);
end
h = plot(cons.x-xshift, cons.y-yshift, 'color', colors);
