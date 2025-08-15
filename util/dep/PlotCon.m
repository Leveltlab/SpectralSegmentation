function [h, cons] = PlotCon(PP, colors, xshift, yshift, idx)
% [h, cons] = PlotCon(PP, colors, xshift, yshift, idx);
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
% - idx ([n x 1] double) which ROIs to plot
% 
% 
% Ouput:
%  - h ([scalar] graphics handle)
%  - cons (struct): with fields x and y
% 
% Leander de Kraker
% 2023-4-3
% 

arguments
    PP struct
    colors
    xshift (1,1) = 0;
    yshift (1,1) = 0;
    idx = 1:PP.Cnt
end
if size(idx, 2)==1 && size(idx, 1)>1
    idx = idx';
end

cons = struct();
cons.x = [];
cons.y = [];
for r = idx
    cons.x = cat(2, cons.x, [PP.Con(r).x, NaN]);
    cons.y = cat(2, cons.y, [PP.Con(r).y, NaN]);
end
h = plot(cons.x-xshift, cons.y-yshift, 'color', colors);
