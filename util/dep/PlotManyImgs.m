function [hAxes, hImgs, xidx, yidx] = PlotManyImgs(imgs, varargin)
% [hAxes, hImgs] = PlotManyImgs(imgs, pos, distr, buffer, titles, cut, selected)
% [hAxes, hImgs] = PlotManyImgs(imgs)
% 
% Plot many images in seperate 
% Input: 
%   - imgs (cell array with 2D or 3D doubles inside): the images to plot.
%   - pos([1x4 double]): Vector saying position limits inside figure:
%               [from left, from bottom, width, height].
%   - distr([1 x 2 double] or string): How to distribute the subplots:
%               [n columns, n rows], or let function figure it out as:
%               'horizontal', 'vertical', 'rect' (default).
%   - buffer([1x4 double]): Space between individual subplots:
%               [top, left, right, bottom].
%   - titles (cell array with strings): titles for the plots.
%   - cut([1x4 double]): Vector saying how to crop every img in px:
%               [vertical top, to vertical down, horizontal left, to horizontal right]
%   - selected ([n x 1 double]): which of the 3rd dimension in each cell of
%                                imgs to select: default: 1 (or 3 for RGB)
% 
% 
% Output:
%   - hAxes (array of handles to the axes)
%   - hImgs (array of handles to the images)
%   - xidx (the x position index of the axes)
%   - yidx (the y position index of the axes)% 
% 
% Leander de Kraker
% 2022-10-7, updated 2023-11-22
% 

nimgs = numel(imgs);

if (nargin >= 2) && ~isempty(varargin{1}) % outer position inside figure
    pos = varargin{1};
else
    pos = [0.1000 0.1100 0.7750 0.8150];
end
if (nargin >= 3) && ~isempty(varargin{2}) % Subplot distribution
    distr = varargin{2};
else
    distr = 'rect';
end
if (nargin >= 4) && ~isempty(varargin{3}) % space between individual subplots
    buffer = varargin{3};
else
    buffer = [0.1 0.1 0 0];
end
if (nargin >= 5) && ~isempty(varargin{4}) % titles
    titles = varargin{4};
    doTitles = true;
else
    doTitles = false;
end
if (nargin >= 6) && ~isempty(varargin{5}) % cropping/ zooming the images
    cut = varargin{5}; 
else
    cut = [];
end
if (nargin >= 7) && ~isempty(varargin{6}) % 3rd dimensions of image matrix to select
    selected = varargin{6};
else
    if size(imgs{1}, 3) == 3
        selected = [1 2 3];
    else
        selected = 1;
    end
end


if ischar(distr)
    if strcmpi(distr, 'rect')
        distr = [ceil(sqrt(nimgs)),  round(sqrt(nimgs))];
    elseif strcmpi(distr, 'horizontal')
        distr = [1, nimgs];
    elseif strcmpi(distr, 'vertical')
        distr = [nimgs, 1];
    else
        warning('non valid distr input')
    end
end

% CALCULATE SUBPLOT POSITIONS
[xidx, yidx] = find(ones(distr));
yidx = distr(2)+1-yidx;
positionx = linspace(pos(1)+buffer(2)/distr(1), pos(1)+pos(3).*(1+buffer(2)), distr(1)+1);
positiony = linspace(pos(2)+buffer(4)/distr(2), pos(2)+pos(4).*(1+buffer(4)), distr(2)+1);

positionh = [positiony(2)-positiony(1)] .* (1-(buffer(1)+buffer(4)));
positionw = [positionx(2)-positionx(1)] .* (1-(buffer(3)+buffer(2)));

positionx = positionx(xidx);
positiony = positiony(yidx);

positionx = positionx(1:nimgs);
positiony = positiony(1:nimgs);
xidx = xidx(1:nimgs);
yidx = yidx(1:nimgs);

positions = [positionx', positiony', repmat([positionw, positionh], [nimgs, 1])];



hAxes = gobjects(nimgs, 1);
hImgs = gobjects(nimgs, 1);
for i = 1:nimgs
    if ~isempty(imgs{i})
        if ~isempty(cut)
            img = imgs{i}(cut(1):cut(2), cut(3):cut(4), selected);
        else
            img = imgs{i}(:,:,selected);
        end
    end
    
    hAxes(i) = subplot('Position', positions(i,:));
    
    if ~isempty(imgs{i})
        hImgs(i) = imagesc(img);
    end
    if doTitles
        title(titles{i})
    end
end

    
    