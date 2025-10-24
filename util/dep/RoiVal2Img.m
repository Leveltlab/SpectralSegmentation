function [rgb, him, habar, hacon, hacen] = RoiVal2Img(Mask, PP, vals, weights, colors, colorway, doPlot, red, dotSize, doCons, doDots)
% Color ROIs based on a property (vals). ROIs can be given a specific color 
% via indexing, or the color can be based on the value of a property.
% ROIs can also be labeled as red rois, which colors their contour/ center
% 
% rgb = RoiVal2Img(Mask, PP, vals, weights, colors, colorway, red, dotSize, doDots, doCons);
% [rgb, him, habar, hacon, hacen]  = RoiVal2Img(Mask, PP, vals);
% rgb = RoiVal2Img(Mask, PP, vals, weights, [], [], [], [], true, false);
% 
% input:
%   Mask. double[w x h]: image of ROI ID
%   PP. struct. with fields: Con.x, Con.y, P (ROI 'center' coordinates)
%                            Cnt, number of ROIs
%   vals. double[nrois x 1]: value that belongs to each individual ROI
% optional. Any of the optional inputs may be kept empty []
%   weights. double[nrois x 1]: weight of the value of each ROI
%                               (alpha/opacity). Between 0 and 1.
%   colors. double[ncolors x 3]
%   colorway. string: 'index' | 'value'. If non given. 'index' if only
%           integers in values, 'value' if there are non-integer values
%             'index' assumes that vals says directly which color of the
%             colormap to pick for everyROI.
%             'value' assumes that vals is a range of values, which still
%             have to be converted to a color 
%   red. Which ROIs are red
%   dotSize double [1x1] or [nrois x 1] double. size of ROI 'center' dots
%   doCons: boolean. Plot the ROI contours? default = true;
%   doDots: boolean. Plot the 'center' Dots? default = true;
%   
% 
% output: A plot in the current axes
%   rgb: the color image
%   him: handle to image
%   habar: handle to colorbar
%   hacon: handle to contour lines
%   hacen: handle to ROI center dots
%
% 
% Leander de Kraker
% 2021-12-6
% 2025-10-15: input parsing arguments, optional plot, 10-100x faster code
%
arguments (Input)
    Mask double
    PP struct
    vals double
    weights double = ones(PP.Cnt, 1)
    colors double = jet
    colorway = 'value'
    doPlot (1,1) logical = true
    red logical = false(PP.Cnt, 1)
    dotSize (1,1) double = 1
    doCons (1,1) logical = true
    doDots (1,1) logical = true
end
arguments (Output)
    rgb double
    him   (1,1) matlab.graphics.Graphics
    habar (1,1) matlab.graphics.Graphics
    hacon
    hacen
end


ncolors = size(colors, 1);
nrois = PP.Cnt;
dims = size(Mask);

% Weight of each ROI (opacity)
if any(weights>1)
    fprintf('weights should not be larger than 1!\n')
    weights(weights>1) = 1;
end
if any(weights<0)
    fprintf('weights should not be smaller than 0!\n')
    weights(weights<0) = 0;
end
if isempty(weights)
    weights =  ones(nrois, 1);
end
% transpose to column vector if needed
if size(weights,1)==1 && size(weights,2)~=1
    weights = weights';
end

if isempty(colors)
    colors = jet;
end
if isempty(colorway)
    colorway = 'value';
end


% Asign color to ROI
if strcmp(colorway, 'index')
    colorsroi = colors(vals, :);
    clims = [1 ncolors+1];
else % 'value'
    % Convert value range to the number of colors in the color range
    clims = [min(vals), max(vals)];
    validx = vals - min(vals); % Set min to 0
    validx = validx / max(validx) * (ncolors-1);% 
    validx = validx + 1;
    validx = round(validx);
    colorsroi = colors(validx, :);
end

colorsroi = colorsroi .* repmat(weights, [1,3]);
colorsroi = [0 0 0; colorsroi];

rgb = zeros(dims(1), dims(2), 3);
for c = 1:3
    % Using the roi ID as index for colormap
    ci = reshape(colorsroi(Mask+1,c), dims);
    rgb(:,:,c) = ci;
end

if doPlot
    % Plot rgb image
    him = imagesc(rgb);
    hold on
    
    % colorbar
    habar = colorbar;
    colormap(colors)
    clim(clims)
    if strcmp(colorway, 'index')
        ticks = habar.Ticks;
        toKeep = floor(ticks)==ticks;
        toKeep(end) = false;
        habar.Ticks = ticks(toKeep);
    end
    
    % Plot contours
    colorsCon = ones(nrois, 3);
    colorsCon(red, 2:3) = 0; % make the red cells red
    if doCons
        hacon = gobjects(nrois, 1);
        conW = ones(nrois, 1) * 1; % linewidth for ordinary ROIs
        conW(red) = 2; % linewidth for red cells
        for i = 1:nrois
            hacon(i) = plot(PP.Con(i).x, PP.Con(i).y,...
                        'color', colorsCon(i,:), 'linewidth', conW(i));
        end
    else
        hacon = gobjects(1);
    end
    
    % Plot ROI 'centers'
    if doDots
        hacen = scatter(PP.P(1,:), PP.P(2,:), dotSize, colorsCon, 'filled');
    else
        hacen = gobjects(1);
    end
else
    him = gobjects(1);
    habar = gobjects(1);
    hacon = gobjects(1);
    hacen = gobjects(1);
end

