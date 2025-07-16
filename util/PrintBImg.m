function [img, s, m] = PrintBImg(img, varargin)
% rgb = PrintBImg(img, lims, colors, doplot, saveName, {scaleBarX, scaleBarY}, doPerc, doLineShift, pos);
% 
% Apply brightness limits and colormap and Try to shift uneven lines so 
% they align best with the even lines of img. And Save the img in its 
% actual resolution.
% 
% 
% input. All optional except first. All may be left empty except first
%   - img (2D double).
%   - limsPerc ([2 x 1] double). minimum and maximum brightness values.
%                              Outside these values the brightness is cut
%   - colors ([n x 3 double]). colormap to apply to img
%   - doplot (boolean). Plot it in MATLAB for your eyes?
%   - saveName (string). filename for picture (with folder) .
%                      if empty, image will not be saved
%   - scalebar ([2 x 1] cell with doubles idx inside).
%   - doPerc (boolean).: Calculate percentile using the lims input?
%   - doLineShift (boolean): Calculate a line shift correction?
%   - pos ([1 x 4] double): position of the figure on the screen
% 
% Values above this percentile will be max color value. 2p data brightness 
% can be highly skewed % this causes most pixels to be dark. This is not
% nice. So select lower than 100% to% cut off the brightness and thereby 
% making the image brighter/ more contrast.
%
% percentile of values below which the values will be the darkest color.
%
% Leander de Kraker
% 2023-8-17
% 

dims = size(img);

if exist('varargin', 'var') && nargin>1 && ~isempty(varargin{1})
    lims = varargin{1};
else
    lims = [0.1, 99.9];
end
if exist('varargin', 'var') && nargin>2 && ~isempty(varargin{2})
    colors = varargin{2};
else
    colors = cmapL('greenFancy', 256);
end
if exist('varargin', 'var') && nargin>3 && ~isempty(varargin{3})
    doplot = varargin{3};
else
    doplot = true;
end
if exist('varargin', 'var') && nargin>4 && ~isempty(varargin{4})
    saveName = varargin{4};
else
    saveName = [];
end
if exist('varargin', 'var') && nargin>5 && ~isempty(varargin{5}) && varargin{5}~=0
    [scaleBarx, scaleBary]= BurnScaleBarIdx(dims, varargin{5});
    doScalebar = true;
else    
    doScalebar = false;
end
if exist('varargin', 'var') && nargin>6 && ~isempty(varargin{6})
    doPerc = varargin{6};
else
    doPerc = true;
end
if exist('varargin', 'var') && nargin>7 && ~isempty(varargin{7})
    doLineShift = varargin{7};
else
    doLineShift = false;
end
if exist('varargin', 'var') && nargin>8
    pos = varargin{8};
else
    pos = GetCenteredFigPos(dims([2 1])); % figure that does not get saved
end




if doPerc
    lims = prctile(img(:), lims); % color axis limits
end

if doLineShift
    [img, s,  m]      = ShiftLinesImg(img);
    % [img, s(2), m(2)] = ShiftLinesImg(img);
else
    s = [NaN, NaN]; m = [NaN, NaN];
end

img = SetLimits(img, lims, colors);

% Apply scalebar
if doScalebar
    img(scaleBary, scaleBarx, :) = 1;
end

% For our eyes only
if doplot
    figure('Units', 'Pixels', 'InnerPosition', pos)
    subplot('Position', [0 0 1 1])
    image(img);
end

% save image
if ~isempty(saveName)
    imwrite(img, saveName)
end

