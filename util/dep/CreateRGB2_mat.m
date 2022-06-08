function [RGB] = CreateRGB2_mat(data, colors, varargin)
% Create a colorfull image by overlaying 2D images from a 3D matrix data
% using a SPECIFIC COLORMAP colors
% CreateRGB2_mat uses a matrix, while CreateRGB2 uses a cell matrix as
% input.
% 
% RGB = CreateRGB2_mat(data, colors, normalization, whitebalance)
% input:
%   - data: 3D matrix [height x width x nslices]
%   - colors: [nslices x 3] RGB color values for every slice of data
%   - optional - 'normalize': true or false. Normalizes each data slice
%   -          - 'wbalance' : true of false. White balances final RGB image
% 
% output:
%   - RGB: 3D double [height x width x 3]
% 
% 
% See also CreateRGB2, CreateRGB
% 
% Leander de Kraker
% 2020-9-17

% normalize?
if (nargin > 3) && (varargin{1} == true)
    norma = true;
else
    norma = false;
end
if (nargin == 5) && (varargin{2} == true)
    wbalance = true;
else
    wbalance = false;
end

dims = size(data);
RGB = zeros(dims(1), dims(2), 3);
if length(dims)>2
    nslices = dims(3);
else
    nslices = 1;
end


% Remove -inf and replace with lowest actual value
minval = min(data(:));
maxval = max(data(:));
if minval == -inf
    vals = mink(data,2);
    data(data<vals(2)) = vals(2);
    minval = vals(2);
    fprintf('CreateRGB2mat: removing -inf value from data with second lowest\n')
end
if maxval == inf
    vals = maxk(data,2);
    data(data>vals(2)) = vals(2);
    maxval = vals(2);
    fprintf('CreateRGB2mat: replacing inf value from data with second highest\n')
end

% Normalize every slice to amplitude range [0 - 1] if norma is requested
minvals = min(min(data));
maxvals = max(max(data));
if norma
    data = (data - minvals) ./ (maxvals-minvals);
else % Otherwise only lower the lowest value to 0, but keeping range
    data = data - minvals;
end


% Create the RGB image of all the slices
for i = 1:nslices
    for c = 1:3
        if colors(i,c)>0
            % RGB update = old colorchannel values + imgage i * colormap  value
            RGB(:,:,c) = RGB(:,:,c) + data(:,:,i) .* colors(i,c);
        end
    end
end

minval = min(RGB(:));
maxval = max(RGB(:));
% normalize so Matlab can plot the image correctly
RGB = (RGB-minval) ./ range([minval, maxval]);

if wbalance
    % Normalize colors (per color channel)(results in unrequested colors...
    % ...but it's "white balanced")
    maxval = squeeze(max(squeeze(max(RGB))));
    minval = squeeze(min(squeeze(min(RGB))));
    for i = 1:3
        % Normalize
        if range([minval(i), maxval(i)])>0
            RGB(:,:,i) = (RGB(:,:,i)-minval(i)) ./ range([minval(i), maxval(i)]);
        end
    end
end

end


