function [RGB, cvals] = CreateRGB(data, colors, varargin)
% Create an RGB image (length x width x 3) array from a cell array
% 
% input:
%   - data: cell array with 2D matrixes
%   - colors: In which color channels to put the slices.
%   - optional - 'binary': make the image binary, anything above 0 becomes
%   1, all other values are 0, per dataslice
%
% output: 
%   - RGB: (length x width x 3) array
%   - cvals: the color value of each slice of data
%
% examples: 
% When you want to create an RGB images of only two of four arrays:
% data = 1x2 cell
% colors = 'rg gb'
% 
% imagesc(CreateRGB(data, colors))
% 
% Puts the insides of cell 1 in the red and green color channels, and the
% insides of cell 4 in the green and blue color channels. And plots the
% resulting image
%
% data = 1x4 cell
% colors = 'r rg gb b'
% Puts the four images in red, red&green(yellow), green&blue(cyan) 
% and blue color channels.
%
% Leander de Kraker
% 2018-2-15
% Last update: 2019-7-10 !slices variable is now not an input anymore!
%       
%


if nargin == 3
    if varargin{1} == 'binary'
        for i = 1:length(data)
            data{i}(data{i}>0) = 1;
            data{i}(data{i}<1) = 0;
        end
    end
end

dims = (cellfun(@size, data,'UniformOutput',false));
dims = cat(1,dims{:});
maxdim(1) = max(dims(:,1));
maxdim(2) = max(dims(:,2));
RGB = zeros(maxdim(1), maxdim(2), 3);
slices = 1:length(data);
cvals = zeros(length(slices),3)+0.4;
if ~iscell(colors)
    colors = strsplit(colors);
end

% Checking input quality
if length(colors) ~= length(slices)
    warning('specify as many color channel allocations as data slices')
    return
end

% 
c = {'r', 'g', 'b'};

% Put slices in requested color channels
for i = slices
    for channel = 1:3
        % Put in the colorchannel if requested
        if ismember(c{channel}, colors{i})
            RGB(1:dims(i, 1), 1:dims(i, 2), channel) = ...
                RGB(1:dims(i,1), 1:dims(i,2), channel) + data{i};
            cvals(i,channel) = 1;
        end
    end
end

% Normalize colors
maxvals = squeeze(max(squeeze(max(RGB))));
minvals = squeeze(min(squeeze(min(RGB))));
for i = 1:3
    % Delete -infs and infs if found
    if minvals(i) == -inf
        vals = unique(RGB(:,:,i));
        minvals(i) = vals(2); % second smallest value (not -inf)
        idx = RGB(:,i) < minvals(i); % indexes of -inf values in this colorchannel
        RGB(idx,i) = minvals(i);
    end
    if maxvals(i) == inf
        vals = unique(RGB(:,:,i));
        maxvals(i) = vals(end-1);
        idx = RGB(:,i) > maxvals(i); % indexes of -inf values in this colorchannel
        RGB(idx,i) = maxvals(i);
    end
    
    % Normalize
    RGB(:,:,i) = (RGB(:,:,i)-minvals(i)) ./ range([minvals(i), maxvals(i)]);
end

end

