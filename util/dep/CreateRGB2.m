function [RGB] = CreateRGB2(data, colors, varargin)
% Create RGB 2
% RGB = CreateRGB2(data, colors, normalization, whitebalance)
% input:
%   - data: 1D cell array with 2D matrixes
%   - colors: In which color values that correspond to the data matrixes.
%   - optional - 'normalize': true or false. Normalizes each data slice
%   -          - 'wbalance' : true of false. White balances final RGB image
% 
% output:
%   - RGB: 2D double [height x width x 3] of largest 2D double in data 
% 
% Create a colorfull image by overlaying 2D images from a cell array using
% a requested SPECIFIC COLORMAP
% 
% 
% Leander de Kraker
% 2019-3-1

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

dims = (cellfun(@size, data,'UniformOutput',false));
dims = cat(1,dims{:});
maxdim(1) = max(dims(:,1));
maxdim(2) = max(dims(:,2));
RGB = zeros(maxdim(1), maxdim(2), 3);
nslices = length(data);

for i = 1:nslices
    
    if norma
        % normalize intensities of the channel
        data{i} = (data{i}-min(data{i}(:)))./ range(data{i}(:));
    else
        data{i} = data{i}-min(data{i}(:));
    end
    % 
    for c = 1:3
        if colors(i,c)>0
            % RGB update = old colorchannel values + imgage i * colormap  value
            RGB(1:dims(i,1),1:dims(i,2),c) = RGB(1:dims(i,1),1:dims(i,2),c) + ...
                                             data{i} .* colors(i,c);
        end
    end
end


% Normalize colors
minvals = min(RGB(:));
maxvals = max(RGB(:));
% Delete -infs and infs if found
if minvals == -inf
    vals = unique(RGB);
    minvals = vals(2); % second smallest value (not -inf)
    idx = RGB < minvals; % indexes of -inf values in this colorchannel
    RGB(idx) = minvals;
end
if maxvals == inf
    vals = unique(RGB);
    maxvals(i) = vals(end-1);
    idx = RGB > maxvals; % indexes of -inf values in this colorchannel
    RGB(idx) = maxvals;
end

% normalize
RGB = (RGB-minvals) ./ range([minvals, maxvals]);

if wbalance
    % Normalize colors (per color channel)(results in unrequested colors...
    % ...but it's "white balanced")
    maxvals = squeeze(max(squeeze(max(RGB))));
    minvals = squeeze(min(squeeze(min(RGB))));
    for i = 1:3
        % Normalize
        if range([minvals(i), maxvals(i)])>0
            RGB(:,:,i) = (RGB(:,:,i)-minvals(i)) ./ range([minvals(i), maxvals(i)]);
        end
    end
end

end


