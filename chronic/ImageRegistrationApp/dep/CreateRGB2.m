function [RGB] = CreateRGB2(data, colors, varargin)
% Create a colorfull image by overlaying 2D images from a cell array using
% a requested specific colormap
% 
% RGB = CreateRGB2(data, colors, normalization, whitebalance, selected)
% 
% input:
%   - data: ([n x 1] cell array), with 2D or 3D doubles.
%   - colors ([n x 3] double): color values to give to corresponding data.
%   Optional input: can be left empty [] or disregarded completely
%   - normalize: true or false. Normalizes each data slice
%   - wbalance : true or false. White balances final RGB image. Can make
%                     the images not correspond to their intended colors.
%   - selected ([n x 1] double): In case each cell contains a 3D matrix
%                         select which image of the 3rd dimension to use.
%  
%  
% output:
%   - RGB: 3D double [height x width x 3] of largest 2D double in data 
% 
% 
% See also CreateRGB, CreateRGB2_mat
% 
% Leander de Kraker
% 2019-3-1
%

% normalize?
if (nargin > 3) && ~isempty(varargin{1})
    norma = varargin{1};
else
    norma = false;
end
if (nargin >= 4) && ~isempty(varargin{2})
    wbalance = varargin{2};
else
    wbalance = false;
end
if (nargin == 5)
    selected = varargin{3};
else
    selected = 1;
end

nslices = length(data);
dims = (cellfun(@size, data,'UniformOutput',false));
dims = cat(1,dims{:});
if size(dims, 2)==3
    for i = 1:nslices
        data{i} = data{i}(:,:,selected);
    end
end

maxdim(1) = max(dims(:,1));
maxdim(2) = max(dims(:,2));
RGB = zeros(maxdim(1), maxdim(2), 3);


for i = 1:nslices
    
    if norma % normalize intensities of the channel
        data{i} = (data{i}-min(data{i}(:)))./ range(data{i}(:));
    else % Set the data to minimum 0
        data{i} = data{i}-min(data{i}(:));
    end
    
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


