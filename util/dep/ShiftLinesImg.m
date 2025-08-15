function [img, bestS, bestMethod, p] = ShiftLinesImg(img, varargin)
% [img, bestS, bestMethod, p] = ShiftLinesImg(img, varargin)
% Some bi-directionally scanned 2P images are misaligned between the even
% and uneven rows.
% This script tries to estimate the best way of re-aligning the rows and
% applies this to the image
% 
% Input: 
%   - img (2D double)
%   OPTIONAL  may be left empty []
%   - shift to apply (scalar double).
%   - method to use (scalar double): 1 = shifting, 2= via resizing
%   - plotting (boolean): to plot or not to plot (default = false)
%   - trans (string): to correct vertical lines, or horizontal lines
%                   'horizontal' (default) | 'vertical'
% 
% Leander de Kraker
% 2023-8-6
% 
%%

plotting = false;
bestS = [];
bestMethod = [];
trans = false;

if exist('varargin', 'var')
    if nargin>1 && ~isempty(varargin{1})
        bestS = varargin{1};
    end
    if nargin>2 && ~isempty(varargin{2})
        bestMethod = varargin{2};
    else
        if isempty(bestS) % Try all methods for registration
            bestMethod = [1 2]; 
        else % Assume shift is done with simple shifting
            bestMethod = 1;
        end
    end
    if nargin==4 && ~isempty(varargin{3})
        plotting = varargin{3};
    end
    if nargin==5 && ~isempty(varargin{4})
        trans = varargin{4};
        if strcmpi(trans, 'horizontal')
            trans = false;
        else
            trans = true;
        end
    end
end

if trans
    img = img';
end

dims = size(img);
% If the image is a 3D double, just use the first image...
if length(dims)==3
    img = img(:,:,1);
end

if isempty(bestS) % Do check of how best to shift the lines
    buf = min(30, floor(dims(2)/2)); % buffer to ignore
    imgTest = img(:,buf:end-buf);
    dimsTest = size(imgTest);
    s = -5:5;
    b = max(s);
    bt = b*2;
    n = length(s);
    p = zeros(n, 2);
    
    imgRest = imgTest(2:2:end, bt:end-bt);
    dimsRest = size(imgRest);
    dimsTake = [floor(dimsTest(1)./2)*2, dimsTest(2)]; % for check
    dimsScheck = [floor(dimsTest(1)/2), dimsTest(2)];
    
    % shifting simply
    if ismember(1, bestMethod)
        for i = 1:n
            imgS = imgTest(1:2:dimsTake(1), s(i)+bt:end-bt+s(i));
            p(i,1) = corr(imgRest(:), imgS(:));
            
            if plotting
                imgRecon = zeros(dimsTake(1), dimsRest(2));
                imgRecon(2:2:end, :) = imgRest;  imgRecon(1:2:end, :) = imgS;
                subplot(4,1,[1 3]); imagesc(imgRecon(1:end/3,:))
                title(sprintf('shifting %d', s(i)))
                subplot(4,1,4); imagesc(s, [1 2], p');
                title(p(i,1))
            end
        end
    end
    
    % Shifting via resizing
    imgRest = imgTest(2:2:end, :);
    if ismember(2, bestMethod)
        % Resizing to the right
        for i = 1:b+1
            imgS = imresize(imgTest(1:2:dimsTake(1), 1:end+s(i)), dimsScheck);
            p(i,2) = corr(imgRest(:), imgS(:));
            
            if plotting
                imgRecon = zeros(dimsTake); imgRecon(2:2:end, :) = imgRest; imgRecon(1:2:end, :) = imgS;
                subplot(4,1,[1 3]); imagesc(imgRecon(1:end/3,:))
                title(sprintf('resizing to right %d', -s(i)))
                subplot(4,1,4); imagesc(s, [1 2], p'); colorbar
                title(p(i,2)); pause
            end
        end
        % Resizing to the left
        for i = 1:b
            imgS = imresize(imgTest(1:2:dimsTake(1), 1+i:end), dimsScheck);
            p(i+b+1, 2) = corr(imgRest(:), imgS(:));
            
            if plotting
                imgRecon = zeros(dimsTake); imgRecon(2:2:end, :) = imgRest; imgRecon(1:2:end, :) = imgS;
                subplot(4,1,[1 3]); imagesc(imgRecon(1:end/3,:))
                title(sprintf('resizing to left %d', s(b+i+1)))
                subplot(4,1,4); imagesc(s, [1 2], p'); colorbar
                title(p(i+b+1,2)); pause
            end
        end
    end
    
    [bestp, bestpidx] = max(p(:));
    [bestS, bestMethod] = ind2sub(size(p), bestpidx);
    bestS = s(bestS);
end

% Apply the shifting
if bestS~=0
    if bestMethod == 1 % Using shifting
        if bestS<0
            img(1:2:end, -bestS+1:end) = img(1:2:end, 1:end+bestS);
        else
            img(1:2:end, 1:end-bestS) = img(1:2:end, 1+bestS:end);
        end
    else % Using shifting-resizing
        dimsS = [ceil(dims(1)/2), dims(2)];
        if bestS<0
            img(1:2:end, :) = imresize(img(1:2:end, 1:end+bestS), dimsS);
        else
            img(1:2:end, :) = imresize(img(1:2:end, bestS:end), dimsS);
        end
    end
end

if trans
    img = img';
end
