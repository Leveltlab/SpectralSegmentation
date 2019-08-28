function [contours, zones, mask2] = BufferMask(Mask,bufSize)
% Creates a contour and a mask around each ROI in the given mask.
% The contour/ masks will be bufSize pixels wide around the original mask. 
% The contours/ masks will not cover ROIs from the input.
%
% input: 
%   - Mask (2D double)
%   - bufSize (digid)
% output: 
%   - contours: a struct with contours for each ROI, and the size 
%       and number of these contours.
%           - C: contour line coordinates
%           - A: number of pixels inside the contours
%	- zones: a struct with a mask for each ROI and its offset
%           - mask: A small mask for each ROI
%           - x: the x coordinates that mask has in the original mask
%           - y: the y coordinates that mask has in the original mask
%           - A: number of pixels inside the contours
%   - mask2: New expanded mask (2D double). Many ROIs should be
%           overlapping, which is not possible in this mask. Only use this
%           when the ROI numbers do not matter (for quick images etc.)
%
% Leander de Kraker
% 2018-2
%

nMasks = max(max(Mask));
dimsM = size(Mask);
Mask = padarray(Mask,[bufSize+1,bufSize+1]);
contours = struct('C', struct('v',[],'x',[],'y',[],'c',[]), 'A', cell(1, nMasks));
zones = struct('mask',cell(1, nMasks),...
    'x',cell(1, nMasks),'y',cell(1, nMasks),'A',cell(1, nMasks));
mask2 = zeros(dimsM);

% Create circular filter
[columns, rows] = meshgrid(1:bufSize*2, 1:bufSize*2);
radius = bufSize;
filt = (rows - bufSize).^2 + (columns - bufSize).^2 <= radius.^2;

for m = 1:nMasks
    [x, y] = find(Mask==m); % Where to expand this mask
    xbeg = min(x) - bufSize - 1;
    xend = max(x) + bufSize + 1;
    ybeg = min(y) - bufSize - 1;
    yend = max(y) + bufSize + 1;
    
    zone = Mask(xbeg:xend, ybeg:yend); % Take small window of dataset
    zone(zone ~= m) = 0; % Only focus on ROI m
    
    % Create mask using circular filter
    zone = filter2(filt, zone);
    zone(zone>0) = m;
    
    % remove parts that fall inside ROIs
    zone(Mask(xbeg:xend, ybeg:yend)>0) = 0; 
    
    
    % Save contours
    % Number of pixels that will be in this contour
    contours(m).A = length(find(zone==m));
    % Get contour of this mask
    C = contourc(zone, [m m]);
    C = getcontourlines(C);
    for i = 1:length(C) % correct x & y positions of this contour
        C(i).x = C(i).x + ybeg - bufSize-2; % (x - y = correction?)
        C(i).y = C(i).y + xbeg - bufSize-2;
    end
    contours(m).C = C;
    

    % Save search zone masks
    % the coordinates of this zone
    zones(m).x = xbeg-bufSize-1:xend-bufSize-1;
    zones(m).y = ybeg-bufSize-1:yend-bufSize-1;
    % If the zone is actually outside unpadded image bounds, cut the zone
    if zones(m).x(1) <= 0 % Cut from x beginning
        cut = 1 -  zones(m).x(1);
        zones(m).x(1:cut) = []; 
        zone(1:cut,:) = [];
    end
    if zones(m).x(end) > dimsM(1) % Cut from x ending
        cut = zones(m).x(end) - dimsM(1);
        zones(m).x(end-cut:end) = []; 
        zone(end-cut:end,:) = [];
    end
    if zones(m).y(1) <= 0 % Cut from y beginning
        cut = 1 - zones(m).y(1);
        zones(m).y(1:cut) = [];
        zone(:,1:cut) = [];
    end
    if zones(m).y(end) > dimsM(2) % Cut from y ending
        cut = zones(m).y(end) - dimsM(2);
        zones(m).y(end-cut:end) = [];
        zone(:,end-cut:end) = [];
        
    end
    zones(m).mask = zone;
    zones(m).A = length(find(zone==m));
    
    
    % Update the mask2 (without overwriting older masks with zeros)
    maskTemp = mask2(zones(m).x, zones(m).y);
    maskTemp(maskTemp==0) = zone(maskTemp==0);
    mask2(zones(m).x, zones(m).y) = maskTemp;
end
end