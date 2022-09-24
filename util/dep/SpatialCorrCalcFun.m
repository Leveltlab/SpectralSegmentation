function [cor, roiidx, r2] = SpatialCorrCalcFun(sbxt, freq, mask, roiNr, originPos, varargin)
% [cor, roiidx, rvar] = SpatialCorrCalcFun(sbxt, freq, mask, roiNr, originPos, true)
% Calculate the spatial corr of a ROI: correlation of the signal of
% individual pixels with the signal of the median of 9 central ROI pixels
% 
% Can also calculate the square of the correlation!
% 
% inputs: 
%   sbxt: transposed memorymapped sbxt from transmemap(transfile), 
%   freq: 1x1 double, acuisition frequency of 2P data
%   mask: 2D double [height x width], the positions/ id's of ROIs
%   roiNr: which ROI to calculate the spatialcorr etc from
%   originPos: 1D vector [2 x 1] or 2D double [2 x nROIs] which x & y 
%           coordinates are the 'origin position'
%   calcRvar (optional): boolean: true = return rvar
%                                 false = don't calculate rvar (default)
% 
% outputs:
%   cor: 2D double, contains the correlation values, same size as mask
%   roiidx: the 1D indexes of where the ROI is located
%   rvar: variance of correlation values, after median filtering
% 
%
% Leander de Kraker
% 2020-1-20
% 

% Get whether to calculate rvar or not
if nargin == 6
    calcR2 = varargin{1};
else
    calcR2 = false;
end

% Check if we got the actual coordinates for this ROI, or all coordinates
if size(originPos, 1)>1 && size(originPos, 2)>1
    % Check if we need to transpose the coordinates...
    if size(originPos, 2)==2 && size(originPos, 1)>2
        originPos = originPos';
    end
    originPos = originPos([1 2], roiNr); 
end

% And round the coordinates to integers
originPos = round(originPos); 


maskT = mask';
roiidx = find(maskT==roiNr); % The 2D indexes of the ROI, from transposed mask
sigs = sbxt.Data.y(:,roiidx); % signal of all the pixels in the ROI

% subsample the signal to ~1 sample / sec if the frequency is higher
freq = floor(freq);
if freq > 1.5
    sigSubsampled = zeros(ceil(size(sigs,1)/freq), size(sigs,2));      
    for q = 1:size(sigs,2)
        sigSubsampled(:,q) = decimate(double(sigs(:,q)), freq);
    end
    sigs = sigSubsampled;
    clearvars sigSubsampled
end

% Get the origin point of the ROI with 8 surrounding pixels: seedpixels
seed = zeros(size(maskT), 'logical');
seedy = originPos(1)-1:originPos(1)+1;
seedx = originPos(2)-1:originPos(2)+1;
seedy(seedy<1 | seedy>size(maskT,1)) = []; % remove indexes that are outside mask
seedx(seedx<1 | seedx>size(maskT,2)) = [];
seed(seedy, seedx) = true;
seedidx = find(seed);
[~, seedidx] = intersect(roiidx, seedidx); % exlude pixels outside ROI

% take median of signal seedpixels to reduce noise
seedSig = median(sigs(:,seedidx), 2); 

% Fill the spatial correlation image with the correlation values
cor = zeros(size(maskT));
cor(roiidx) = corr(seedSig, sigs, 'rows', 'pairwise');


% Calculate the mean R squared
if calcR2
    thr = 0.5 * max(abs(cor(:))); 
    if thr < 0.5
        %smooth with median filter, since with lower cutoff more variable
        %correlations are included leading to noisy contours
        m = floor((1-thr)*3)* 2 + 1; % let median filter depend on threshold height
        corFiltered = medfilt2(cor, [m m]);
    else
        corFiltered = cor;
    end
    r2 = mean(corFiltered(roiidx).^2);
else
    r2 = nan;
end

% Transpose to the original size
cor = cor';
roiidx = find(mask==roiNr); % The 2D indexes of the ROI

