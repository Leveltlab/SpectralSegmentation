function [cor, roiidx, rvar] = SpatialCorrCalcFun(sbxt, freq, mask, roiNr, originPos, varargin)
% Calculate the spatial corr of a ROI: correlation of the signal of
% individual pixels with the signal of the median of 9 central ROI pixels
% 
% Can also calculate the variance of the correlation!
% 
% This function can be executed to a SPSIG file via SpatialCorrCalcRun.m
% 
% inputs: 
%   sbxt: sbxt from transemap(transfile), 
%   freq: 1x1 double, acuisition frequency of 2P data
%   mask: 2D double [height x width], the positions/ id's of ROIs
%   roiNr: which ROI to calculate the spatialcorr etc from
%   originPos: 1D vector [2 x 1] or 2D double [2 x nROIs] which x & y 
%           coordinates are the 'origin position'
%   calcRvar (optional): boolean: 0=return rvar, 1=don't calculate rvar
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
    calcRvar = varargin{1};
end

maskT = mask';
roiidx = find(maskT==roiNr); % The 2D indexes of the ROI, from transposed mask
sigs = sbxt.Data.y(:,roiidx); % signal of all the pixels in the ROI

% subsample the signal to ~1 sample / sec
freq = floor(freq);
sigSubsampled = zeros(ceil(size(sigs,1)/freq), size(sigs,2));      
for q = 1:size(sigs,2)
    sigSubsampled(:,q) = decimate(double(sigs(:,q)), freq);
end

% Get the origin point of the ROI with 8 surrounding pixels: seedpixels
seed = zeros(size(maskT), 'logical');
seed(originPos(1)-1:originPos(1)+1, originPos(2)-1:originPos(2)+1) = true;
seedidx = find(seed);
[~, seedidx] = intersect(roiidx, seedidx); % exlude pixels outside ROI

% take median of signal seedpixels to reduce noise
seedSig = median(sigSubsampled(:,seedidx), 2); 

% Fill the spatial correlation image with the correlation values
cor = zeros(size(maskT));
cor(roiidx) = corr(seedSig, sigSubsampled, 'rows', 'pairwise');


% Calculate the mean variance of correlation
if calcRvar
    thr = 0.5 * max(abs(cor(:))); 
    if thr < 0.5
        %smooth with median filter, since with lower cutoff more variable
        %correlations are included leading to noisy contours
        m = floor((1-thr)*3)* 2 + 1; % let median filter depend on threshold height
        corFiltered = medfilt2(cor, [m m]);
    else
        corFiltered = cor;
    end
    rvar = mean(corFiltered(roiidx).^2);
else
    rvar = nan;
end

% Transpose to the original size
cor = cor';
roiidx = find(mask==roiNr); % The 2D indexes of the ROI

