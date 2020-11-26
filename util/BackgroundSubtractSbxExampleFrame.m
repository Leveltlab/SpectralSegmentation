
clearvars info
clearvars -global info
[fileNameIn, filePathIn] = uigetfile('.sbx');
fileIn = [filePathIn, fileNameIn(1:end-4)];
%% 
fRadius = 60;
method = 'Circular Average'; % possiblities 'Gaussian' | 'Circular Average'
shift = 35 .* 2^8; % Shift background prevents underexposure
padTaker = 15; % How deep into the img edge to take data to average for padding values
smoothDim = 1;
smoothSe = 1.75;

switch method
    case 'Circular Average'
        filt = fspecial('disk', fRadius);
        
    case 'Gaussian Average'
        sigma = fRadius/3;
        filt = fspecial('gaussian', fRadius*2+1, sigma);
        
    otherwise
        warning('no valid method')
end

% [subplot margin top&side],[figure bottomspace,topspace],[leftspace,rightspace]
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.01], [0.04 0.04], [0.05 0.05]);


% Example
frametoRead = 100;
img = sbxread(fileIn, frametoRead, 1);
imgDim = size(img);

% Replicate the outer values via padding, so the filter doesn't use a bunch
% of 0s when calculating the background values at the edges
padTaker = 10; % how much pixels to average for padding?

clearvars imgPadded2 imgPadded
% imgPadded = img;
imgPadded = zeros(imgDim+fRadius*2);

% Fill above and besides the image, with padding that has the averaged value the outer edges of img
imgPadded(1:fRadius, fRadius+1:end-fRadius)         = repmat(mean(img(1:padTaker, :)), [fRadius, 1]); % top
imgPadded(end-fRadius+1:end, fRadius+1:end-fRadius) = repmat(mean(img(end-padTaker:end, :)), [fRadius, 1]); % bottom
imgPadded(fRadius+1:end-fRadius, 1:fRadius)         = repmat(mean(img(:,1:padTaker), 2), [1, fRadius]); % left
imgPadded(fRadius+1:end-fRadius, end-fRadius+1:end) = repmat(mean(img(:,end-padTaker:end), 2), [1, fRadius]); % right
% Fill in the corners
imgPadded(1:fRadius, 1:fRadius)                 = mean(imgPadded(fRadius+1:fRadius+padTaker,1)); % top left corner
imgPadded(1:fRadius, end-fRadius+1:end)         = mean(imgPadded(fRadius+1:fRadius+padTaker,end)); % top right corner
imgPadded(end-fRadius+1:end, 1:fRadius)         = mean(imgPadded(end-fRadius-padTaker:end-fRadius, 1)); % bottom left corner
imgPadded(end-fRadius+1:end, end-fRadius+1:end) = mean(imgPadded(end-fRadius-padTaker:end-fRadius, end)); % bottom right corner
% Put in the actual image
imgPadded(fRadius+1:end-fRadius, fRadius+1:end-fRadius) = img; 
    
% imgPadded = padarray(img, [fRadius fRadius], 'replicate');

if false
    background = conv2(imgPadded, filt, 'same');
    % Remove padding of the img
    background = background(fRadius+1:end-fRadius, fRadius+1:end-fRadius);
else
    background = medfilt2(imgPadded, [fRadius fRadius]);
    background = double(background(fRadius+1:end-fRadius, fRadius+1:end-fRadius));
end

imgCorrected = uint16(double(img) - (background - shift));

% Banding diminishing smooth
if smoothDim == 1
    % Simple smoothing to reduce CMOS sensor vertical banding noise
    imgCorrected = smoothG(imgCorrected', smoothSe)';
elseif smoothDim == 2
    % Simple smoothing to reduce CMOS sensor horizontal banding noise
    imgCorrected = smoothG(imgCorrected, smoothSe);
end


figure('Position',[270 383 1610 595])
subplot(1,3,1)
imagesc(imgPadded); colorbar
% caxis(prctile(imgPadded(:), [60 95]));
subplot(1,3,2)
imagesc(background); colorbar
% caxis(prctile(background(:), [20 90]));
title(sprintf('%s, size %d', method, fRadius*2))
subplot(1,3,3)
imagesc(imgCorrected); colorbar



%% Fix vertical banding
background2 = movmean(imgCorrected, 10, 2);
background2 = double(imgCorrected) + 10 - background2;
background3 = mean(background2);
imgCorrected2 = imgCorrected - uint16(background3);
figure
subplot(2,2,1)
imagesc(background2)
subplot(2,2,2)
imagesc(imgCorrected2)
caxis(prctile(imgCorrected2(:), [45 55]));
subplot(2,1,2)
plot(background3)
xlim([0 length(background3)])
