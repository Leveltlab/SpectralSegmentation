% Experiment the effect of background subtraction parameters on Sbx files.
% 
% 
% Leander de Kraker
% 2020-12
%

clearvars info
clearvars -global info
[fileNameIn, filePathIn] = uigetfile('.sbx');
fileIn = [filePathIn, fileNameIn(1:end-4)];
%% 
fRadius = [2 15]; % Radius of the filter
subsample = 3;  % rescaling the data images
method = 'Donut'; % possiblities 'Gaussian' | 'Circular' | 'Donut'
doMeanFilt = true; % Mean filtering (convolution), or median filtering (ordfilt)? true | false
shift = 35 .* 2^8; % Shift background prevents underexposure
padTaker = 10; % How deep into the img edge to take data to average for padding values
% Which dimension to do gaussian smoothing?
smoothDim = 0; %  0: no smoothing. 1: horizontal smoothing for vertical bands. 2: vertical smoothing
smoothSe = 1.75; % Smoothing gaussian standard deviation (size)
freq = 15; % Sampling frequency of the data
frametoRead = 4900:5150; % Which frames to read
nframetoRead = length(frametoRead);

% sigy = [40, 40, 40]; % Signal positions to extract
% sigx = [230, 241 252]; 
sigy = [220, 220, 220]; % Signal positions to extract
sigx = [121, 135 150]; 
sigSize = 11; % Size of signal extraction patch
npos = length(sigx);

switch method
    case 'Circular'
        filt = fspecial('disk', fRadius);
        
    case 'Gaussian'
        sigma = fRadius/3;
        filt = fspecial('gaussian', fRadius*2+1, sigma);
        
    case 'Donut'
        filtInner = fspecial('disk', fRadius(1));
        filt = fspecial('disk', fRadius(2));
        filtInner = filtInner ./ ...
            (filtInner(fRadius(1)+1, fRadius(1)+1) ./ filt(fRadius(2), fRadius(2)));
        filtSize = size(filt, 1);
        fBuf = fRadius(2)-fRadius(1);
        idx =fBuf+1:filtSize-fBuf;
        filt(idx, idx) = filt(idx,idx)-filtInner;
        filt = filt ./ sum(filt(:));
        fRadiusGap = fRadius(1);
        fRadius = fRadius(2);
        
    otherwise
        warning('no valid method')
end

% [subplot margin top&side],[figure bottomspace,topspace],[leftspace,rightspace]
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.01], [0.04 0.04], [0.05 0.05]);

imgOri = sbxread(fileIn, frametoRead(1), 1);
imgOriDim = size(imgOri);
img = imresize(imgOri, 1/subsample);
imgDim = size(img);

sigSub = zeros(nframetoRead, npos);
sigRaw = zeros(nframetoRead, npos);
backgroundMax = zeros(imgOriDim);
imgCorrectedMax = zeros(imgOriDim);
imgMax = zeros(imgOriDim);
clearvars imgPadded2 imgPadded

% i = 1;
for i = 1:nframetoRead
    % Example
    imgOri = sbxread(fileIn, frametoRead(i), 1);
    img = imresize(imgOri, 1/subsample);
    
    if doMeanFilt
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

        background = conv2(imgPadded, filt, 'same');
        % Remove padding of the img
        background = background(fRadius+1:end-fRadius, fRadius+1:end-fRadius);
        background = imresize(background, imgOriDim);
    else
%         background = double(ordfilt2(imgPadded, 5, filt));
%         % Remove padding of the img
%         background = background(fRadius+1:end-fRadius, fRadius+1:end-fRadius);
%         background = double(medfilt2(imgPadded, filt));
        background = double(ordfilt2(imgOri, 5, filt, 'symmetric'));
        % Remove padding of the img
%         background = imresize(background, imgOriDim);
    end
    
%     background = max(background, double(imgOri));
    backgroundMax = max(backgroundMax, background);
    imgCorrected = uint16(double(imgOri) - (background - shift));

    % Banding diminishing smooth
    if smoothDim == 1
        % Simple smoothing to reduce CMOS sensor vertical banding noise
        imgCorrected = smoothG(imgCorrected', smoothSe)';
    elseif smoothDim == 2
        % Simple smoothing to reduce CMOS sensor horizontal banding noise
        imgCorrected = smoothG(imgCorrected, smoothSe);
    end
    
    imgMax = max(double(imgOri), imgMax);
    imgCorrectedMax = max(double(imgCorrected), imgCorrectedMax);
    
    for p = 1:npos
        sigSub(i, p) = mean(mean(imgCorrected(sigx(p):sigx(p)+sigSize, sigy(p):sigy(p)+sigSize)));
        img = imresize(img, imgOriDim);
        sigRaw(i, p) = mean(mean(img(sigx(p):sigx(p)+sigSize, sigy(p):sigy(p)+sigSize)));
    end
end

colors = lines(npos);
style = {'-', '--'};
figure('Position',[270 67 1610 911])
ha = subplot(3,3,[1 4]);
% imagesc(imgPadded); colorbar
imagesc(CreateRGB({img, imgMax}, 'g r'))
% caxis(prctile(imgPadded(:), [60 95]));
ha(2) = subplot(3,3,[2 5]);
imagesc(CreateRGB({background, backgroundMax}, 'g r'))
% caxis(prctile(background(:), [20 90])); colorbar
title(sprintf('%s, size %d, image subsample %.1f', method, fRadius*2, subsample))
ha(3) = subplot(3,3,[3 6]);
imagesc(CreateRGB({imgCorrected, imgCorrectedMax}, 'g r'))
% imagesc(imgCorrected); colorbar
linkaxes(ha, 'xy')
h = subplot(3,3,[8 9]);
hold on
for i = 1:npos
    subplot(3,3,[3 6])
	rectangle('position', [sigy(i),sigx(i),sigSize,sigSize],...
              'edgecolor', colors(i,:), 'LineWidth', 2)
    subplot(3,3,[8 9]);
    yyaxis right
    plot(frametoRead, sigSub(:,i),'color', colors(i,:), 'LineStyle', style{1})
    ylabel('corrected')
    hold on
    yyaxis left
    ylabel('raw')
    plot(frametoRead, sigRaw(:,i),'color', colors(i,:), 'LineStyle', style{2})
    hold on
end
set(h.YAxis, 'Color', [0 0 0]) % black double y axes
subplot(3,3,7)
imagesc(filt);
title('The filter')
colorbar


%% Fix vertical banding
% background2 = movmean(imgCorrected, 10, 2);
% background2 = double(imgCorrected) + 10 - background2;
% background3 = mean(background2);
% imgCorrected2 = imgCorrected - uint16(background3);
% figure
% subplot(2,2,1)
% imagesc(background2)
% subplot(2,2,2)
% imagesc(imgCorrected2)
% caxis(prctile(imgCorrected2(:), [45 55]));
% subplot(2,1,2)
% plot(background3)
% xlim([0 length(background3)])
