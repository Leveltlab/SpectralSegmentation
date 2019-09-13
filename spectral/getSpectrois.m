%
% Get ROIs from cross spectral components
% Edit the ROIs with the RoiRejecterGUI
% 
%

 %% Load and Process Spectral Images:  load('SPic.mat')
[fn, pn] = uigetfile('*_SPSIG.mat');
filename = [pn fn];
fprintf('\nloading...')
load(filename)
Transfile = [filename(1:end-9) 'Trans.dat'];
[sbxt, ~, freq] = transmemap(Transfile);
fprintf('\nloaded\n')

% Activate much tighter subplots
% [subplot margin top&side],[figure bottomspace,topspace],[leftspace,rightspace]
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.01], [0.06 0.05], [0.06 0.05]);

% obtain average spectral density for each image: decays exponentially
imgStack = log(SPic(:,:,2:end));
Sax(1) = []; %first spectral component is the average power over al components

vals = sort(unique(imgStack(:))); % remove -inf values
if vals(1)==-inf
    imgStack(imgStack==-inf) = vals(2);
end
imgStackT = permute(imgStack,[2 1 3]); % transpose the SPic variable so it's same as BImg


Slow = find(Sax <= 0.2);
Spect = imgStack(:,:,Slow);

%Fast = find(Sax > 0.5 & Sax <= 1.0);
%Sfast = mean(Imgstack(:,:,Fast),3);
%Spect = cat(3, Spect, Sfigure, fast);
Spect = permute(Spect, [2 1 3]); %transposes the images,
Spect = setminlevel(Spect); %set minimum level to zero

BImg = max(Spect, [], 3);
figure; subplot(1,1,1)
imagesc(BImg), caxis(prctile(BImg(:), [0.01 99.9])), hold on, colorbar


%% check how the spectral looks (press any key to advance, ctrl+c to stop)
nSpec = size(imgStackT,3);

% figure
% subplot(1,1,1)
% for i = 1:nSpec
%     imagesc(squeeze(imgStack(:,:,i)))
%     title(sprintf('i = %d. frequency %.2fHz', i, Sax(i)))
%     colorbar
%     pause
% end

% A spectral image for all the different frequencies
colors = flip(jet(nSpec));
norma = false; % normalize each frequency?
figure
subplot(1,1,1)
SPicTCell = squeeze(num2cell(imgStackT, [1 2]));
RGBspectral = CreateRGB2(SPicTCell, 1:nSpec, colors, norma);
imagesc(RGBspectral)
title('The intensity of the spectral density at all different frequencies')
% Set the correct colorbar for this image
colormap(colors)
h = colorbar;
hTicks = linspace(Sax(1),Sax(end),length(h.Ticks));
h.Ticks = h.Ticks; % prevent ticks from changing
h.TickLabels = num2str(hTicks', '%2.2f');
ylabel(h, 'spectral frequency', 'FontSize',12)


%% Plot the maximum and average projection for comparison with spectral

if exist('BImgA','var') && ~exist('BImgMax','var')
    BImgMax = BImgA;
end

% Plot and save picture
hand2 = figure;
colors = cmapL([1     1    1;...
                1     1    0.75;...
                0.75  1    0.5;...
                0.25  0.75 0.5;...
                0     0.25 0.5;...
                0     0    0], 256);
subplot(1,1,1)
imagesc(log1p(BImgMax), prctile(log1p(BImgMax(:)), [0.01 99.8]))
colormap(colors);
title('log scaled maximum fluorescence image')
h = colorbar;
ylabel(h, 'maximum fluorescence', 'FontSize',12)
% Save picture
print(hand2, [fn(1:end-18) 'img_maximumFluorescence.png'],'-dpng', '-r400')

if exist('BImgAverage','var') % Plot average fluorescence
    colormap(colors);
    hand3 = figure;
    subplot(1,1,1)
    imagesc(log1p(BImgAverage), prctile(log1p(BImgAverage(:)), [0.01 99.8]))
    colormap(colors);
    title('log scaled average fluorescence image')
    h = colorbar;
    ylabel(h, 'average fluorescence', 'FontSize',12)
    % Save picture
    print(hand3, [fn(1:end-18) 'img_averageFluorescence.png'],'-dpng', '-r400')
end
            
vals = prctile(BImgMax(:),99.5);  % Cut off the highest 0.5% of values
BImgMax(BImgMax>vals) = vals;

figure;
imagesc(CreateRGB({Spect(:,:,1), Spect(:,:,end), log1p(BImgMax)}, 'r b g'))
title('red: lowest frequency, blue highest low frequency, green: max projection')


 %% ROI selection based on cross spectral density
 %find ROIs with higest local maxima in low frequency component images

PP = [];
PP.Cnt = 0; 
dim = size(Spect); %width, height, z 
Mask = zeros(dim(1:2)); 
SpatialCorr = zeros(dim(1:2));

spar = Spectroiparm(); %reads roi segmentation parameters from file
for i = 1:dim(3)
    tic
    fprintf('iteration %2d/%2d: ', i, dim(3))
    Img = Spect(:,:,i); 
    Img(Mask>0) = 0;
    figure(2), hold off,imagesc(Img), colormap gray, hold on %, caxis([h0 h1])     
    [PP, Mask, SpatialCorr] = roisfromlocalmax(Img, PP, Mask, spar, sbxt, freq, SpatialCorr);   

    str = sprintf('number of ROIs found: %5d. time elapsed = %.2fminutes\n', PP.Cnt, toc/60);
    fprintf(str)
    
    Con = PP.Con;
    for k = 1:PP.Cnt
        plot(Con(k).x, Con(k).y, 'r')
    end
    title(str)
    pause(0.01)
end

%swap x and y values, otherwise positions are not correct in Displayrois
P = PP.P; % plot(x,y), plots at the position (y, x) in an image
PP.P(1,:) = P(2,:);
PP.P(2,:) = P(1,:);

nSpec = size(imgStack,3);
PP.SpecProfile = zeros(nSpec, PP.Cnt); % For each ROI save the spectral density values
Mask3D = repmat(Mask, [1,1,nSpec]);

for i = 1:PP.Cnt
    %also redefine peak maxima based on projected maxima in BImg
    ROIi = Mask == i; 
    npixels = sum(ROIi(:)); % number of pixels for this ROI
    PP.P(3,i) = max(BImg(ROIi)); %max
    PP.P(4,i) = min(BImg(ROIi)); %min
    
    % Retrieve the average spectral profile of the ROI
    specProfile = reshape(imgStackT(Mask3D==i), [npixels, nSpec]);
    PP.SpecProfile(:,i) = mean(specProfile);
    
    %average pixel correlation for each ROI
    PP.P(5,i) = mean(SpatialCorr(Mask==i));
end

% Calculate which frequency had the highest spectral density values for
[~, peakFreq] = max(PP.SpecProfile);
PP.peakFreq = Sax(peakFreq);

%% Display  ROIS
 
%contours
hold on
Con = PP.Con;
for k = 1:PP.Cnt
    vx = Con(k).x;
    vy = Con(k).y;
    plot(vx, vy, 'r')
end


%SpatialCorr
figure, imagesc(SpatialCorr)
hold on
Con = PP.Con;
for k = 1:PP.Cnt
    vx = Con(k).x;
    vy = Con(k).y;
    plot(vx, vy, 'r')
end

% mask
cm = flipud(parula(max(Mask(:))));
cm(1,:) = [0.5, 0.5, 0.5];
figure, imagesc(Mask), colormap(cm)
for i = 1:PP.Cnt
    [iy,ix] = find(Mask == i);
    mx = median(ix);
    my = median(iy);
    text(mx, my, num2str(i), 'color', 'w')
end


%% Reject and add ROIs, & look at signal with UI!
RoiRejecterGUI(Mask, PP, SPic, Transfile, Sax)


%% Calculate power spectral density function for each roi and then save
ln = length(Sax);
Ps = zeros(ln,PP.Cnt);
Mt = Mask';
SPImg = reshape(imgStack, [], ln);
for i = 1:PP.Cnt
    ichan = find(Mt(:)== i);
    Ps(:,i) = mean(SPImg(ichan,:)); %and averge this selection 
end

% Save rois
% save(filename, 'PP', 'Mask', 'BImg', 'Ps', 'spar', 'SpatialCorr', '-append')

save(filename, 'PP', 'Mask', 'BImg', 'Ps', 'spar', '-append')
fprintf('saved.\n')


