% Example of RoiVal2Img
% 
% 
% 

%% Koen data example
cd('F:\2pdata\Koen\PLChandelier\Analyses\PreTrainingGray\')
load('Noor_20201030_005_normcorr_SPSIG_Res.mat', 'Res', 'info')
load('Noor_20201030_005_normcorr_SPSIG', 'PP')

Mask = info.Mask;
red = [info.rois(:).red];

%% Just plot and color the red ROIs

figure
RoiVal2Img(Mask, PP, ones(PP.Cnt, 1),[],[],[], red);

%% Plot property with many different values

vals = [info.rois(:).correye];
weights = rand(1, PP.Cnt, 1); % With random weights per ROI
colors = jet(255);
colorway = 'value';
doPlot = true;
doDots = true;
doCons = false;
dotSize = 10;

figure
[~, him, habar, hacon, hacen] = RoiVal2Img(Mask, PP, vals, weights, colors, colorway, doPlot, red, dotSize, doDots, doCons);
title('corr eye')
hold off

%% simple call
RoiVal2Img(Mask, PP, vals);

%% Plot ROIs randomly clustered into n groups
clc
n = 2;
vals = randi(n, PP.Cnt, 1); % either a 1 or 2
weights = [];
colors = cmapL('painbow', n);
colorway = 'index';
doPlot = true;
doCons = true;
dotSize = 0;


figure
RoiVal2Img(Mask, PP, vals, weights, colors, colorway, doPlot, red, dotSize, doCons);
title(sprintf('%d random groups', n))
hold off


%% Prem data example
cd('C:\Users\Leander\Documents\Baan\PICLe\good\14.90\Timor\The12\')
load('Timor_20170821_003_split3_normcorr_SPSIG', 'PP', 'Mask')

%%
vals = PP.P(4,:);
red = [];
weights = rand(1, PP.Cnt, 1); % With random weights per ROI
colors = jet(255);
colorway = 'value';
doPlot = true;
doDots = true;
doCons = false;
dotSize = 10;

figure
[~, him, habar, hacon, hacen] = RoiVal2Img(Mask, PP, vals, weights, colors, colorway, doPlot, red, dotSize, doDots, doCons);
title('corr eye')
hold off

%% Big fake example

n = 3;
ntot = n^2;
dims = size(Mask);
mask = zeros(dims*n);
pp = PP;
nrois = PP.Cnt;
counter = 0;
for i = 1:n
    for j = 1:n
        counter = counter + 1;
        maskij = Mask;
        maskij(maskij>0) = maskij(maskij>0) + nrois*counter;
        mask((dims(1)*(i-1)+1):(dims(1)*i), (dims(2)*(j-1)+1):(dims(2)*j)) = maskij;
        imagesc(mask)
        pause
    end
end
pp.Cnt = pp.Cnt * ntot;
pp.Con = repmat(pp.Con, [ntot, 1]);
pp.P = repmat(pp.P, [1, ntot]);
pp.A = repmat(pp.A, [1, ntot]);

%% 
vals = pp.P(4,:);
weights = rand(1, pp.Cnt, 1); % With random weights per ROI
colors = jet(255);
colorway = 'value';
doPlot = true;
doDots = false;
doCons = false;
dotSize = 10; 
RoiVal2Img(mask, pp, vals, weights, colors, colorway, doPlot, red, dotSize, doDots, doCons);
