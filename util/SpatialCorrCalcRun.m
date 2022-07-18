% Calculate spatial corr image and inner correlation variance:
% correlation of signal in existing ROIs with the signal from the spectral
% maximum / origin point of the ROI
% 
% To calculate these variables the transposed motioncorrected raw data is
% necessary, together with the PP and Mask variable from the SPSIG.mat file
% 
% 
% Leander de Kraker
% 2020-1-20
% 

% CALCULATE variance of correlation for PP update?? % % %
calcRvar = true;

% Load the mask and contour data, and the old spatialcorr
[fnSpsig, pnSpsig] = uigetfile('*SPSIG.mat','load mask data (SPSIG file)');
nameSpsig = [pnSpsig, fnSpsig];
load(nameSpsig, 'PP', 'Mask', 'SpatialCorr')

% Load transposed sbx data
[fnSbx, pnSbx] = uigetfile('*DECTRANS.dat', 'load raw data (sbx file)');
nameSbx = [pnSbx, fnSbx];
[sbxt, dim, freq] = transmemap(nameSbx);

% % enlarge the ROI masks for extra viewing pleasure!
% [~,~,maskExtra] = BufferMask(Mask, 2); 
% Mask = Mask + maskExtra;

% Checking whether or not SpatialCorr already exists
if exist('SpatialCorr', 'var')
    fprintf('spatialcorr already exists\n')
    oldExists = true; % Old spatial corr already exists!
else
    oldExists = false;
end

% Calculate the new SpatialCorr
SpatialCorrNew = zeros(size(Mask));
rvar = zeros(PP.Cnt, 1);
figure
imagesc(SpatialCorrNew);
colormap(jet)
tic
for i = 1:PP.Cnt
    [Corri, idx, rvari] = SpatialCorrCalcFun(sbxt, freq, Mask, i, PP.P([1 2], i), calcRvar);
    SpatialCorrNew(idx) = Corri(idx); % 1D indexing in 2D matrix
    rvar(i) = rvari;
    if mod(i, 10)==1
        imagesc(SpatialCorrNew)
        colorbar
        hold on 
        for j = 1:i
            plot(PP.Con(j).x, PP.Con(j).y, 'color', [1 0.2 1])
        end
        hold off
        title(sprintf('done with ROI %d/%d: elapsed time: %.2fminutes\n', i, PP.Cnt, toc/60))
        pause(0.01)
    end
end
toc

% Compare with old existing spatialCorr
if oldExists
    figure
    imagesc(CreateRGB({SpatialCorr, SpatialCorrNew}, 'r gb'))
    hold on
    plot(PP.P(1,:), PP.P(2,:), '.g') % ROI center/ origin points
    for i = 1:PP.Cnt % plot ROI contours
        plot(PP.Con(i).x, PP.Con(i).y, 'color', [1 1 1 0.25])
    end
    title('old spatialCorr = r, new spatialCorr = cyan')
end


%% Save the new SpatialCorr
SpatialCorr = SpatialCorrNew;

fprintf('\nSaving the data to the SPSIG file\n')
if calcRvar
    PP.Rvar = rvar;
    save(nameSpsig, 'SpatialCorr', 'PP', '-append')
else
    save(nameSpsig, 'SpatialCorr', '-append')
end
fprintf('\nSaved %s\n', fnSpsig)


