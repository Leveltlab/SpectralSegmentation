%% process for spectral analysis
% Getting signals by processing the cross spectral power of adjacent pixels

%first you need to transpose the sbx file with this utility
StackTranspose

%then you need to calculate the cross spectral power for each pixel with 8
%surrounding pixels
spectral

% Create representative images using usual calcium fluorescence (maximum
% and average projection). Gets saved into SPSIG file
BackgroundImgSbx

% then do segmentation (neurons, or boutons), display and
% compare roi sets
getSpectrois

% retrieve the signals associated with these rois from the aligned data
%and deconvolve signals
retrievesignals