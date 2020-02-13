%% Process for spectral analysis
% Define rois through a two step process, first using the cross spectral 
% power of adjacent pixels and then constraining predefined rois using 
% cross correlation analysis.
% Note: image sequences should be in sbx format (if not run tif2sbx)

%% First you need to transpose the sbx file with this utility
StackTranspose

%% Then you need to calculate the cross spectral power for each pixel with 8
%surrounding pixels
spectral

%% Create representative images using usual calcium fluorescence (maximum
% and average projection). Gets saved into SPSIG file
BackgroundImgSbx;

%% Then do segmentation (neurons, or boutons)
getSpectrois

%% Remove or add ROIS manually, display and compare roi sets
RoiRejecterGUI()

%% Retrieve the signals associated with these rois from the aligned data
retrievesignals

%% deconvolve signals
DeconvolveSignals
