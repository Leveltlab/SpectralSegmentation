%% Pipeline for spectral segmentation analysis
%
% Via this script all major steps of the spectral segmentation pipeline 
% for ROI selection and matching can be executed.
% All of the functions can be called without input, the functions will
% ask which file to analyse.
% Be aware that your image sequence needs to be aligned in advance!

% sbx data is required. TIFF can be converted to sbx
%% Convert (multiple) tiff files that have multiple frames (full movie) to as many sbxfiles
TiffStack2Sbx

%% Convert multiple tiff files that have multiple frames to one sbx file
TiffStack2OneSbx

%% Convert a folder full of tiff files, each tif file consisting of 1 frame
TiffImages2Sbx

%% Convert H5 files to sbx
Hdf52Sbx


%% For 1 photon/ miniscope data background subtraction is necessary
BackgroundSubtractSbx


%% Motion correct sbx data
Showsbx 
% Press the A button for the Alignment/motion correction
% Rightclick on main image to adjust color scale

%% Transpose the motion corrected sbx file to improve reading speed
% Creates _Trans.dat file
StackTranspose

%% Decimate the Trans.dat file to ~1 Hz, and convert data to double class
% Creates _DecTrans.dat file
DecimateTrans

%% Calculate the cross spectral power for each pixel with 8 surrounding 
% pixels. Requires DecTrans.dat file to function
% Creates _SPSIG.mat file
spectral

%% Create calcium fluorescence images (maximum and average projection)
% Projections get saved into SPSIG file (and outputted into the workspace)
FluorescenceImgSbx;

%% Print pngs of fluorescence images & colored spectral (optional)
% Uses images in the SPSIG file created by FLuorescenceImgSbx
PrintBImgs

%% ROI segmentation (neurons, or boutons)
% ROIs, ROI properties and used spectral images get saved into SPSIG file
getSpectrois

%% Remove or add ROIS manually, display and compare roi sets
% If save button is pressed, edits are made in the SPSIG file
RoiManagerGUI

%% Retrieve the signals associated with these rois from the aligned data
% Signals get added to the _SPSIG file
retrievesignals

%% Create a fast Simple Estimate of Activity by Leander (SEAL)
% Seal signal gets added to the _SPSIG file
SealSignals

%% deconvolve signals
% MLspike spike estimation signals get added to the _SPSIG file
DeconvolveSignals

%% Match ROIs of multiple recordings together
% Run this script step by step.
ChronicMatching

