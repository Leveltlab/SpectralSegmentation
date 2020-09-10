%% Process for spectral analysis
% Define rois through a two step process, first using the cross spectral 
% power of adjacent pixels and then constraining predefined rois using 
% cross correlation analysis.

%% sbx data is required. TIfF can be converted to sbx
tif2sbx

%% Motion correct sbx data
Showsbx 
% Rightclick on main image to adjust color scale
% Press the A button for the Alignment/motion correction


%% Transpose the motion corrected sbx file to improve reading speed
% Creates _Trans.dat file
StackTranspose

%% Decimate the Trans.dat file to ~1 Hz, and convert data to double
% Creates _DecTrans.dat file
DecimateTrans

%% Then you need to calculate the cross spectral power for each pixel with 8
%surrounding pixels. Requires DecTrans.dat file to function
% Creates _SPSIG.mat file
spectral

%% Create calcium fluorescence images (maximum and average projection)
% Projections get saved into SPSIG file
BackgroundImgSbx;

%% Then do segmentation (neurons, or boutons)
% ROIs, ROI properties and used spectral images get saved into SPSIG file
getSpectrois

%% Remove or add ROIS manually, display and compare roi sets
% If save button is pressed, edits are made in the SPSIG file
RoiManagerGUI

%% Retrieve the signals associated with these rois from the aligned data
% Signals get added to the _SPSIG file
retrievesignals

%% deconvolve signals
% MLspike spike estimation signals get added to the _SPSIG file
DeconvolveSignals

%% Match ROIs of multiple recordings together
ChronicMatching

