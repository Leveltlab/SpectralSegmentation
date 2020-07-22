Spectral Segmentation
===========

This Github contains open source calcium imaging processing tools to aquire neural traces from Calcium imaging data.  
It is designed to be able to detect active neural elements (soma/ dendrites/ axons) in calcium imaging data.
This pipeline includes code adapted from normcorre, to apply motion correction to sbx files ().
This is the pipeline that is used by the Leveltlab in the Netherlands Institute of Neuroscience (NIN). 
The Levelt lab uses the normcorre motion correction from Simons foundation.

<img src="https://github.com/Leveltlab/SpectralSegmentation/blob/master/docs/Figure_Analysis_Pipeline.png" width="700" align="center">

Why use the spectral segmentation toolbox
===========
- Easy to use pipeline.
- The cross-spectral power images that result from the spectral analysis show active neural elements very clearly. 
The cross-spectral power has the advantage over the usual correlation of the fluorescence signal in that it has power over multiple frequencies. 
Different frequencies contain different signal sources; noise can be visible in some frequencies, while neurons can be visible in others (see figure below).
- Automatic ROI detection, with intuitive properties to select before-hand like ROI size, minimum roundness, minimum signal correlation.
- Powerful ROI editor. Create new ROIs, split existing ROIs, delete ROIs manually or based on properties like size, roundness and spectral power.
Visualize and explore the fluorescence data in many ways.
- To increase speed and decrease memory load fluorescence data is transposed and memory mapped.
- Some scripts increase speed by using parallel processing.
- Chronic tracking of ROIs over multiple days to months.

<img src="https://github.com/Leveltlab/SpectralSegmentation/blob/master/docs/Figure_spectral-Delier_20190916_002-1-Declutter.png">

<img src="D:\Documents\GitHub\SpectralSegmentation\docs\Fig_4_v6_DataExamplesproperCropped_KS.png" width="350">


Readme's
======
Roi editor GUI (RoiRejecterGUI) manual

Pipeline manual

Paper
======
We're publishing a paper for this pipeline jippie (WIP)
(reference)


Matlab dependencies
======
The following matlab toolboxes are used in the Spectral Segmentation pipeline.
They are required for some of the functions and script in the toolbox.

- Image Processing Toolbox.
- Signal Processing Toolbox.
- Statistics and Machine Learning Toolbox.
- Statistics Toolbox.
- (optional) MATLAB Parallel Server. 
- (optional) Parallel Computing Toolbox. 
- (optional) Polyspace Bug Finder. 


Developers
======
- Chris v.d. Togt.  **Netherlands Insititue of Neuroscience (NIN)**
- Leander de Kraker.  **Netherlands Insititue of Neuroscience (NIN)**


3rd party toolboxes used
=======
- Motion correction code [*normcorre*](https://github.com/flatironinstitute/NoRMCorre) is adapted to process .sbx files. (not in github yet).
- ROI calcium signal signal spike deconvolution from [*CaIman*](https://github.com/flatironinstitute/NoRMCorre). (not in github yet).


3rd party functions
=======
The following functions are used in the spectral segmentation toolbox, they are inside the dependency folders.

- [*subtightplot*](https://www.mathworks.com/matlabcentral/fileexchange/39664-subtightplot): 
	Felipe G. Nievinski (2020)
- [*xcorr2\_fft*](https://www.mathworks.com/matlabcentral/fileexchange/53570-xcorr2\_fft-a-b):
	Alessandro Masullo (2020).
	
Troubleshooting/ questions
========
Github l.de.kraker @ nin.knaw.nl
	
Licence
========

	
	
	
	
	
	