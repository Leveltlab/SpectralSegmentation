Spectral Segmentation
===========

This Github contains open source calcium imaging processing tools to detect regions of interest (ROIs) in Calcium imaging data and extract their fluorescence traces.  
It is designed to be able to detect active neural elements (soma/ dendrites/ axons) in both one- or two-photon calcium imaging data.
This pipeline includes code adapted from [normcorre](https://github.com/flatironinstitute/NoRMCorre) (reference to Normcorre),
 to apply motion correction to sbx files ().<br />
Spike estimation can be done on the fluorescence traces, with code from [*MLspike*](https://github.com/MLspike/spikes)
that is inside this pipeline.<br />
This is the pipeline that is used by the Leveltlab in the Netherlands Institute of Neuroscience (NIN).

<img src="https://github.com/Leveltlab/SpectralSegmentation/blob/master/docs/Figure_Analysis_Pipeline.png" width="700" align="center">

Why use the spectral segmentation toolbox
===========
- Easy to use pipeline.
- The cross-spectral power images that result from the spectral analysis show active neural elements very clearly. 
The cross-spectral power has the advantage over correlation in the time domain of the fluorescence signal in that it has power over multiple different frequencies. 
Different frequencies contain different signal sources; noise can be visible in some frequencies, while neurons can be visible in others (see figure below).
- Automatic ROI detection, with intuitive properties to select before-hand like ROI size, minimum roundness, minimum signal correlation.
- Powerful ROI editor. Create new ROIs, split existing ROIs, delete ROIs manually or based on properties like size, roundness and spectral power.
Visualize and explore the fluorescence data in many ways.
- To increase speed and decrease memory load fluorescence data is transposed and memory mapped.
- Some scripts increase speed by using parallel processing.
- Chronic tracking of ROIs over multiple days to months.

<img src="https://github.com/Leveltlab/SpectralSegmentation/blob/master/docs/Figure_spectral-Delier_20190916_002-1-Declutter.png">
This image shows results of the spectral analysis on two-photon microscopy in a mouse V1.
The average- and maximum fluorescence fluorescence can show fluorescence that is not related to activity of neurons. For example,
 higher baseline fluorescence of a neuron that does not show activity (cyan box, with its signal shown on the right), or lower baseline 
 fluorescence in blood vessels (green box). However, the cross-spectral power shows different active neurons, depending on the frequency.
 The red box shows a neuron that was not visible in the average- and maximum fluorescence projection, but it is clearly visible in some of the cross-spectral power images.
 Therefore, these cross-spectral power images are very usefull for ROI selection.

Readme's
======
Pipeline manual

Roi editor GUI [RoiRejecterGUI manual](https://github.com/Leveltlab/SpectralSegmentation/blob/master/docs/RoiRejecterGUI_Manual_V3.pdf)

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



3rd party functions & toolboxes
======
- Motion correction code [*normcorre*](https://github.com/flatironinstitute/NoRMCorre) is adapted to process .sbx files. (not in github yet).
- ROI calcium signal signal spike estimation with [*MLspike*](https://github.com/MLspike/spikes). (not in github yet).

The following functions are used in the spectral segmentation toolbox, they are inside the dependency folders.
- [*subtightplot*](https://www.mathworks.com/matlabcentral/fileexchange/39664-subtightplot): 
	Felipe G. Nievinski (2020)
- [*xcorr2\_fft*](https://www.mathworks.com/matlabcentral/fileexchange/53570-xcorr2\_fft-a-b):
	Alessandro Masullo (2020).
	
	
Developers
======
- Chris v.d. Togt.  **Netherlands Insititue of Neuroscience (NIN)**
- Leander de Kraker.  **Netherlands Insititue of Neuroscience (NIN)**


Troubleshooting/ questions
======
Github issues tracker?
Chris email?
	
Licence
======

	
	
	
	
	
	