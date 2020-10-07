SpecSeg: Cross spectral power-based segmentation of neurons and neurites in chronic calcium imaging datasets
===========

This Github contains 'SpecSeg', an open source calcium imaging processing toolbox to detect regions of interest (ROIs) in calcium imaging datasets.
ROI segmentation is based on cross-spectral power over a range of frequencies and has been tested and optimized for a variety of neuronal compartments
(cell bodies, dendrites, axons) and (imaging) techniques (chronic glass windows, GRIN lenses, 1P, 2P). 
The SpecSeg toolbox contains user-friendly graphical interfaces to detect, adjust ROIs and visualize their activity traces.
This pipeline includes code adapted from [normcorre](https://github.com/flatironinstitute/NoRMCorre) ([Pnevmatikakis & Giovannucci 2016](https://doi.org/10.1016/j.jneumeth.2017.07.031)),
 to apply motion correction to sbx files. To be able to execute the motion correction, download the code from 
 [normcorre](https://github.com/flatironinstitute/NoRMCorre) and add the code into the matlab path below SpectralSegmentation.<br /> 
Spike estimation can be done on the fluorescence traces, with code that uses [*MLspike*](https://github.com/MLspike/spikes) ([Deneux et al. 2016](https://doi.org/10.1038/ncomms12190)).<br />
This is the pipeline that is used by the Leveltlab in the Netherlands Institute of Neuroscience (NIN).<br /><br />
The pipeline is easily executable manually per step with the script [SpectralPipeline](https://github.com/Leveltlab/SpectralSegmentation/blob/master/SpectralPipeline.m).<br />
Or via the script that runs the pipeline automatically for as many files as requested [AutomatedPipeline](https://github.com/Leveltlab/SpectralSegmentation/blob/master/AutomatedAnalysis.m), until RoiManagerGUI which requires manual input.<br />


<img src="https://github.com/Leveltlab/SpectralSegmentation/blob/master/docs/Figure_Analysis_Pipeline.png" width="700" align="center">


Why use the spectral segmentation toolbox
===========
- User-friendly.
- The cross-spectral power images that result from the spectral analysis show excellent separation of active neural elements from the background. 
The cross-spectral power has the advantage over correlation in the time domain of the fluorescence signal in that it has power over multiple different frequencies. 
Different frequencies contain different signal sources; noise can be visible in some frequencies, while neurons can be visible in others (see figure below).
- Automatic ROI detection, with a user interface for a-priori selection of properties such as ROI size, minimum roundness, minimum signal correlation, etc.
- Powerful ROI editor. Create new ROIs, split existing ROIs, delete ROIs manually or based on properties like size, roundness and spectral power.
Visualize and explore the fluorescence data in many ways.
- Designed for and tested on a range of techniques and neuronal compartments 
(cell bodies, dendrites, axons, 1P, 2P, chronic glass windows, chronic GRIN lenses, cortical and subcortical brain regions)
- To increase speed and decrease memory load fluorescence data is transposed and memory mapped.
- Some scripts increase speed by using parallel processing.
- Chronic tracking of ROIs over multiple sessions (imaged days to months apart).

<img src="https://github.com/Leveltlab/SpectralSegmentation/blob/master/docs/Figure_spectral-Delier_20190916_002-1-Declutter.png">
This image shows results of the spectral analysis on two-photon calcium imaging in mouse V1.
The average- and maximum fluorescence can show fluorescence that is not related to activity of neurons. For example,
 higher baseline fluorescence of a neuron that is not clearly active (cyan box, with its signal shown on the right), or lower baseline 
 fluorescence in blood vessels (green box). However, the cross-spectral power images shows different active neurons, depending on the frequency.
 The red box indicates a neuron that was not visible in the average- and maximum fluorescence projection, but it is clearly visible in some of the cross-spectral power images.
 Therefore, these cross-spectral power images are very useful for ROI selection.


Readme's
======
[Pipeline manual](https://github.com/Leveltlab/SpectralSegmentation/blob/master/docs/Spectral-Segmentation_Manual.pdf)

[RoiManagerGUI manual](https://github.com/Leveltlab/SpectralSegmentation/blob/master/docs/RoiManagerGUI_Manual_V3.pdf)


Paper
======
We are submitting a manuscript on this pipeline (WIP)
(reference)


Matlab dependencies
======
The following Matlab toolboxes are used in the Spectral Segmentation pipeline.
They are required for some of the functions and scripts in the toolbox.

- Image Processing Toolbox.
- Signal Processing Toolbox.
- Statistics and Machine Learning Toolbox.
- Statistics Toolbox.
- (optional) MATLAB Parallel Server. 
- (optional) Parallel Computing Toolbox. 
- (optional) Polyspace Bug Finder. 


Troubleshooting/ questions
======
For questions, troubleshooting or comments, 
please contact Chris van der Togt (c.vandertogt@nin.knaw.nl) or Leander de Kraker (l.de.kraker@nin.knaw.nl).
	
	
Developers
======
- Chris v.d. Togt.  **Netherlands Insititue of Neuroscience (NIN)**
- Leander de Kraker.  **Netherlands Insititue of Neuroscience (NIN)**


License
=======
This software is published under creative common license **CC-BY-NC 4.0**


	
3rd party functions & toolboxes
======
Download the code from these toolboxes to be able to use them in the Spectral Segmentation toolbox.
- Motion correction code [*NoRMCorre*](https://github.com/flatironinstitute/NoRMCorre) is adapted to process .sbx files. Place the NoRMCorre path below the SpectralSegmentation path in Matlab.
- ROI calcium signal signal spike estimation with [*MLspike*](https://github.com/MLspike/spikes).

The following 3rd party functions are used in the spectral segmentation toolbox, they are inside the dependency folders.
- [*subtightplot*](https://www.mathworks.com/matlabcentral/fileexchange/39664-subtightplot): 
	Felipe G. Nievinski (2020).
- [*xcorr2\_fft*](https://www.mathworks.com/matlabcentral/fileexchange/53570-xcorr2\_fft-a-b):
	Alessandro Masullo (2020).
- [*figtitle*](https://www.mathworks.com/matlabcentral/fileexchange/42667-figtitle):
	Chad Greene (2020).
	
	
References
======
- Deneux, T., Kaszas, A., Szalay, G., Katona, G., Lakner, T., Grinvald, A., Rózsa, B. & Vanzetta, I. (2016) 
	Accurate spike estimation from noisy calcium signals for ultrafast three-dimensional imaging of large neuronal populations in vivo. 
	*Nature Communications* 7, 12190. doi: [10.1038/ncomms12190](https://doi.org/10.1038/ncomms12190) <br />
	
- Eftychios A. Pnevmatikakis and Andrea Giovannucci, (2017)
	NoRMCorre: An online algorithm for piecewise rigid motion correction of calcium imaging data, 
	*Journal of Neuroscience Methods*, vol. 291, pp 83-94, 2017; 
	doi: [/10.1016/j.jneumeth.2017.07.031](https://doi.org/10.1016/j.jneumeth.2017.07.031)
	
	
