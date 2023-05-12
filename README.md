SpecSeg: Cross spectral power-based segmentation of neurons and neurites in chronic calcium imaging datasets
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7113362.svg)](https://doi.org/10.5281/zenodo.7113362)
[![View SpectralSegmentation on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://nl.mathworks.com/matlabcentral/fileexchange/118135-spectralsegmentation)
===========

This Github contains 'SpecSeg', an open source calcium imaging processing toolbox to detect regions of interest (ROIs) in calcium imaging datasets.
ROI segmentation is based on cross-spectral power over a range of frequencies and has been tested and optimized for a variety of neuronal compartments
(cell bodies, dendrites, axons) and (imaging) techniques (chronic glass windows, GRIN lenses, 1P, 2P). 
The SpecSeg toolbox contains user-friendly graphical interfaces to detect, adjust ROIs and visualize their activity traces.
This pipeline includes code adapted from [normcorre](https://github.com/flatironinstitute/NoRMCorre) ([Pnevmatikakis & Giovannucci 2016](https://doi.org/10.1016/j.jneumeth.2017.07.031)),
 to apply motion correction to sbx files. To be able to execute the motion correction, download the code from 
 [NoRMCorre](https://github.com/flatironinstitute/NoRMCorre) and add the code into the matlab path below SpectralSegmentation.<br /> 
Spike estimation can be done on the fluorescence traces, with code that uses [*MLspike*](https://github.com/MLspike/spikes) ([Deneux et al. 2016](https://doi.org/10.1038/ncomms12190)). 
MLspike requires the [brick](https://github.com/MLspike/brick) toolbox. Add MLspike and brick folders to the Matlab path.<br />
This is the pipeline that is used by the Leveltlab in the Netherlands Institute of Neuroscience (NIN).<br /><br />
The pipeline is easily executable manually per step with the script [SpectralPipeline](https://github.com/Leveltlab/SpectralSegmentation/blob/master/SpectralPipeline.m).<br />
Or via the script that runs the pipeline automatically for as many files as requested [AutomatedPipeline](https://github.com/Leveltlab/SpectralSegmentation/blob/master/AutomatedAnalysis.m), until RoiManagerGUI which requires manual input.<br />


<img src="https://github.com/Leveltlab/SpectralSegmentation/blob/master/docs/Figure_Analysis_Pipeline.png" width="700" align="center">


Readme's
======
[Pipeline manual](https://github.com/Leveltlab/SpectralSegmentation/blob/master/docs/Spectral-Segmentation_Manual.pdf)

[RoiManagerGUI manual](https://github.com/Leveltlab/SpectralSegmentation/blob/master/docs/RoiManagerGUI_Manual_V4.pdf). User interface for data exploration and ROI editing


Why use the spectral segmentation toolbox
===========
- User-friendly.
- Automatic ROI detection, with a user interface for a-priori selection of properties such as ROI size, minimum roundness, minimum signal correlation, etc.
- Powerful ROI editor. Create new ROIs, split existing ROIs, delete ROIs manually or based on properties like size, roundness and spectral power.
Visualize and explore the fluorescence data in many ways.
- Designed for and tested on a range of techniques and neuronal compartments 
(cell bodies, dendrites, axons, 1P, 2P, chronic glass windows, chronic GRIN lenses, cortical and subcortical brain regions)
- To increase speed and decrease memory load fluorescence data is transposed and memory mapped.
- Some scripts increase speed by using parallel processing.
- The cross-spectral power images that result from the spectral analysis show excellent separation of active neural elements from the background. 
The cross-spectral power has the advantage over correlation in the time domain of the fluorescence signal in that it has power over multiple different frequencies. 
Different frequencies contain different signal sources; noise can be visible in some frequencies, while neurons can be visible in others (see figure below).
- Chronic tracking of ROIs over multiple sessions (imaged days to months apart).

<img src="https://github.com/Leveltlab/SpectralSegmentation/blob/master/docs/Figure_spectral-Delier_20190916_002-1-Declutter.png">
This image shows results of the spectral analysis on two-photon calcium imaging in mouse V1.
The average- and maximum fluorescence can show fluorescence that is not related to activity of neurons. For example,
 higher baseline fluorescence of a neuron that is not clearly active (cyan box, with its signal shown on the right), or lower baseline 
 fluorescence in blood vessels (green box). However, the cross-spectral power images shows different active neurons, depending on the frequency.
 The red box indicates a neuron that was not visible in the average- and maximum fluorescence projection, but it is clearly visible in some of the cross-spectral power images.
 Therefore, these cross-spectral power images are very useful for ROI selection.


Installation
======
The SpecSeg pipeline only requires MATLAB to be able run on most operating systems.
1. Download the SpecSeg code from this Github.
2. To be able to use motion correction, download [NoRMCorre](https://github.com/flatironinstitute/NoRMCorre).
3. To be able to do spike estimation on ROI signals, download [*MLspike*](https://github.com/MLspike/spikes) and [brick](https://github.com/MLspike/brick).
4. In MATLAB, add the folders from these repositories to your MATLAB path.

    Transposing datasets (StackTranspose.m) works via a MATLAB mex function.
     The mex function is compiled in MATLAB from the C++ code [Zorder.cpp](https://github.com/Leveltlab/SpectralSegmentation/blob/master/spectral/dep/Zorder.cpp).
     For Windows 64-bit, and Mac 64-bit Zorder is compiled and should work automatically.
    However, for other operating systems it needs to be compiled for your system in MATLAB.<br />
5. (if necessary). Compile Zorder.cpp by setting MATLAB's current directory to SpectralSegmentation/Spectral/Dep and executing: *mex Zorder.cpp*
For this to work a compiler must be installed in MATLAB, get information on that by executing *mex -setup* in MATLAB.
To find a suitable compiler go to: https://www.mathworks.com/support/compilers.
If you get the C++ compiler by installing Visual Studio, make sure to install:
 [*Desktop Development with C++*](https://www.mathworks.com/matlabcentral/answers/443349-how-do-i-install-visual-studio-2017-or-2019-for-use-with-matlab-simulink).


Paper
======
The pipeline is published in Cell Reports Methods:

***Leander de Kraker, Koen Seignette, Premnath Thamizharasu, Bastijn J.G. van den Boom, Ildefonso Ferreira Pica, Ingo Willuhn, Christiaan N. Levelt, Chris van der Togt,
SpecSeg is a versatile toolbox that segments neurons and neurites in chronic calcium imaging datasets based on low-frequency cross-spectral power,
Cell Reports Methods, 2022, 100299, ISSN 2667-2375, https://doi.org/10.1016/j.crmeth.2022.100299.***


Matlab dependencies
======
The following Matlab toolboxes are used in the Spectral Segmentation pipeline.
They are required for some of the functions and scripts in the toolbox.

- Image Processing Toolbox.
- Signal Processing Toolbox.
- Statistics and Machine Learning Toolbox.
- Statistics Toolbox.
- (optional) Parallel Computing Toolbox. 
- (optional) Polyspace Bug Finder.


Troubleshooting/ questions
======
For any questions, troubleshooting or comments, 
please contact Chris van der Togt (c.vandertogt@nin.knaw.nl) or Leander de Kraker (l.de.kraker@nin.knaw.nl).
Or open an [issue](https://github.com/Leveltlab/SpectralSegmentation/issues) in this Git.
	
Developers
======
- Chris v.d. Togt.  **Netherlands Insititue of Neuroscience (NIN)**
- Leander de Kraker.  **Netherlands Insititue of Neuroscience (NIN)**

cite as: Leander de Kraker, & Chris van der Togt. (2022). Leveltlab/SpectralSegmentation: Public Release 1.0.1 (V1.01). Zenodo. https://doi.org/10.5281/zenodo.6993003

License
=======
This software is published under creative common license **CC-BY-NC 4.0**


	
3rd party functions & toolboxes
======
Download the code from these toolboxes to be able to use them in the Spectral Segmentation toolbox.
- Motion correction code [*NoRMCorre*](https://github.com/flatironinstitute/NoRMCorre) is adapted to process .sbx files. Place the NoRMCorre path below the SpectralSegmentation path in Matlab.
- ROI calcium signal signal spike estimation with [*MLspike*](https://github.com/MLspike/spikes).

The following 3rd party functions are used in the spectral segmentation toolbox, they are already present in the Spectral Segmentation toolbox.
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
	
- de Kraker, L., Seignette, K., Thamizharasu, P., Boom, B.J.G., Ferreira Pica, I., Willuhn, I., Levelt, C N. & Togt, C. (2022) 
	SpecSeg is a versatile toolbox that segments neurons and neurites in chronic calcium imaging datasets based on low-frequency cross-spectral power
	doi: [10.1016/j.crmeth.2022.100299](https://doi.org/10.1016/j.crmeth.2022.100299)
	
	
- Kraker L., Seignette, K., Thamizharasu, P., Boom, B.J.G., Pica, I.F., Levelt, C.N. & Togt, C. (2020) 
	Specseg: cross spectral power-based segmentation of neurons and neurites in chronic calcium imaging datasets. Cell Reports Methods, 2(10), 100299
	doi: [10.1101/2020.10.20.345371](https://doi.org/10.1101/2020.10.20.345371)

