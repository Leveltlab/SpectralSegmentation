# MATLAB Image registration and GUI

## Method:
</br>
Image registration is performed pairwise, i.e. at each iteration one image is selected as a fixed image, and another is chosen as the moving image, which is shifted (Transformed) to match up onto the fixed image as best as possible.
The current method for registration is built around affine transformations (Translation, Rotation, Scaling, Shearing) and takes part in three steps.

-  Creating a "Pyramid level" representation of the images:</br>reducing their quality and initally registering the "easier" (and smaller) images before iteratively using the best found solution on the next larger image
-  Solving an initial cross correlation on the lowest Pyramid level:</br> Initially selecting the best found solution by cross-correlating the two smallest images together. This is usually done with only translation, It is possible to change the configuration to include rotation and scaling
-  Iterative solving:</br> After the initial cross correlation, we optimise based on your chosen metric by translating, rotating, scaling or shearing the images around. This is iteratively done, with the best found solution being passed down the pyramid back to the original images

## Features:

As a user you have full control over the optimisation algorithm and metric. A configuration file is added to save your preferences as well as provide standard settings, but there is a choice between different metrics and optimsation algorithms, as well as their parameters. </br>

When selecting more than two images, two options are possible:

-  With Pivot:</br> One image is selected to remain the fixed image, iterate through the list of images, solving between the fixed image and the n<sup>th</sup> image. The best found solution is passed onto the next image down the list
-  Without pivot: </br> The list of images is descended iteratively, choosing the current and next image as the fixed and moving image respectively

I.e. for a list of 5 images to be registered:
  
| With Pivot | Without Pivot |
| --- | --- |
| <table><tr><th>Fixed</th><th>Moving</th></tr><tr><td>Image 1</td><td>Image 2</td></tr><tr><td>Image 1</td><td>Image 3</td></tr><tr><td>Image 1</td><td>Image 4</td></tr><tr><td>Image 1</td><td>Image 5</td></tr></table> | <table><tr><th>Fixed</th><th>Moving</th></tr><tr><td>Image 1</td><td>Image 2</td></tr><tr><td>Image 2</td><td>Image 3</td></tr><tr><td>Image 3</td><td>Image 4</td></tr><tr><td>Image 4</td><td>Image 5</td></tr></table> |

## Outputs:

The scripts to run the optimisation are also provided as standalone scripts called by the app, and so can be used without needing to run the app. These scrips return affinetform2d objects, which describe an affine transformation. Functions to apply these transforms to your images are also provided and the app
natively does this as seen on the "Image" tab. 

Options to either save the transform object, or the transformed images are possible.

## GUI:

The scripts can either be run independently, or via the provided GUI. The GUI enables easy integration of the configuration settings, with customisability to allow for different image formats, optimisation settings and more.
An integrated image tab allows you to verify the optimisation as well as selectively refine unsatisfactory registrations.

## Coming Soon:

- Automatic integrated ROI selection to enable tracking of neuron activity between recordings.
- Automatic analysis on tracked ROIs including correlation of activity between recordings
