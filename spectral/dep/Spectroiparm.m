function [spar, succesBool]= Spectroiparm(spar)
% Set spar (Spectral PARameters). Without the need for a SPSIG file to show
% what the parameters mean (use SpectParArm for that).
% 
% Chris v.d. Togt
% Leander de Kraker
%
arguments
    spar = [];
end

succesBool = false;
    
% Try to find spar in current folder
if isempty(spar) & exist('spar.mat', 'file')
     load spar.mat
end
% Check if spar contains all the required fields
valid = false;
if ~isempty(spar)
    valid = SparCheckValid(spar);
end
% 
if ~valid
    spar(1).cutOffHzMin = 0.0;
    spar.cutOffHzMax = 0.2; %upper value(Hz) of range of spectral components to use for roi selection
    spar.border = 15;       %border width of image to ignore
    spar.areasz = [25 250]; %minimal and maximal size of the rois in number of pixels
    spar.roundedness = 0.9; %roundedness between 0 and 1.0;
    spar.voxel = 50;        %size of area to find roi contour in spectral image
    spar.cutOffCorr = 0.5;  %Minimum signal correlation within ROI:
    %ROIs drawn from a spectral image, in many cases, represent overlapping
    %cell bodies or dendrites. To separate these overlapping structures
    %pixel correlations are calculated within the ROI drawn from a spectral
    %image. cutoffcorr determines at which correlation strength a new contour will be selected. 
    spar.useFluorescenceImg = true; % Also use maximum and average fluorescence projection in automatic search
    spar.doPlot = false;
end

if spar.useFluorescenceImg
    defansFluo = 'yes';
else
    defansFluo = 'no';
end
if spar.doPlot
    defansDoPlot = 'yes';
else
    defansDoPlot = 'no';
end

prompt = {'cutOffHzMin: spectral cut-off Hz minimum', 'cutOffHzMax: spectral cut-off Hz maximum',...
          'border: width (pixels) around image that does not allow ROI centers [minimum = 5]',...
          'areasz: area size of ROIs [min max]', ...
          'ROI roundedness, between 0 (any shape) and 1 (perfectly round)',...
          'voxel: ROI search field size (px)',...
          'cutOffCorr: signal correlation cut-off',...
          'Also use fluorescence projection images (Max and average) in ROI search (yes or no)',...
          'Plot all individual ROI detection attempts (yes or no)'};
dlgtitle = 'Enter parameters for ROI segmentation';
defans = {num2str(spar.cutOffHzMin), num2str(spar.cutOffHzMax), num2str(spar.border), num2str(spar.areasz), ...
          num2str(spar.roundedness), num2str(spar.voxel), num2str(spar.cutOffCorr), defansFluo, defansDoPlot};

answer = inputdlg(prompt, dlgtitle, 1, defans);
if ~isempty(answer)
    spar.cutOffHzMin = str2double(answer{1}); 
    spar.cutOffHzMax = str2double(answer{2}); 
    spar.border = str2double(answer{3}); 
    spar.areasz = str2num(answer{4});    %#ok<ST2NM>
    spar.roundedness = str2double(answer{5});  
    spar.voxel = str2double(answer{6});
    spar.cutOffCorr = str2double(answer{7});
    
    if strcmpi(answer{8}, 'yes')
        spar.useFluorescenceImg = true;
    else
        spar.useFluorescenceImg = false;
    end
    
    if strcmpi(answer(9), 'yes')
        spar.doPlot = true;
    else
        spar.doPlot = false;
    end
    
    save('spar.mat', 'spar') % Save spar to current folder for future quick loading
    succesBool = true;
end
