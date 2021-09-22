function spar = Spectroiparm()
% Set spar (Spectral PARameters). Without the need for a SPSIG file to show
% what the parameters mean (use SpectParArm for that).
% 
% Chris v.d. Togt
% Leander de Kraker
%

global DISPLAY
DISPLAY = false;  %toggle roiselection display on or off

valid = false;
if exist('spar.mat', 'file')
     load spar.mat
     flds = fieldnames(spar);
     aflds = { 'cutOffHzMin', 'cutOffHzMax', 'border', 'areasz',...
               'roundedness', 'voxel', 'cutOffCorr', 'useFluorescenceImg'};
     if sum(ismember(aflds, flds)) == 8
         valid = true; 
     end
end
if ~valid    
    spar.cutOffHzMin = 0.0;
    spar.cutOffHzMax = 0.2; %upper value(Hz) of range of spectral components to use for roi selection
    spar.border = 15;       %border width of image to ignore
    spar.areasz = [25 250]; %minimal and maximal size of the rois in number of pixels
    spar.roundedness = 0.9; %roundedness between 0 and 1.0;
    spar.voxel = 50;        %size of area to find roi contour in spectral image
    spar.cutOffCorr = 0.5;  %Fraction of maximal correlation in ROI:
    spar.useFluorescenceImg = false; % Also use maximum and average fluorescence projection in automatic search
    %ROIs drawn from a spectral image, in many cases, represent overlapping
    %cell bodies or dendrites. To separate these overlapping structures
    %pixel correlations are calculated within the ROI drawn from a spectral
    %image. cutoffcorr determines at which correlation strength a new contour will be selected. 

end
% 
% figure()
% subplot(2,1,1)
if spar.useFluorescenceImg
    defansFluo = 'yes';
else
    defansFluo = 'no';
end

prompt = {'cutOffHzMin: spectral cut-off Hz minimum', 'cutOffHzMax: spectral cut-off Hz maximum',...
          'border: width (pixels) around image that does not allow ROI centers [minimum = 5]',...
          'areasz: area size of ROIs [min max]', ...
          'ROI roundedness, between 0 (any shape) and 1 (perfectly round)',...
          'voxel: ROI search field size (px)',...
          'cutOffCorr: signal correlation cut-off',...
          'Also use fluorescence projection images (Max and average) in ROI search (yes or no)',...
          'DISPLAY (on or off)'};
dlgtitle = 'Enter parameters for ROI segmentation';
defans = {num2str(spar.cutOffHzMin), num2str(spar.cutOffHzMax), num2str(spar.border), num2str(spar.areasz), ...
          num2str(spar.roundedness), num2str(spar.voxel), num2str(spar.cutOffCorr), defansFluo, 'off'};

answer = inputdlg(prompt, dlgtitle, 1, defans);
if ~isempty(answer)
    
    spar.cutOffHzMin = str2double(answer{1}); 
    spar.cutOffHzMax = str2double(answer{2}); 
    spar.border = str2double(answer{3}); 
    spar.areasz = str2num(answer{4});    %#ok<ST2NM>
    spar.roundedness = str2double(answer{5});  
    spar.voxel = str2double(answer{6});
    spar.cutOffCorr = str2double(answer{7});
    
    if strcmpi(answer{8}, 'no')
        spar.useFluorescenceImg = false;
    else
        spar.useFluorescenceImg = true;
    end
    
    if strcmpi(answer(9), 'on')
        DISPLAY = true;
    end
    
    save('spar.mat', 'spar')
end
