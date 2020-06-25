function spar = Spectroiparm()

global DISPLAY
DISPLAY = 0;  %toggle roiselection display on or off

valid = false;
if exist('spar.mat', 'file')
     load spar.mat
     flds = fieldnames(spar);
     aflds = { 'cutOffHzMin', 'cutOffHzMax', 'border', 'areasz', 'roundedness', 'voxel', 'cutOffCorr'};
     if sum(ismember(aflds, flds)) == 6
         valid = true; 
     end
end
if ~valid    
    spar.cutOffHzMin = 0.0;
    spar.cutOffHzMax = 0.2;    %upper value(Hz) of range of spectral components to use for roi selection
    spar.border = 15;       %border width of image to ignore
    spar.areasz = [25 250]; %minimal and maximal size of the rois in number of pixels
    spar.roundedness = 0.9; %roundedness between 0 and 1.0;
    spar.voxel = 50;        %size of area to find roi contour in spectral image
    spar.cutOffCorr = 0.5;  %Fraction of maximal correlation in ROI:
    %ROIs drawn from a spectral image, in many cases, represent overlapping
    %cell bodies or dendrites. To separate these overlapping structures
    %pixel correlations are calculated within the ROI drawn from a spectral
    %image. cutoffcorr determines at which correlation strength a new contour will be selected. 

end
% 
% figure()
% subplot(2,1,1)

prompt = {'cutoffhzMin', 'cutoffhzMax', 'border', 'areasz [min max]', ...
          'roundedness', 'voxel', 'cutoffcorr', 'DISPLAY(on or off)'};
dlgtitle = 'Enter parameters for roi segmentation';
defans = {num2str(spar.cutOffHzMin), num2str(spar.cutOffHzMax), num2str(spar.border), num2str(spar.areasz), ...
          num2str(spar.roundedness), num2str(spar.voxel), num2str(spar.cutOffCorr), 'off'};

answer = inputdlg(prompt, dlgtitle, 1, defans);
if ~isempty(answer)
    
    spar.cutOffHzMin = str2double(answer{1}); 
    spar.cutOffHzMax = str2double(answer{2}); 
    spar.border = str2double(answer{3}); 
    spar.areasz = str2num(answer{4});    %#ok<ST2NM>
    spar.roundedness = str2double(answer{5});  
    spar.voxel = str2double(answer{6});
    spar.cutOffCorr = str2double(answer{7});
    
    if strcmpi(answer(7), 'on')
        DISPLAY = 1;
    end
    
    save('spar.mat', 'spar')
end
