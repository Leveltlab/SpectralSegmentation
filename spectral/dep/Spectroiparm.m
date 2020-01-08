function spar = Spectroiparm()

global DISPLAY
DISPLAY = 0;  %toggle roiselection display on or off

valid = false;
if exist('spar.mat', 'file')
     load spar.mat
     flds = fieldnames(spar);
     aflds = { 'cutOffHz', 'border', 'areasz', 'roundedness', 'voxel', 'cutoffcorr'};
     if sum(ismember(aflds, flds)) == 6
         valid = true; 
     end
end
if ~valid    
    spar.cutOffHz = 0.2;    %upper value(Hz) of range of spectral components to use for roi selection
    spar.border = 15;       %border width of image to ignore
    spar.areasz = [25 250]; %minimal and maximal size of the rois in number of pixels
    spar.roundedness = 0.4; %roundedness between 0 and 1.0;
    spar.voxel = 50;        %size of area to find roi contour in spectral image
    spar.cutoffcorr = 0.5;  %Fraction of maximal correlation in ROI:
    %ROIs drawn from a spectral image, in many cases, represent overlapping
    %cell bodies or dendrites. To separate these overlapping structures
    %pixel correlations are calculated within the ROI drawn from a spectral
    %image. cutoffcorr determines at which correlation strength a new contour will be selected. 

end

prompt = {'cutoffhz', 'border', 'areasz [min max]', ...
          'roundedness', 'voxel', 'cutoffcorr', 'DISPLAY(on or off)'};
dlgtitle = 'Enter parameters for roi segmentation';
defans = {num2str(spar.cutOffHz), num2str(spar.border), num2str(spar.areasz), ...
          num2str(spar.roundedness), num2str(spar.voxel), num2str(spar.cutoffcorr), 'off'};

answer = inputdlg(prompt, dlgtitle, 1, defans);
if ~isempty(answer)
    
    spar.cutOffHz = str2double(answer{1}); 
    spar.border = str2double(answer{2}); 
    spar.areasz = str2num(answer{3});    %#ok<ST2NM>
    spar.roundedness = str2double(answer{4});  
    spar.voxel = str2double(answer{5});
    spar.cutoffcorr = str2double(answer{6});
    
    if strcmpi(answer(6), 'on')
        DISPLAY = 1;
    end
    
    save('spar.mat', 'spar')
end
