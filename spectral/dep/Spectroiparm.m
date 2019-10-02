function spar = Spectroiparm()

global DISPLAY
DISPLAY = 0;  %toggle roiselection display on or off

valid = false;
if exist('spar.mat', 'file')
     load spar.mat
     aflds = { 'cutOffHz', 'border', 'areasz', 'minimaldistance', 'pixelthreshold', 'significance', 'fractionmax', 'roundedness', 'voxel'};
     if sum(ismember(aflds, flds)) == 9
         valid = true; 
     end
end
if ~valid    
    spar.cutOffHz = 0.2; %upper value(Hz) of range of spectral components to use for roi selection
    spar.border = 15; %edge of image to ignore
    spar.areasz = [50 500];    %minimal and maximal size of the rois in number of pixels
    spar.minimaldistance = 30; %minimal distance between roi centers
    spar.pixelthreshold = 0.90; %max == 1.0 !! pixel value threshold on pixel value range to find contours
    spar.significance = 0.95; %maximum shoud be greater than significance value, within range of the pixels
    spar.fractionmax = 0.25; %higest fraction of maxima (pnt) to use
    spar.roundedness = 0.4;  %roundedness between 0 and 1.0;
    spar.voxel = 60;
end

prompt = {'cutoffhz', 'border', 'areasz [min max]', 'minimaldistance', 'pixelthreshold', 'significance',...
          'fractionmax', 'roundedness', 'voxel', 'DISPLAY(on or off)'};
dlgtitle = 'Enter parameters for roi segmentation';
defans = {num2str(spar.cutOffHz), num2str(spar.border), num2str(spar.areasz), num2str(spar.minimaldistance), ...
        num2str(spar.pixelthreshold), num2str(spar.significance), num2str(spar.fractionmax), ...
        num2str(spar.roundedness), num2str(spar.voxel), 'off'};

answer = inputdlg(prompt, dlgtitle, 1, defans);
if ~isempty(answer)
    
    spar.cutOffHz = str2double(answer{1}); 
    spar.border = str2double(answer{2}); 
    spar.areasz = str2num(answer{3});    %#ok<ST2NM>
    spar.minimaldistance = str2double(answer{4}); 
    spar.pixelthreshold = str2double(answer{5});
    spar.significance = str2double(answer{6});
    spar.fractionmax = str2double(answer{7});
    spar.roundedness = str2double(answer{8});  
    spar.voxel = str2double(answer{9});
    
    if strcmpi(answer(10), 'on')
        DISPLAY = 1;
    end
    
    save('spar.mat', 'spar')
end
