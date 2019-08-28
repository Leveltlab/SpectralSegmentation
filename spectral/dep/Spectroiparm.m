function spar = Spectroiparm()


if exist('spar.mat', 'file')
    load spar.mat
else
    
    spar.border = 15; %edge of image to ignore
    spar.areasz = [60 1000];    %minimal and maximal size of the rois in number of pixels
    spar.minimaldistance = 30; %minimal distance between roi centers
    spar.pixelthreshold = 0.90; %max == 1.0 !! pixel value threshold on pixel value range to find contours
    spar.significance = 0.95; %maximum shoud be greater than significance value, within range of the pixels
    spar.fractionmax = 0.25; %higest fraction of maxima (pnt) to use
    spar.roundedness = 0.4;  %roundedness between 0 and 1.0;
    spar.voxel = 250;
end

prompt = {'border', 'areasz [min max]', 'minimaldistance', 'pixelthreshold', 'significance',...
          'fractionmax', 'roundedness', 'voxel'};
dlgtitle = 'Enter parameters for roi segmentation';
defans = {num2str(spar.border), num2str(spar.areasz), num2str(spar.minimaldistance), ...
        num2str(spar.pixelthreshold), num2str(spar.significance), num2str(spar.fractionmax), ...
        num2str(spar.roundedness), num2str(spar.voxel)};

answer = inputdlg(prompt, dlgtitle, 1, defans);
if ~isempty(answer)
    
    spar.border = str2double(answer{1}); 
    spar.areasz = str2num(answer{2});    %#ok<ST2NM>
    spar.minimaldistance = str2double(answer{3}); 
    spar.pixelthreshold = str2double(answer{4});
    spar.significance = str2double(answer{5});
    spar.fractionmax = str2double(answer{6});
    spar.roundedness = str2double(answer{7});  
    spar.voxel = str2double(answer{8});
    
    save('spar.mat', 'spar')
end
