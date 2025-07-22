% Remove a recording from a chronic file
% because you might want to do that for some reason
% 
% Leander de Kraker
% 2025-4-2
% 

clear
clc

%%
[filenameChronic, pathnameChronic] = uigetfile('*chronic.mat');
load([pathnameChronic, filenameChronic]);

%% State which recording to delete

recToDel = 3;

%% Delete all the things

BImgs(recToDel) = [];
Masks(recToDel) = [];
PPs(recToDel) = [];
filedates(recToDel) = [];
filenames(recToDel) = [];
filenamesShort(recToDel) = [];
if exist('recLabels', 'var')
    recLabels(recToDel) = [];
end

transformed.vert(recToDel, :) = [];
transformed.hori(recToDel, :) = [];
transformed.dims(recToDel, :) = [];
transformed.rotation(recToDel, :) = [];

% Simply recalculate the overlaps because that's the easiest
inRoi = CalcRoiOverlap(Masks);
linked = LinkOverlappingRois(inRoi, thres);

% Delete from the linkMat
linkMat(:,recToDel) = [];
nLinks = sum(linkMat~=0,2);
matchToDel = nLinks==1;
linkMat(matchToDel,:) = [];
nLinks(matchToDel) = [];

% Delete from the optional variables
% imgs (MoveLikeChronic)
if exist('imgs', 'var')
    fields = fieldnames(imgs);
    for i = 1:length(fields)
        imgs.(fields{i})(recToDel) = [];
    end
end
% red related variables (ChronicRedRois)
if exist('redConfirmed', 'var')
    redLinkMat(:,recToDel) = [];
    redLinkMat(matchToDel,:) = [];
    redROIsChronic(recToDel) = [];
    redConfirmed(matchToDel) = [];
end

nfiles = length(filenames);

%% Save new chronic file