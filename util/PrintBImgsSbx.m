% Very quickly save fluorescence images with nice aspect ratio using an sbx 
% file only
% 
% This file wants to know scaleUm (scale per pixel) and aspRat (aspect ratio).
% Run AddScaleInfo.m to get those fields into the info
% 
% 
% Leander de Kraker
% 2024-12-20
%

clear
clc

%% 1. Navigate to mouse folder and collect sbx files
removeDuplicates = false;
searchDepth = 1;
[filepaths, filenames] = CollectResFiles(removeDuplicates, '.sbx', searchDepth);
nfiles = length(filenames);

clear removeDuplicates searchDepth

%%  Where to save the pictures?

% % Save in designated folder
picFolders = '\\mvp2nin\Projects\PLInterneurons\17.20.12\Saturn\pictures_fluorescence\';
picFolders = repmat({picFolders}, nfiles);

% % Save in original folder
% picFolders = filepaths;

%% Other important settings

% Scalebar settings
scaleBarUm = 100; % How big you want your scalebar (0 for no scalebar)

percentageToUse = 5; % minumum 2 percent of recording used
maxChunks = 30; % maximum number of chunks to take (takes more time)
chunkSize = 40; % number of frames per chunk

limsPerc = [0.1, 99.983]; % percentile to cut brightness values of
doPlotForOurEyes = false;
doLineShift = true; % cannot hurt

%% 

strChannelNames = {'green', 'red'};

picNamesAverage = cell(nfiles, 2);
picNamesMax = cell(nfiles, 2);
dims = zeros(nfiles, 2);
nColors  = zeros(nfiles, 1);
nchunks = zeros(nfiles, 1);
nframes  = zeros(nfiles, 1);
scales = zeros(nfiles, 1);
aspRats = zeros(nfiles, 1);
perms = zeros(nfiles, 4);
nFramesRead = zeros(nfiles, 1);
percentageRead = zeros(nfiles, 1);
s = nan(nfiles, 1); % line shift applied on images
m = nan(nfiles, 1); 
limsAverage = nan(nfiles, 2);
limsMax = nan(nfiles, 2);

colors = {cmapL('greenFancy', 256);...
          cmapL([1 1 1; 1 0.2 0; 0 0 0], [50 206])};

%%
tic
for i = 1:nfiles
    %% Calculate the fluorescence images    
    clearvars info
    clearvars -global info
    
    filei = [filepaths{i}, filenames{i}(1:end-4)];
    imgTest = sbxread(filei, 2,1);
    global info
    
    nframes(i) = info.max_idx;
    aspRats(i) = info.pixelAspectRatio;
    scales(i) = info.scaleUm;
    perms(i,1:length(info.Perm)) = info.Perm;
    permBad = perms(i,1)==1;
    if permBad
        imgTest = permute(imgTest, [2 3 1]);
    end
    dims(i, :) = size(imgTest, [1 2]);
    nColors(i) = size(imgTest, 3);
    
    percentageOneChunk = chunkSize/nframes(i)*100;
    % take at least 1 chunk but not more than the max chunks to get to the
    % requested percentage of data used for the image
    nchunks(i) = max(1, min(maxChunks, ceil(percentageToUse/percentageOneChunk)));
    toRead = floor(linspace(2, nframes(i)-chunkSize-1, nchunks(i)));% 
    nFramesRead(i) = nchunks(i)  * chunkSize;
    percentageRead(i) = nFramesRead(i) / nframes(i) * 100;
    
    % Calculate average projection
    img = zeros(dims(i, 1), dims(i, 2), nColors(i));
    imgMax = zeros(dims(i, 1), dims(i, 2), nColors(i));
    for j = 1:nchunks(i)
        imgj = sbxread(filei, toRead(j), chunkSize);
        if permBad
            imgj = permute(imgj, [2 3 1 4]);
        end
        imgj = mean(double(imgj), 4);
        img = img + imgj;
        imgMax = max(imgMax, imgj);
    end
    img = img ./ nchunks(i);

    
    % Correct aspect ratioi = 
    if aspRats(i)~=1
        img = imresize(img, [dims(i,1)*aspRats(i), dims(i,2)]);
        imgMax = imresize(imgMax, [dims(i,1)*aspRats(i), dims(i,2)]);
    end
    
    scaleBarSize = scaleBarUm * scales(i);
    for c = 1:nColors(i)
        imgc = img(:,:,c);
        % the edges can have artifacts that fuck up brightness measurement
        imgBrightnessCheck = imgc(50:end-50, 150:end-50); 
        limsAverage(i,:) = prctile(imgBrightnessCheck(:), limsPerc);
        doPerc = false; % the values for lims represent percentile
        
        % picNames{i,c} = sprintf('%sFile%d_-_%s_-_%d-frames_average_-_%s_-_%d-%dfluorescence.png',...
        %                 picFolders{i}, i, filenames{i}(1:end-4), nFramesRead, strChannelNames{c}, round(lims(1)), round(lims(2)));
        % 
        picNames{i,c} = sprintf('%sFile%d_-_%s_-_%d-frames_average_-_%s.png',...
                        picFolders{i}, i, filenames{i}(1:end-4), nFramesRead(i), strChannelNames{c});
        
        [rgb, s(i,:), m(i,:)] = PrintBImg(imgc, limsAverage(i,:), colors{c},...
                        doPlotForOurEyes, picNames{i,c}, scaleBarSize,...
                        doPerc, doLineShift);
        
    end

    imgBrightnessCheck = imgMax(50:end-50, 150:end-50); 
    limsMax(i,:) = prctile(imgBrightnessCheck(:), limsPerc);
    picNamesMax{i} = strrep(picNames{i,1}, 'average', 'max');
    [rgbMax] = PrintBImg(imgMax, limsMax(i,:), colors{1},...
                    doPlotForOurEyes, picNamesMax{i}, scaleBarSize,...
                    doPerc, doLineShift);
    
    
    fprintf('Done with %d/%d. ETA %.1f minutes\n', i, nfiles, ETA(i, nfiles, toc/60))
end

%% save some variables
save([picFolders{1}, 'imageDataProperties.mat'], 'filepaths', 'filenames', 'nfiles',...
                                                 's', 'm', 'limsAverage', 'scales', 'perms', 'aspRats',...
                                                 'nFramesRead', 'percentageRead')

