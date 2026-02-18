%
% Automatic mask placing from a recording that has a mask to other
% maskless recordings who have to get a Mask
% 
% This code also puts the variable BImg (spectral image) in the files that 
%   get a Mask when BImg is absent in that file.
% If BImg is absent in a file that has to get a mask an attempt is made to
%   create the BImg, if that attempt failes the maximum projection 
%   (BImgMax) will be loaded from that dataset, BImgMax will then be saved
%   as BImg into that file. Because some scripts require the variable BImg
% 
%
% If results are not good: use PlaceMask.m for manual placement!
% 
% Leander de Kraker
% 2020-3-6
%

%% 1. Load recording that has the mask that all other recordings should have
clear all
close all

[fnSource, pnSource] = uigetfile('*SPSIG.mat', 'get file that has the mask');
load([pnSource fnSource], 'Mask', 'BImg', 'PP', 'BImgMax', 'BImgAverage')
sourceMask = Mask;
sourcePP = PP;

useBImgMax = false; % Force the use of BImgMax
if ~useBImgMax
	sourceBImg = BImg;
    fprintf('using BImg (spectral) when possible\n')
else
    sourceBImg = BImgMax;
    fprintf('going to use BImgMax only (fluorescence max projection)\n')
end

%% 2. Load all the recordings that should get the mask

filenames = {};
filepaths = {};
selecting = true; % true as long as files are being selected
i = 0;
while selecting
    i = i + 1;
    if i == 1
        str = sprintf('load spectrally analyzed file %d', i);
    else
        str = sprintf('load spectrally analyzed file %d. previous: %s. Press cancel when done'...
            , i, filenames{i-1});
    end
    [filenames{i}, filepaths{i}] = uigetfile('*SPSIG.mat', str);
    
    if filenames{i} == 0 % Cancel is pressed probably: stop with selecting
        filenames(i) = [];
        filepaths(i) = [];
        selecting = false;
    end
end
filenames = filenames';
filepaths = filepaths';
nfiles = length(filenames); % The number of files that have been selected
fprintf('selected %d files\n\n', nfiles)
clearvars i selecting str


%% 3. Load the Masks and (create) BImgs

% Load the cutoffHz from the spar in this folder
sparfile = [pnSource, 'spar.mat'];
if exist(sparfile, 'file')~=0
    load(sparfile)
    if isfield(spar, 'cutOffHz')
        cutOffHzMax = spar.cutOffHz;
        cutOffHzMin = 0;
    elseif isfield(spar, 'cutOffHzMax')
        cutOffHzMax = spar.cutOffHzMax;
        cutOffHzMin = spar.cutOffHzMin;
    end
else
    cutOffHzMax = 0.4;
    cutOffHzMin = 0;
end

BImgs = cell(nfiles, 1);
hasMask = false(nfiles,1);
hasBImg = false(nfiles,1);
for i = 1:nfiles
    clearvars BImg BImgMax Mask SPic Sax Spect Slow % clear variables from previous iterations
    load([filepaths{i} filenames{i}], 'Mask', 'BImg', 'BImgMax', 'SPic', 'Sax')
    
    % Check if the recording has a mask already
    if exist('Mask', 'var')
        fprintf('!file %d, %s already has a mask with %d ROIs! (source mask has %d ROIs)\n',...
            i, filenames{i}(1:end-4), length(unique(Mask(:)))-1, max(sourceMask(:)))
        hasMask(i) = true;
    end
    
    % Get a BImg for this recoring
    if useBImgMax % force use of BImgMax if requested
        if exist('BImgMax', 'var')
            BImgs{i} = BImgMax;
            fprintf('..used max projection\n')
        else
            warning('useBImgMax was set to true, but there is no BImgMax to use!')
        end
        
        if exist('BImg', 'var') % check whether to save BImgMax as BImg
            hasBImg(i) = true; % if BImg is already present then no need to save it
        end
	else
        % Create BImg if not present, otherwise use BImgMax
        if ~exist('BImg', 'var') 
            fprintf('creating BImg for file %d\n', i)
            if sum(isnan(SPic(:))) < (length(SPic)/5) % If there are not too many nans
                imgStack = log(SPic(:,:,2:end));
                Sax(1) = [];
                selectedHz = (Sax <= cutOffHzMax) & (Sax >= cutOffHzMin);
                Spect = imgStack(:,:,selectedHz);
                Spect = permute(Spect, [2 1 3]); %transposes the images,
                Spect = setminlevel(Spect); %set minimum level to zero
                BImg = max(Spect, [], 3);
            else % Try to use maximum projection instead of spectral
                fprintf('file %d has %d NaNs in the %d SPic values.. ùsing BImgMax (maximum fluorescence)..\n',...
                    i, sum(isnan(SPic(:))), length(SPic(:)))
                if exist('BImgMax', 'var')
                    BImg = BImgMax;
                    fprintf('..used max projection\n')
                else
                    warning('no way to create image for file %d. quiting', i)
                    fprintf('Advice: Run BackgroundImgSbx.m to create BImgMax for this recording\n\n')
                    return
                end
            end
        else
            hasBImg(i) = true;
            fprintf('file %d already has BImg\n', i)
        end

        BImgs{i} = BImg;
    end
end
fprintf('done\n\n')


%% 4. Register the BImgs and apply register 

sourceBImgReg = cell(nfiles,1);
transforms = struct('x',[],'y',[],'rotAng',[]);
SpatialCorrs = cell(nfiles,1);
Masks = cell(nfiles,1);
missingROIs = cell(nfiles,1);
nMissing = zeros(nfiles, 1);
PPs = sourcePP;answer = 'everything fine'; % in case ROIs are missing
tic
for i = 1:nfiles
    % Register the source BImg to the recording to place the mask on
    [sourceBImgReg{i}, transforms(i)] = Register2Imgs(BImgs{i}, sourceBImg);
    transforms(i).x = -transforms(i).x;
    transforms(i).y = -transforms(i).y;
%     transforms(i).rotAng = -transforms(i).rotAng;
    fprintf('done with file %d: time elapsed %.0f seconds.\n', i, toc)
    
    % Apply registration to the Mask
    Masks{i} = RegistrationApply(sourceMask, transforms(i));
    PPs(i)   = RegistrationApplyPP(sourcePP, transforms(i), size(Masks{i}));
    
    % Check if all ROIs are still present
    presentROIs = unique(Masks{i}(:));
    presentROIs(presentROIs==0) = [];
    missingROIs{i} = 1:PPs(i).Cnt;
    missingROIs{i}(presentROIs) = [];
    nMissing(i) = length(missingROIs{i});
    if nMissing(i)>0 && strcmp(answer, 'everything fine')
        answer = questdlg({sprintf('%d ROIs are missing and this WILL cause problems. You want the same ROIs in every recording.%s', nMissing(i),...
              ' Given the original set of ROIs this is impossible according to current registration');...
              '';...
              'Delete ROIs from ORIGINAL recording?? Overwriting that file (including signals)';...
              'If you choose delete & overwrite, rerun PlaceMaskAuto (from 3rd section onward if you have the variables in workspace)'},...
              sprintf('Missing ROIs found in rec %d', i),...
              'Delete & overwrite', 'ignore problems', 'Delete & overwrite');
    end

    % Plot
    figure
    imagesc(CreateRGB({sourceBImgReg{i}, BImgs{i}, Masks{i}>0},'r g b'))
    hold on
    for roi = 1:PPs(i).Cnt
        plot(PPs(i).Con(roi).x, PPs(i).Con(roi).y, 'w')
    end
    title(sprintf('%d, correlation %.3f\nred=registered sourceBImg, green=rec %d, blue=registered Mask',...
        i, corr(sourceBImgReg{i}(:), BImgs{i}(:)),i))
    pause(0.01)

    if ~strcmp(answer, 'Delete & overwrite')
        % Calculate the SpatialCorr and PP.Rvar specifically for this recording
        [sbxt, ~, freq] = transmemap([filepaths{i} filenames{i}(1:end-9) 'Trans.dat']);
        fprintf('calculating spatialcorr and PP.Rvar...')
        SpatialCorrs{i} = zeros(size(Masks{i}));
        for roi = 1:PPs(i).Cnt
            [cor, roiidx, PPs(i).Rvar(roi)] = SpatialCorrCalcFun(sbxt, freq, Masks{i}, roi, PPs(i).P([1 2],:), true);
            SpatialCorrs{i}(roiidx) = cor(roiidx);
        end
        fprintf(' done\n')
        
    %     % Also update PPs.Rvar
    %     [PPs(i).SpecProfile, PPs(i).peakFreq] = SpecProfileCalcFun(imgStackT{i}, Masks{i}, 1:PP.Cnt, data(i).Sax);
    end
end
fprintf('done!\n')

if strcmp(answer, 'Delete & overwrite')
    % Remove signals from original template file
    missingROIsAll = unique(cell2mat(missingROIs));
    RemoveSigs([pnSource fnSource], missingROIsAll)
    % Replace Mask & PP from original template file
    [Mask, PP] = RemoveROIs(sourceMask, sourcePP, missingROIsAll);
    save(fileTemplate, 'Mask', 'PP', '-append')
    fprintf('Overwritten Mask and PP from source file %s\n', fnSource)
end

% clearvars BImgAverage BImg BImgMax i PP cor roiidx
clearvars BImgAverage BImgMax i PP cor roiidx

%% 5. Save to the Masks and other info to the files

for i = 1:nfiles
    clearvars Mask PP SpatialCorr
    if hasMask(i)
        fprintf('This recording already has a mask that will now be overwritten.\n')
    end
    
    Mask = Masks{i};
    PP = PPs(i);
    SpatialCorr = SpatialCorrs{i};
    fprintf('Saving file %d...', i)
    save([filepaths{i} filenames{i}], 'Mask', 'PP', 'SpatialCorr','-append')

    if ~hasBImg(i) % If BImg was not present in file, save it too
        fprintf('. including new BImg...')
        BImg = BImgs{i};
        save([filepaths{i} filenames{i}], 'BImg','-append')
    end
    fprintf(' saved\n')
end
