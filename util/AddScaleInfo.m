% Add scaleUm and pixelAspectRatio to SPSIG and/or to info of sbx files
% This script is made for GAIA scanner which images 512 x 796 pixels. and
% has a field of view (FOV) of 1000um when zoomlevel is 1.
% 
% It also recalculated roundedness of ROIs with new pixel aspect ratio.
% 
% - scaleUm: physical width of one pixel in micrometer (um)
% - pixelAspectRatio: How much higher a pixel is than it is wide.
%     example: a physically square image recorded with dimensions 
%        16x9 pixels (width x height) results in pixel aspect ratio 
%        16/9 = 1.7778. The height of the pixel covers more physical space.
% 
% Leander de Kraker
% 2024-11-4
% 

% Select sbx files, or select spsig files. Later script will try to find
% accompanying files
%% 1. Navigate to mouse folder and collect sbx files
removeDuplicates = false;
searchDepth = 1;
[filepaths, filenames] = CollectResFiles(removeDuplicates, '\*.sbx', searchDepth);
nfiles = length(filenames);

clear removeDuplicates searchDepth
%% 1. Or navigate to mouse folder and collect SPSIG files
[filepaths, filenames] = CollectResFiles(false, '\*SPSIG.mat', 1);
nfiles = length(filenames);


%% 1. Or find SPSIG files via chronic file
[chronicName, chronicPath] = uigetfile('*chronic.mat');
load([chronicPath, chronicName], 'filepaths', 'filenames', 'filedates', 'nfiles')
filenames = strrep(filenames, '.mat', '_Res.mat');


%% 2. Collect accompanying files

clearvars global
clearvars -global info
clearvars info

if strcmp(filenames{1}(end-3:end), '.sbx') % Add SPSIG files and check their existence
    filenamesSbx = strrep(filenames, '.sbx', '.mat');
    filenamesSPSIG = strrep(filenames, '.sbx', '_SPSIG.mat');

    presentSbx = true(nfiles, 1);
    presentSPSIG = false(nfiles, 1);
    for i = 1:nfiles
        presentSPSIG(i) = exist([filepaths{i}, filenamesSPSIG{i}], 'file')==2;
    end
else
    filenamesSPSIG = filenames;
    filenamesSbx = strrep(filenamesSPSIG, '_SPSIG.mat', '.mat');
    
    presentSPSIG = false(nfiles, 1);
    presentSbx = false(nfiles, 1);
    for i = 1:nfiles
        presentSbx(i) = exist([filepaths{i}, filenamesSbx{i}], 'file')==2;
    end
end
fprintf('(%d/%d) files had accompanying file\n', sum(presentSbx|presentSPSIG), nfiles)

%% 3. Add scale info to the files


doSaving = true; % Actually update the files? 
doPlot = false; % plot and print stats on recalculation of roundedness

% The theoretical pixel aspect ratio seemed nice as checked by
% BeadsCalibration.m and RoundednessCheck.m
pixelAspectRatioExpected = 796/512;
umAt1xZoom = 1000; % ! Assuming this number of micrometers when zoom=1  

changed  = true(nfiles, 1); % Has new info been found?
scaleUms = nan(nfiles, 1); % pixel sizes filled into each file
asRats   = nan(nfiles, 1); % aspect ratios filled into each file
freqs    = nan(nfiles, 1); % sampling rate Hz
dims      = nan(nfiles, 2); % size of image according to sbx info
zoomLevels = nan(nfiles, 1); % zoomlevels according to sbx info
scaleUmOldSbx = nan(nfiles, 1);
asRatsOldSbx  = nan(nfiles, 1);
scaleUmOldSpsig = nan(nfiles, 1);
asRatsOldSpsig  = nan(nfiles, 1);
% Preallocate question dialog variables
defInput = {};
defInput{4} = pixelAspectRatioExpected;
FOVum = cell(nfiles, 1);
squareFOV = cell(nfiles, 1);

for i = 1:nfiles
    fprintf('\nFile %d: %s\n', i, filepaths{i})
    % Forget previous iteration a bit.
    clear info
    
    % Get the zoom
    zoomFound = false;
    if presentSbx(i) 
        load([filepaths{i}, filenamesSbx{i}], 'info')
        dims(i,:) = [info.d1, info.d2];
        
        if isfield(info, 'scaleUm')
            scaleUmOldSbx(i) = info.scaleUm;
        end
        if isfield(info, 'pixelAspectRatio')
            asRatsOldSbx(i) = info.pixelAspectRatio;
        end
        
        if isfield(info, 'config') && isfield(info.config, 'magnification')
            zoomLevels(i) = str2double(info.config.magnification_list(info.config.magnification, :));
            scaleUms(i) = info.d1 / (umAt1xZoom / zoomLevels(i)); % original number of pixels / µm full width of FOV
            zoomFound = true;
        end
    end
    
    if presentSPSIG(i)
        spsig = load([filepaths{i}, filenamesSPSIG{i}], 'scaleUm', 'freq', 'pixelAspectRatio', 'PP');
        freqs(i) = spsig.freq;
        defInput{1} = spsig.freq;
        if isfield(spsig, 'pixelAspectRatio')
            asRatsOldSpsig(i) = spsig.pixelAspectRatio;
        end
        
        if isfield(spsig, 'scaleUm')  % 
            scaleUmOldSpsig(i) = spsig.scaleUm;
            fprintf('scaleUm variable already found in SPSIG file\n')
            if isempty(defInput{2})
                defInput{2} = spsig.scaleUm;
            end
        end
    end
    
    asRats(i) = pixelAspectRatioExpected;
    
    % If not sure about zoom ask. (can overwrite pixel aspect ratio)
    if ~zoomFound
        sAnswer = struct();
        [hzi, scaleUmi, FOVum{i}, pixelAspectRatioi, squareFOV{i}, defInput] = RequestRecInfo(defInput);
        sAnswer = RequestRecInfoProcess(sAnswer, hzi, scaleUmi, FOVum{i}, pixelAspectRatioi, squareFOV{i});
        asRats(i) = sAnswer.pixelAspectRatio;
        scaleUms(i) = sAnswer.scaleUm;
    end
    
    % Recalculate ROI roundedness based on pixelaspect ratio
    sameRoundedness = false;
    if presentSPSIG(i) && isfield(spsig, 'PP')
        oldRoundedness = spsig.PP.Roundedness;
        for r = 1:spsig.PP.Cnt
            xr = spsig.PP.Con(r).x;
            yr = spsig.PP.Con(r).y;
            spsig.PP.Roundedness(r) = perimarea(xr, yr*asRats(i));
        end
        % if doPlot
        %     figure
        %     values = [oldRoundedness', spsig.PP.Roundedness'];
        %     groups = {'old'; 'new'};
        %     dotSize = NaN;
        %     toDraw = 1:spsig.PP.Cnt;
        %     colors = cmapL('viridis', spsig.PP.Cnt);
        %     boxplot_w_paireddata(values, groups, dotSize, toDraw, colors);
        %     ylabel('ROI roundedness')
        %     title({sprintf('%s', strrep(filenamesSPSIG{i}, '_', '\_'));...
        %            sprintf('pixel aspect ratio from %.2f to %.2f', asRatsOldSpsig(i), asRats(i));...
        %            sprintf('roundedness: %.2f to %.2f', mean(oldRoundedness), mean(spsig.PP.Roundedness))})
        % end
        fprintf('Old roundedness mean %.2f (+-%.2fstd). pixel aspect %.3f\n',...
                mean(oldRoundedness), std(oldRoundedness), asRatsOldSpsig(i))
        fprintf('New roundedness mean %.2f (+-%.2fstd). pixel aspect %.3f\n',...
                mean(spsig.PP.Roundedness), std(spsig.PP.Roundedness), asRats(i))
        [~, p, ~, stats] = ttest2(oldRoundedness, spsig.PP.Roundedness);
        fprintf('t(%d)=%.2f, p=%.4f\n', stats.df, stats.tstat, p)
        sameRoundedness = all(oldRoundedness == spsig.PP.Roundedness);
    end
    % If not info has been changed don't save
    changed(i) = (asRats(i) ~= asRatsOldSpsig(i) & asRats(i) ~= asRatsOldSbx(i))...
                | (scaleUms(i) ~= scaleUmOldSpsig(i) & scaleUms(i) ~= scaleUmOldSbx(i))...
                | ~sameRoundedness;
    if ~changed(i); fprintf('scale and pixel aspect ratio unchanged, not saving\n'); end
    
    
    % Save to the sbx and SPSIG files
    if doSaving && changed(i)
        if presentSbx(i)
            load([filepaths{i}, filenamesSbx{i}], 'info')
            info.scaleUm = scaleUms(i);
            info.pixelAspectRatio = asRats(i);
            save([filepaths{i}, filenamesSbx{i}], 'info', '-append')
            fprintf('Saved updated info to Sbx\n')
        end

        if presentSPSIG(i)
            spsig.scaleUm = scaleUms(i);
            spsig.pixelAspectRatio = asRats(i);
            save([filepaths{i}, filenamesSPSIG{i}], '-struct', 'spsig', '-append')
            fprintf('Saved variables to SPSIG\n')
        end
    end
end
