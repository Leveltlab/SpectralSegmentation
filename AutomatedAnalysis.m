function timedData = AutomatedAnalysis(do, files, spar)
% Execute all preprocessing steps from raw sbx file to RoiManagerGUI step.
% 
% Doing the normcorr motion correction, transposing, spectral analysis and
% automated ROI creation for all splits, for multiple files in one script
% Optional steps are motion correction, background subtraction, ROI getting
% 
% 
% Usage:
% 1. Load recordings by running this first section.
% 2(option 1). Set analysis settings either via the script by changing 
%              the code via "Settings via script" (section 2)
%  (option 2). Set all relevant analysis settings via walkthrough pop-up
%              prompts via "Get settings via input dialogs" (section 3)
% 3. Run analysis by running section "Process all the files" (section 5)
% 
% 
% Leander de Kraker
arguments
    do = [];
    files = [];
end

% Set the files to analyse
if isempty(files) % Select files via pop up if nothing is given
    [files, nfiles] = SelectFiles();

elseif ischar(files) || isstring(files) % if text is given
    [filesTxtFolder, filesTxtName, filesTxtExt] = fileparts(files);
    
    % If the input was the location to a text file, read it to extract the foldernames filenames
    if strcmpi(filesTxtExt, 'txt')
        filesTxt = files;
        files = readtable(filesTxt, 'ReadVariableNames', false, 'TextType', 'string', 'Delimiter', '\n');
        [files.path, files.name] = fileparts(files.path);
    elseif ismember(filesTxtExt, '.sbx') % If it was a mat file, load the mat file
        analyseThisFile = true;
    elseif ismember()
    end
end
if isstruct(files)
    
end

if isempty(do)
    SetAnalysisSettings
end

global info

%% Settings via script

% CROPPING SETTINGS % % % % % (for motion correction!)
% Most important part of cropping is removing the lines (on the left
% and right side) which have value 653357 and will thus overpower all 
% relevant data during motion correction.
do.motion.crop.x = 3:793; % HORIZONTAL CROP
do.motion.crop.y = 5:510;  % VERTICAL CROP

% Do line shifting correction? % % % % % % % % %
do.LineShift = false;
do.LineShiftTrans = false; % false=correct for horizontal lines, true=correct for vertical lines
do.lineShiftValue = 0;

% Run motion normcorre correction? % % % % % % % % % % % % %
do.NoRMCorre = false; % (if false, cropping settings aren't used)
do.motion.alignMethod = 'Rigid'; % options: 'Nonrigid' | 'Rigid'

% Convert *_eye.mat file into .mp4 file?
do.EyeConvert = true;
do.EyeInvert = true; % Invert black-white?
do.EyeDeFlicker = false; % Diminish flickering brightness between subsequent frames?

files.pathOut = files.path;

% BACKGROUND SUBTRACTION parameters % % % % % Necessary for 1 photon data
do.BackgroundSubtract = false;
do.BackgroundSub.FilterRadius = 55; % The radius of the filter for background subtraction
do.BackgroundSub.filterMethod = 'Circular Average'; % Options Gaussian Average | Circular Average
do.BackgroundSub.Shifter = 9000; % The data values has to be increased to prevent underexposure/ over correction
% 1D gaussian smoothing to slightly decrease vertical or horizontal banding noise. 
do.BackgroundSub.SmoothDim = 1; % options: 0 (no smoothing), 1 (horizontal smoothing), 2 (vertical smoothing)
do.BackgroundSub.SmoothSe  = 1.75; % size for the smoothing
do.BackgroundSub.plotter = true; % plot last frame

% Also run getSpectrois to detect the ROIs? % % % % % % % %
do.getROIs = false;
if do.getROIs
    % Arm the spar: Spectral roi finder PARameters
    global spar
    spar = Spectroiparm();
end

% Save timing info? % % % % % % % % % % % % % % % % % % % % 
do.timed = true;
if do.timed
    % File in which timing info is stored
    do.timedFile = 'E:\data\spectralTimed4.mat';
    if ~isfile(do.timedFile)
        warning('%s does not exist!! please say which file to add the timing data to!', do.timedFile)
        return
% %         or run following code to create a new timedFile there: 
        timedData = struct('computerName',{},'fileSize',[],'filePath', [], 'fileName',{},...
                         'nSplits', [], 'error', {}, 'memory', {}, 'alignMethod', {},...
                         'lineShiftT', [], 'normcorrT',[],'transposeT',[],'decimatT',[],...
                         'spectralT',[],'fluorescenceT',[], 'backSubT', [], 'getRoisT', [],...
                         'nRois', [], 'date', {});
        save(do.timedFile, 'timedData')
    end
end


%% Get settings via input dialogs
dlgdims = [1 85];
clearvars info; clearvars -global info

% Do line shift? %
do.LineShift = questdlg('Do line shifting? to correct bi-directional scanning misalignment', 'line shift?',...
                       'yes', 'no', 'no');
if strcmp(do.LineShift, 'yes')
    do.LineShift = true;
    do.LineShiftTrans = questdlg('correct for horizontal or vertical lines?', 'orientation of lines?',...
                                'horizontal', 'vertical', 'horizontal');
    if strcmp(do.LineShiftTrans, 'horizontal')
        do.LineShiftTrans = false;
    else
        do.LineShiftTrans = true;
    end
elseif strcmp(do.LineShift, 'no')
    do.LineShift = false;
else
    fprintf('Canceled selecting settings\n')
    return
end
if do.LineShift
    im = sbxread([files.path{1}, files.name{1}], 5, 100);
    do.lineShiftValue = ShiftLinesCheck(im, do.LineShiftTrans, true);
end


% Do NoRMCorre? %
%
do.NoRMCorre = questdlg('do NoRMCorre motion correction?', 'run NoRMCorre?', 'yes', 'no', 'yes');
if strcmp(do.NoRMCorre, 'yes')
    do.NoRMCorre = true;
elseif strcmp(do.NoRMCorre, 'no')
    do.NoRMCorre = false;
else
    fprintf('canceled selecting settings\n')
    return
end
do.motion.alignMethod = '';
if do.NoRMCorre
    % Get crop for NoRMCorre
    if exist('filepaths', 'var')
        load([files.path{1} files.name{1} '.mat'], 'info');
        definput = {sprintf('1:%d', info.sz(2)), sprintf('1:%d', info.sz(1))};
        clearvars info
    else
        definput = {'1:796', '1:512'};
    end
    
    prompt = {'x, horizontal crop', 'y, vertical crop'};
    dlgtitle = 'cropping settings for motion correction';
    
    answer = inputdlg(prompt, dlgtitle, dlgdims, definput);
    if isempty(answer)
        fprintf('canceling selecting settings\n')
        return
    end
    do.motion.crop.x = str2num(answer{1});
    do.motion.crop.y = str2num(answer{2});
    
    % Set rigid or non rigid
    do.motion.alignMethod = questdlg('Which registration method to use? (rigid = default)',...
                           'NoRMCorre registration method',...
                           'Rigid','Nonrigid','Rigid');
end


% Save in original file location ? % 
do.saveLocation = questdlg('Save new files in original folder? (yes = default)',...
                          'Output location',...
                          'yes', 'no (select different folders)', 'yes');
files.pathOut = files.path;
if strcmp(do.saveLocation, 'yes')
    do.saveLocation = false;
elseif strcmp(do.saveLocation, 'no (select different folders)') % select folder for every file
    do.saveLocation = true;
    outputpathi = 'lalala';
    i = 1;
    while i <= nfiles && ischar(outputpathi)
        outputpathi = uigetdir(cd, sprintf('output path for file %d: %s',i, files.name{i}));
        if ischar(outputpathi)
            files.pathOut{i} = [outputpathi '\'];
        end
        i = i + 1;
    end
else
    fprintf('Canceling selecting settings\n')
    return
end


% Do eye file format conversion? (keeps the original)
do.EyeConvert = questdlg('Convert "*_eye.mat" file into "*_eye.mp4" (keeps original)',...
                        'Create videofile for eye tracking?',...
                        'yes','no','yes');
if strcmp(do.EyeConvert, 'yes')
    do.EyeConvert = true;
    do.EyeInvert = questdlg('Invert black-white in eye video file?',...
                               'Eye video file settings',...
                               'yes', 'no', 'yes');
    do.EyeDeFlicker = questdlg('Deflicker the eye video?',...
                              'Eye video file settings',...
                               'yes', 'no', 'yes');
    if strcmp(do.EyeInvert, 'yes')
        do.EyeInvert = true;
    else
        do.EyeInvert = false;
    end
    if strcmp(do.EyeDeFlicker, 'yes')
        do.EyeDeFlicker = true;
    else
        do.EyeDeFlicker = false;
    end
elseif strcmp(do.EyeConvert, 'no')
    do.EyeConvert = false;
else
    fprintf('Canceling selecting settings\n')
    return
end
% Do background subtraction? %
do.BackgroundSubtract = questdlg('Do background subtraction (for 1-photon data)',...
                                'background subtraction settings',...
                                'yes','no','no');
if strcmp(do.BackgroundSubtract, 'yes');    do.BackgroundSubtract = true;
elseif strcmp(do.BackgroundSubtract, 'no'); do.BackgroundSubtract = false;
else; fprintf('Cancelling selecting settings\n'); return
end
% Get the specific background subtraction settings
if do.BackgroundSubtract
    prompt = {'filterRadius (actual size of filter = [radius*2, radius*2])', ...
        ['shift data to prevent underexposure/ over correction\n'...
        '(only change if BackgroundSubtract gives exposure warnings)'],...
        'Plot results and background image for final frame? (no or yes)'};
    dlgtitle = 'background subtract settings';
    definput = {'55', '9000', 'no'};
    answer = inputdlg(prompt, dlgtitle, dlgdims, definput);
    if isempty(answer)
        fprintf('stopping settings selection\n')
        return
    end
    do.BackgroundSub.FilterRadius = str2double(answer{1});
    do.BackgroundSub.Shifter = str2double(answer{2});
    if strcmpi(answer{3}, 'yes')
        do.BackgroundSub.plotter = true;
    else
        do.BackgroundSub.plotter = false;
    end
    do.BackgroundSub.filterMethod = questdlg('filtering method for background image estimation', 'background subtract settings',...
                            'Circular Average (advised)',...
                            'Gaussian Average',...
                            'Circular Average (advised)');
    if strcmp(do.BackgroundSub.filterMethod, 'Circular Average (advised)')
        do.BackgroundSub.filterMethod = 'Circular Average';
    elseif strcmp(do.BackgroundSub.filterMethod, 'Gaussian Average')
        % That's the correct string already: do nothing
    else
        fprintf('Stopping setting selection, no valid method for background subtraction\n')
        return
    end
    do.BackgroundSub.SmoothDim = questdlg('1D gaussian smoothing to diminish (vertical or horizontal banding) noise?',...
                         'background subtract setings',...
                         'no', 'horizontal smoothing', 'vertical smoothing', 'vertical smoothing');
    if strcmp(do.BackgroundSub.SmoothDim, 'no')
        do.BackgroundSub.SmoothDim = 0;
        do.BackgroundSub.SmoothSe = [];
    elseif strcmp(do.BackgroundSub.SmoothDim, 'horizontal smoothing')
        do.BackgroundSub.SmoothDim = 1;
    elseif strcmp(do.BackgroundSub.SmoothDim, 'vertical smoothing')
        do.BackgroundSub.SmoothDim = 2;
    else % canceled
        fprintf('Canceling selecting settings\n')
        return
    end
    if do.BackgroundSub.SmoothDim > 0
        answer = inputdlg('size of smoothing (standard deviation of gaussian)',...
                          'background subtract settings',...
                          dlgdims, {'1.7'});
        if isempty(answer)
            fprintf('Canceled selecting settings at final dlg\n')
            return
        end
        do.BackgroundSub.SmoothSe = str2num(answer{1});
    end
end

% get ROIs? %
do.getROIs = questdlg('automatic ROI detection?', 'get ROIs?', 'yes', 'no', 'yes');
if strcmp(do.getROIs, 'yes')
    do.getROIs = true;
    global spar
    spar = Spectroiparm();
elseif strcmp(do.getROIs, 'no')
    do.getROIs = false;
else
    fprintf('Cancelling selecting settings\n')
    return
end

% Timing info % 
do.timed = questdlg('Save processing time/log of the analysis?', 'time/log analysis?',...
    'yes', 'no', 'no');
if strcmp(do.timed, 'yes')
    do.timed = true;
elseif strcmp(do.timed, 'no')
    do.timed = false;
else
    do.timed = false;
end
if do.timed
    % File in which timing info is stored
    do.timedFile = 'D:\2Pdata\spectralTimed.mat'; % Default timedFile name! % % % 
    timedFileLinePos = dbstack;
    if ~isfile(do.timedFile) % Creating new file for timing data
        warning('%s does not exist! please say which file to add the timing data to!', do.timedFile)
        [timedFileName, timedFilePath] = uiputfile('*.mat', 'Where to save timing file', 'spectralTimed.mat');
        do.timedFile = [timedFilePath timedFileName];
        if timedFileName==0 % User didn't select a file so forget about saving timed data
            do.timed = false;
        elseif ~isfile(do.timedFile) % Save new requested file
            if isscalar(timedFileLinePos) && isscalar(timedFileLinePos.line)
                fprintf('\nPlease change default timedFile name (line %d) to the new file for future use:\n%s\n',...
                    timedFileLinePos.line-1, do.timedFile)
            else
                fprintf('\nPlease change default timedFile name to the new file for future use:\n%s\n', do.timedFile)
            end
            timedData = struct('computerName',{},'fileSize',[],'filePath', [], 'fileName',{},...
                         'nSplits', [], 'error', {}, 'memory', {}, 'alignMethod', {},...
                         'lineShiftT', [], 'normcorrT',[],'transposeT',[],'decimatT',[],...
                         'spectralT',[],'fluorescenceT',[], 'backSubT', [], 'getRoisT', [],...
                         'nRois', [], 'date', {});
            save(do.timedFile, 'timedData')
        end
    end
end

clearvars answer prompt dlgtitle dlgdims definput

%% 
% Set the i if you want to run a specific section of the analysis for a
% specific file (after things crash or are unsatisfactory for example)
i = 1;

%% Process all the files

for i = 1:nfiles
    try
        fn = files.name{i};
        pn = files.path{i};
        pno = files.pathOut{i};
        % prepare to load a different sbx file: thoroughly delete info variable
        clearvars info ME
        clearvars -global info
        timedDatai = struct('computerName',{},'fileSize',[],'filePath', [], 'fileName',{},...
                     'nSplits', [], 'error', {}, 'memory', {}, 'alignMethod', {},...
                     'lineShiftT', [], 'normcorrT',[],'transposeT',[],'decimatT',[],...
                     'spectralT',[],'fluorescenceT',[], 'backSubT', [], 'getRoisT', [],...
                     'nRois', [], 'date', {});
        nRois = 0; % number of detected ROIs for this recording
        
        
        %% Shift lines, and replace erroneously high values (>65530)
        if do.LineShift & do.lineShiftValue~=0
            tics.lineShift = tic;
            fprintf('Shifting lines & replacing high value of file %s...\n', fn)
            ShiftLinesSbx({pn fn}, do.LineShiftTrans, do.lineShiftValue, pno)
            timedDatai(1).lineShiftT = toc(tics.lineShift);
            pnUpdate = pno;
            fnUpdate = [fn '_shiftedLines'];
            clearvars info 
            clearvars -global info
        else
            pnUpdate = pn;
            fnUpdate = fn;
        end
    
        %% Load example frame and do some sbx info settings
        % An example frame of the recording (creates new global info variable)
        Img = sbxread([pnUpdate fnUpdate], 0,1);
        
        global info
        
        % set the name of the recording
        if ~isfield(info, 'strfp')
            info.strfp = [pnUpdate fnUpdate];
        elseif ~strcmp(info.strfp, [pnUpdate fnUpdate]) % inconsitent names, maybe it'll all go wrong
            info.strfp
            [pnUpdate fnUpdate]
            fprintf('expected name and name in info are inconsistent\nPROBLEM\n\n')
        else
            fprintf('everything is going well\n')
        end
        
        %number of splits
        if ~isfield(info, 'Slices')
            if isfield(info, 'otparam') && length(info.otparam)>=3
                 info.Slices = info.otparam(3);
            else
    %             % Maybe the info file is super bad, ask how many slices there are
    %             if ~isfield(info, 'Slices')
    %                 answ = inputdlg('How many Slices?', 'Slice info!!...', [1 40],  {'1'});
    %                  info.Slices = str2double(answ{1});
    %             end
    
                % The info file could have no slices because only 1 depth is imaged
                 info.Slices = 1;
            end
        end
        if size(Img,1) < 3
            info.bVers = true; % bvers means data is permuted because of GPU recording..
        else
            info.bVers = false;
        end
        
        nSlices = info.Slices;
        if ~isfield(info, 'bsplit')
            if nSlices > 1 % bsplit true means there are more than 1 split
                info.bsplit = true;
            else
                info.bsplit = false;
            end
        end
        
        
        %% Normcorr    
        % Download normcorre code from flatironinstitute from github.
        % Place normcorre folder in Matlab path
        
        % Show the data before committing the analysis
        Img = sbxread([pnUpdate fnUpdate], 0,100*nSlices);
        %determine pixel value range
        cLimits = prctile(Img(:), [0.5 97]);
        % Activate much tighter subplots
        % [subplot margin top&side],[figure bottomspace,topspace],[leftspace,rightspace]
        subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.01], [0.01 0.04], [0.1 0.01]);
        figure('units','normalized','position',[0.1 0.15 0.25 0.7]);
        handles = gobjects(1,nSlices);
        for j = 1:nSlices
            handles(j) = subplot(nSlices, 1, j);
            if info.bVers
                imagesc(mean(squeeze(Img(1,:,:,j:nSlices:end)),3))
            else
                imagesc(mean(squeeze(Img(:,:,1,j:nSlices:end)),3))
            end
            caxis(cLimits)
            if do.NoRMCorre
                hold on
                rectangle('position',[x(1), y(1), x(end)-x(1), y(end)-y(1)],...
                          'edgecolor',[1 1 0.75],'linewidth',2.5)
            end
            title(sprintf('slice %d',j))
        end
        linkaxes(handles, 'xy')
        colormap(cmapL('greenFancy', 256))
        drawnow
        
        if do.NoRMCorre
            % critical info for NoRMCorre
            info.crop.x = do.motion.crop.x; 
            info.crop.y = do.motion.crop.y;
            info.d1 = length(info.crop.x);
            info.d2 = length(info.crop.y);
            % info.skipFrame = 65537; % Skip frame 65537 (for old recordings which
                       % have a bug that repeats a frame after 65537 (2^16) frames)
            info.Skipframe = -1; % Don't skip any frames
            info.AlignMethod = do.motion.alignMethod;
            
            % START NORMCORR REGISTRATION % % % % %
            tics.normcorr = tic;
            simonalign3(pno);
            timedDatai(1).normcorrT = toc(tics.normcorr);
            
            files.nameNormcorr = cell(nSlices,1);
            if info.bsplit
                for j = 1:nSlices
                    files.pathNormcorr = pno;
                    files.nameNormcorr{j} = [fnUpdate sprintf('_Split%d_normcorr', j)];
                end
            else
                files.pathNormcorr = pno;
                files.nameNormcorr{j}= [fnUpdate '_normcorr'];
            end
            
        else % No NoRMCorre
            nSlices = 1;
            files.pathNormcorr = pn;
            files.nameNormcorr = {fnUpdate};
            timedDatai(1).normcorrT = NaN;
        end
        
        %% BACKGROUND SUBTRACTION % % % %
        % Basically only necessary for 1P & miniscope data, to remove extreme
        % vignetting, out of focus fluorescence, blur
        if do.BackgroundSubtract
            files.nameBackgroundSub = cell(nSlices, 1);
            fprintf('Doing background subtraction...\n')
            tics.backSub = tic;
            
            for j = 1:nSlices
                files.nameBackgroundSub{j} = [files.nameNormcorr{j} '_Sub'];
                BackgroundSubtractSbx([files.pathNormcorr files.nameNormcorr{j}],...
                                      files.nameBackgroundSub{j}, do.BackgroundSub.FilterRadius,...
                                      do.BackgroundSub.filterMethod, do.BackgroundSub.Shifter, do.BackgroundSub.SmoothDim, do.BackgroundSub.SmoothSe, do.BackgroundSub.plotter)
            end
            timedDatai(1).backSubT = toc(tics.backSub);
            % Replacing the normcorr names! So next analysis takes correct file
            files.nameNormcorr = files.nameBackgroundSub;
        else % No background subtraction
            timedDatai(1).backSubT = NaN;
        end
        
        
        %% Eye file file conversion % % % % 
        if do.EyeConvert
            files.nameEye = strsplit(fn, {'_normcorr','_Split','_split','_copy','.sbx'});
            files.nameEye = [files.nameEye{1}, '_eye.mat'];
            % Check presence eye files.
            if exist([pn files.nameEye], 'files.')
                EyeData2Avi(pn, files.nameEye, do.EyeDeFlicker, do.EyeInvert)
            else
                fprintf('No Eye file found!\n')
            end
        end
        
        
        %% TRANSPOSE DATASETS % % % % %
        files.nameTrans = cell(nSlices,1);
        fprintf('Starting Transposing the %d splits!\n', nSlices)
        
        tics.transpose = tic;
        for j = 1:nSlices
            StackTranspose([files.pathNormcorr files.nameNormcorr{j} '.sbx'], pno)
            files.nameTrans{j} = [pno files.nameNormcorr{j} '_Trans.dat'];
        end
        timedDatai(1).transposeT = toc(tics.transpose);
        
        
        %% DECIMATE DATASETS % % % % % %  
        freqDec = 1;
        fprintf('Decimating the %d datasets to %.1fHz!\n', nSlices, freqDec)
        files.nameDecTrans = cell(nSlices,1);
        tics.decimate = tic;
        for j = 1:nSlices
            DecimateTrans(files.nameTrans{j}, freqDec)
            files.nameDecTrans{j} = [pno files.nameNormcorr{j} '_DecTrans.dat'];
        end
        timedDatai(1).decimatT = toc(tics.decimate);
        
        
        %% SPECTRAL % % % % % % % %
        fprintf('Doing spectral analysis on the transposed files\n')
        files.nameSpectral = cell(nSlices, 1);
        
        tics.spectral = tic;
        for j = 1:nSlices
            spectral(files.nameDecTrans{j});
            files.nameSpectral{j} = [pno files.nameNormcorr{j} '_SPSIG.mat'];
        end
        timedDatai(1).spectralT = toc(tics.spectral);
        
        
        %% FLUORESCENCE PROJECTION IMAGES % % % % % 
        clearvars info
        clearvars -global info
        
        tics.fluorescenceImg = tic;
        for j = 1:nSlices
            FluorescenceImgSbx([pno files.nameNormcorr{j}]);
        end
        timedDatai(1).fluorescenceT = toc(tics.fluorescenceImg);
        
        
        %% Automatic ROI creation
        if do.getROIs
            tics.getRois = tic;
            for j = 1:nSlices
                getSpectrois(files.nameSpectral{j}, spar)
            end
            timedDatai(1).getRoisT = toc(tics.getRois);
            load(files.nameSpectral{j},'PP')
            nRois = nRois + PP.Cnt;
        else
            timedDatai(1).getRoisT = NaN;
        end
        
    catch ME
        if do.timed
            timedDatai(1).error = ME;
            timedDatai(1).memory = memory;
            FillTimedData(timedDatai, do.timedFile, pn, fn, nRois, nSlices, do.motion.alignMethod)
            fprintf('Saved timed data for file %d\n', i)
        end
        rethrow(ME)
    end
    %% fill in the timedData tracker investigation file
    if do.timed
        FillTimedData(timedDatai, do.timedFile, pn, fn, nRois, nSlices, do.motion.alignMethod)
        fprintf('Saved timed data for file %d\n', i)
    end
end


