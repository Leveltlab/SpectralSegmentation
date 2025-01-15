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

% Load as many files as user wants
filenames = {};
filepaths = {};
selecting = true;
i = 0;
while selecting
    i = i + 1;
    if i > 1
        prompt = sprintf('file %d. Previous file: %s',i, filenames{i-1});
    else
        prompt = sprintf('file %d. Press cancel when done.', i);
    end
    
    [filenames{i}, filepaths{i}] = uigetfile('*sbx', prompt, 'MultiSelect', 'on');
    
    if ~iscell(filenames{i}) & filenames{i} == 0 % Cancel is pressed probably: stop with selecting
        filenames(i) = [];
        filepaths(i) = [];
        selecting = false;
    elseif iscell(filenames{i}) % multiple files have been selected at once
        nselected = length(filenames{i});
        filenames(i:i+nselected-1) = filenames{i};
        filepaths(i:i+nselected-1) = filepaths(i);
        i = i + nselected - 1;
    end
end
filenames = filenames';
filepaths = filepaths';
nfiles = length(filenames);

for i = 1:nfiles
    filenames{i} = strsplit(filenames{i}, '.sbx');
    filenames{i} = filenames{i}{1};
end
clearvars selecting i


global info

%% Settings via script

% CROPPING SETTINGS % % % % % (for motion correction!)
% Most important part of cropping is removing the lines (on the left
% and right side) which have value 653357 and will thus overpower all 
% relevant data during motion correction.
x = 3:793; % HORIZONTAL CROP
y = 5:510;  % VERTICAL CROP

% Do line shifting correction? % % % % % % % % %
doLineShift = false;
doLineShiftTrans = false; % false=correct for horizontal lines, true=correct for vertical lines
lineShift = 0;

% Run motion normcorre correction? % % % % % % % % % % % % %
doNoRMCorre = true; % (if false, cropping settings aren't used)
alignMethod = 'Rigid'; % options: 'Nonrigid' | 'Rigid'

% Convert *_eye.mat file into .mp4 file?
doEyeConvert = true;
doEyeInvert = true; % Invert black-white?
doEyeDeFlicker = false; % Diminish flickering brightness between subsequent frames?

filepathsOutput = filepaths;

% BACKGROUND SUBTRACTION parameters % % % % % Necessary for 1 photon data
doBackgroundSubtract = false;
filterRadius = 55; % The radius of the filter for background subtraction
filterMethod = 'Circular Average'; % Options Gaussian Average | Circular Average
shifter = 9000; % The data values has to be increased to prevent underexposure/ over correction
% 1D gaussian smoothing to slightly decrease vertical or horizontal banding noise. 
smoothDim = 1; % options: 0 (no smoothing), 1 (horizontal smoothing), 2 (vertical smoothing)
smoothSe  = 1.75; % size for the smoothing
plotter = true; % plot last frame

% Also run getSpectrois to detect the ROIs? % % % % % % % %
getROIs = false;
if getROIs
    % Arm the spar: Spectral roi finder PARameters
    global spar
    spar = Spectroiparm();
end

% Save timing info? % % % % % % % % % % % % % % % % % % % % 
timed = false;
if timed
    % File in which timing info is stored
    timedFile = 'D:\spectralTimed.mat';
    if ~isfile(timedFile)
        warning('%s does not exist!! please say which file to add the timing data to!', timedFile)
        return
% %         or run following code to create a new timedFile there: 
        timedData = struct('computerName',{},'fileSize',[],'filePath', [], 'fileName',{},...
                     'nSplits', [], 'normcorrT',[],'transposeT',[],'decimatT',[],...
                     'spectralT',[],'fluorescenceT',[], 'backSubT', [], 'getRoisT', [],...
                     'nRois', [], 'date', {});
        save(timedFile, 'timedData')
    end
end


%% Get settings via input dialogs
dlgdims = [1 85];
clearvars info; clearvars -global info

% Do line shift? %
doLineShift = questdlg('Do line shifting? to correct bi-directional scanning misalignment', 'line shift?',...
                       'yes', 'no', 'no');
if strcmp(doLineShift, 'yes')
    doLineShift = true;
    doLineShiftTrans = questdlg('correct for horizontal or vertical lines?', 'orientation of lines?',...
                                'horizontal', 'vertical', 'horizontal');
    if strcmp(doLineShiftTrans, 'horizontal')
        doLineShiftTrans = false;
    else
        doLineShiftTrans = true;
    end
elseif strcmp(doLineShift, 'no')
    doLineShift = false;
else
    fprintf('Canceled selecting settings\n')
    return
end
if doLineShift
    im = sbxread([filepaths{1}, filenames{1}], 5, 100);
    lineShift = ShiftLinesCheck(im, doLineShiftTrans, true);
%     lineShift = inputdlg('shift by how much');
end


% Do NoRMCorre? %
%
doNoRMCorre = questdlg('do NoRMCorre motion correction?', 'run NoRMCorre?', 'yes', 'no', 'yes');
if strcmp(doNoRMCorre, 'yes')
    doNoRMCorre = true;
elseif strcmp(doNoRMCorre, 'no')
    doNoRMCorre = false;
else
    fprintf('canceled selecting settings\n')
    return
end
if doNoRMCorre
    % Get crop for NoRMCorre
    if exist('filepaths', 'var')
        load([filepaths{1} filenames{1} '.mat'], 'info');
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
    x = str2num(answer{1});
    y = str2num(answer{2});
    
    % Set rigid or non rigid
    alignMethod = questdlg('Which registration method to use? (rigid = default)',...
                           'NoRMCorre registration method',...
                           'Rigid','Nonrigid','Rigid');
end


% Save in original file location ? % 
doSaveLocation = questdlg('Save new files in original folder? (yes = default)',...
                          'Output location',...
                          'yes', 'no (select different folders)', 'yes');
filepathsOutput = filepaths;
if strcmp(doSaveLocation, 'yes')
    doSaveLocation = false;
elseif strcmp(doSaveLocation, 'no (select different folders)') % select folder for every file
    doSaveLocation = true;
    outputpathi = 'lalala';
    i = 1;
    while i <= nfiles && ischar(outputpathi)
        outputpathi = uigetdir(cd, sprintf('output path for file %d: %s',i, filenames{i}));
        if ischar(outputpathi)
            filepathsOutput{i} = [outputpathi '\'];
        end
        i = i + 1;
    end
else
    fprintf('Canceling selecting settings\n')
    return
end


% Do eye file format conversion? (keeps the original)
doEyeConvert = questdlg('Convert "*_eye.mat" file into "*_eye.mp4" (keeps original)',...
                        'Create videofile for eye tracking?',...
                        'yes','no','yes');
if strcmp(doEyeConvert, 'yes')
    doEyeConvert = true;
    doEyeInvert = questdlg('Invert black-white in eye video file?',...
                               'Eye video file settings',...
                               'yes', 'no', 'yes');
    doEyeDeFlicker = questdlg('Deflicker the eye video?',...
                              'Eye video file settings',...
                               'yes', 'no', 'yes');
    if strcmp(doEyeInvert, 'yes')
        doEyeInvert = true;
    else
        doEyeInvert = false;
    end
    if strcmp(doEyeDeFlicker, 'yes')
        doEyeDeFlicker = true;
    else
        doEyeDeFlicker = false;
    end
elseif strcmp(doEyeConvert, 'no')
    doEyeConvert = false;
else
    fprintf('Canceling selecting settings\n')
    return
end
% Do background subtraction? %
doBackgroundSubtract = questdlg('Do background subtraction (for 1-photon data)',...
                                'background subtraction settings',...
                                'yes','no','no');
if strcmp(doBackgroundSubtract, 'yes');    doBackgroundSubtract = true;
elseif strcmp(doBackgroundSubtract, 'no'); doBackgroundSubtract = false;
else; fprintf('Cancelling selecting settings\n'); return
end
% Get the specific background subtraction settings
if doBackgroundSubtract
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
    filterRadius = str2double(answer{1});
    shifter = str2double(answer{2});
    if strcmpi(answer{3}, 'yes')
        plotter = true;
    else
        plotter = false;
    end
    filterMethod = questdlg('filtering method for background image estimation', 'background subtract settings',...
                            'Circular Average (advised)',...
                            'Gaussian Average',...
                            'Circular Average (advised)');
    if strcmp(filterMethod, 'Circular Average (advised)')
        filterMethod = 'Circular Average';
    elseif strcmp(filterMethod, 'Gaussian Average')
        % That's the correct string already: do nothing
    else
        fprintf('Stopping setting selection, no valid method for background subtraction\n')
        return
    end
    smoothDim = questdlg('1D gaussian smoothing to diminish (vertical or horizontal banding) noise?',...
                         'background subtract setings',...
                         'no', 'horizontal smoothing', 'vertical smoothing', 'vertical smoothing');
    if strcmp(smoothDim, 'no')
        smoothDim = 0;
        smoothSe = [];
    elseif strcmp(smoothDim, 'horizontal smoothing')
        smoothDim = 1;
    elseif strcmp(smoothDim, 'vertical smoothing')
        smoothDim = 2;
    else % canceled
        fprintf('Canceling selecting settings\n')
        return
    end
    if smoothDim > 0
        answer = inputdlg('size of smoothing (standard deviation of gaussian)',...
                          'background subtract settings',...
                          dlgdims, {'1.7'});
        if isempty(answer)
            fprintf('Canceled selecting settings at final dlg\n')
            return
        end
        smoothSe = str2num(answer{1});
    end
end

% get ROIs? %
getROIs = questdlg('automatic ROI detection?', 'get ROIs?', 'yes', 'no', 'yes');
if strcmp(getROIs, 'yes')
    getROIs = true;
    global spar
    spar = Spectroiparm();
elseif strcmp(getROIs, 'no')
    getROIs = false;
else
    fprintf('Cancelling selecting settings\n')
    return
end

% Timing info % 
timed = questdlg('Save processing time/log of the analysis?', 'time/log analysis?',...
    'yes', 'no', 'no');
if strcmp(timed, 'yes')
    timed = true;
elseif strcmp(timed, 'no')
    timed = false;
else
    timed = false;
end
if timed
    % File in which timing info is stored
    timedFile = 'D:\2Pdata\spectralTimed.mat'; % Default timedFile name! % % % 
    timedFileLinePos = dbstack;
    if ~isfile(timedFile) % Creating new file for timing data
        warning('%s does not exist! please say which file to add the timing data to!', timedFile)
        [timedFileName, timedFilePath] = uiputfile('*.mat', 'Where to save timing file', 'spectralTimed.mat');
        timedFile = [timedFilePath timedFileName];
        if timedFileName==0 % User didn't select a file so forget about saving timed data
            timed = false;
        elseif ~isfile(timedFile) % Save new requested file
            timedFile = [timedFilePath, timedFile];
            if isscalar(timedFileLinePos) && isscalar(timedFileLinePos.line)
                fprintf('\nPlease change default timedFile name (line %d) to the new file for future use:\n%s\n',...
                    timedFileLinePos.line-1, timedFile)
            else
                fprintf('\nPlease change default timedFile name to the new file for future use:\n%s\n', timedFile)
            end
            timedData = struct('computerName',{},'fileSize',[],'filePath', [], 'fileName',{},...
                         'nSplits', [], 'normcorrT',[],'transposeT',[],'decimatT',[],...
                         'spectralT',[],'fluorescenceT',[], 'backSubT', [], 'getRoisT', [],...
                         'nRois', [], 'date', {}, 'alignMethod', {});
            save(timedFile, 'timedData')
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
    fn = filenames{i};
    pn = filepaths{i};
    pno = filepathsOutput{i};
    % prepare to load a different sbx file: thoroughly delete info variable
    clearvars info
    clearvars -global info
    
    %% Shift lines, and replace erroneously high values (>65530)
    if doLineShift
        fprintf('Shifting lines & replacing high value of file %s...\n', fn)
        ShiftLinesSbx({pn fn}, doLineShiftTrans, lineShift, pno)
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
    cLimits = prctile(Img(:), [0.02 96]);
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
        if doNoRMCorre
            hold on
            rectangle('position',[x(1), y(1), x(end)-x(1), y(end)-y(1)],...
                      'edgecolor',[1 1 0.75],'linewidth',2.5)
        end
        title(sprintf('slice %d',j))
    end
    linkaxes(handles, 'xy')
    colormap(cmapL('greenFancy', 256))
    drawnow
    
    if doNoRMCorre
        % critical info for NoRMCorre
        info.crop.x = x; 
        info.crop.y = y;
        info.d1 = length(info.crop.x);
        info.d2 = length(info.crop.y);
        % info.skipFrame = 65537; % Skip frame 65537 (for old recordings which
                   % have a bug that repeats a frame after 65537 (2^16) frames)
        info.Skipframe = -1; % Don't skip any frames
        info.AlignMethod = alignMethod;
        
        % START NORMCORR REGISTRATION % % % % %
        normcorrTic = tic;
        simonalign3(pno);
        normcorrtoc = toc(normcorrTic);
        
        filenameNormcorr = cell(nSlices,1);
        if info.bsplit
            for j = 1:nSlices
                filepathNormcorr = pno;
                filenameNormcorr{j} = [fnUpdate sprintf('_Split%d_normcorr', j)];
            end
        else
            filepathNormcorr = pno;
            filenameNormcorr{j}= [fnUpdate '_normcorr'];
        end
        
    else % No NoRMCorre
        nSlices = 1;
        filepathNormcorr = pn;
        filenameNormcorr = {fnUpdate};
        normcorrtoc = NaN;
    end
    
    %% BACKGROUND SUBTRACTION % % % %
    % Basically only necessary for 1P & miniscope data, to remove extreme
    % vignetting, out of focus fluorescence, blur
    if doBackgroundSubtract
        filenameBackgroundSub = cell(nSlices, 1);
        fprintf('Doing background subtraction...\n')
        backSubtic = tic;
        
        for j = 1:nSlices
            filenameBackgroundSub{j} = [filenameNormcorr{j} '_Sub'];
            BackgroundSubtractSbx([filepathNormcorr filenameNormcorr{j}],...
                                  filenameBackgroundSub{j}, filterRadius,...
                                  filterMethod, shifter, smoothDim, smoothSe, plotter)
        end
        backSubtoc = toc(backSubtic);
        % Replacing the normcorr names! So next analysis takes correct file
        filenameNormcorr = filenameBackgroundSub;
    else % No background subtraction
        backSubtoc = NaN;
    end
    
    
    %% Eye file file conversion % % % % 
    if doEyeConvert
        filenameEye = strsplit(fn, {'_normcorr','_Split','_split','_copy','.sbx'});
        filenameEye = [filenameEye{1}, '_eye.mat'];
        % Check presence eye file
        if exist([pn filenameEye], 'file')
            EyeData2Avi(pn, filenameEye, doEyeDeFlicker, doEyeInvert)
        else
            fprintf('No Eye file found!\n')
        end
    end
    
    
    %% TRANSPOSE DATASETS % % % % %
    filenameTrans = cell(nSlices,1);
    fprintf('Starting Transposing the %d splits!\n', nSlices)
    
    transposeTic = tic;
    for j = 1:nSlices
        StackTranspose([filepathNormcorr filenameNormcorr{j} '.sbx'], pno)
        filenameTrans{j} = [pno filenameNormcorr{j} '_Trans.dat'];
    end
    transposetoc = toc(transposeTic);
    
    
    %% DECIMATE DATASETS % % % % % %  
    freqDec = 1;
    fprintf('Decimating the %d datasets to %.1fHz!\n', nSlices, freqDec)
    filenameDecTrans = cell(nSlices,1);
    decTic = tic;
    for j = 1:nSlices
        DecimateTrans(filenameTrans{j}, freqDec)
        filenameDecTrans{j} = [pno filenameNormcorr{j} '_DecTrans.dat'];
    end
    dectoc = toc(decTic);
    
    
    %% SPECTRAL % % % % % % % %
    fprintf('Doing spectral analysis on the transposed files\n')
    filenameSpectral = cell(nSlices, 1);
    
    spectralTic = tic;
    for j = 1:nSlices
        spectral(filenameDecTrans{j});
        filenameSpectral{j} = [pno filenameNormcorr{j} '_SPSIG.mat'];
    end
    spectraltoc = toc(spectralTic);
    
    
    %% FLUORESCENCE PROJECTION IMAGES % % % % % 
    clearvars info
    clearvars -global info
    
    fluorescenceTic = tic;
    for j = 1:nSlices
        FluorescenceImgSbx([pno filenameNormcorr{j}]);
    end
    fluorescencetoc = toc(fluorescenceTic);
    
    
    %% Automatic ROI creation
    nRois = 0;
    if getROIs
        getRoisTic = tic;
        for j = 1:nSlices
            getSpectrois(filenameSpectral{j}, spar)
        end
        getRoistoc = toc(getRoisTic);
        load(filenameSpectral{j},'PP')
        nRois = nRois + PP.Cnt;
    else
        getRoistoc = NaN;
    end
    
    %% fill in the timedData tracker investigation file
    if timed
        load(timedFile)
        fileInfo = dir([pn fn '.sbx']);
        timedData(end+1).fileSize = fileInfo.bytes * 10^-9;
        timedData(end).normcorrT  = normcorrtoc;
        timedData(end).backSubT   = backSubtoc;
        timedData(end).transposeT = transposetoc;
        timedData(end).decimatT   = dectoc;
        timedData(end).spectralT  = spectraltoc;
        timedData(end).fluorescenceT  = fluorescencetoc;
        timedData(end).getRoisT   = getRoistoc;
        timedData(end).nRois      = nRois;
        timedData(end).nSplits    = nSlices;
        timedData(end).computerName = getenv('COMPUTERNAME');
        timedData(end).fileName   = fn;
        timedData(end).filePath   = pn;
        timedData(end).date       = datetime;
        timedData(end).alignMethod= alignMethod;
        save(timedFile, 'timedData')
        fprintf('saved timeddata for file %d\n\n', i)
    end
end


