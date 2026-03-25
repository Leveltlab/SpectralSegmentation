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
    spar = [];
end

% Set the files to analyse
if ischar(files) || isstring(files) % if text is given
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

if isempty(do) || isempty(spar)
    SetAnalysisSettings(files, do, spar)
elseif isempty(files) % Select files via pop up if nothing is given
    [files, nfiles] = SelectFiles();
end

global info


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
            if do.doNoRMCorre
                hold on
                rectangle('position',[x(1), y(1), x(end)-x(1), y(end)-y(1)],...
                          'edgecolor',[1 1 0.75],'linewidth',2.5)
            end
            title(sprintf('slice %d',j))
        end
        linkaxes(handles, 'xy')
        colormap(cmapL('greenFancy', 256))
        drawnow
        
        if do.doNoRMCorre
            % critical info for NoRMCorre
            info.crop.x = do.doNoRMCorre.crop.x; 
            info.crop.y = do.doNoRMCorre.crop.y;
            info.d1 = length(info.crop.x);
            info.d2 = length(info.crop.y);
            % info.skipFrame = 65537; % Skip frame 65537 (for old recordings which
                       % have a bug that repeats a frame after 65537 (2^16) frames)
            info.Skipframe = -1; % Don't skip any frames
            info.AlignMethod = do.doNoRMCorre.alignMethod;
            
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
        if do.doBackgroundSubtract
            files.nameBackgroundSub = cell(nSlices, 1);
            fprintf('Doing background subtraction...\n')
            tics.backSub = tic;
            
            for j = 1:nSlices
                files.nameBackgroundSub{j} = [files.nameNormcorr{j} '_Sub'];
                BackgroundSubtractSbx([files.pathNormcorr files.nameNormcorr{j}],...
                                      files.nameBackgroundSub{j}, do.BackgroundSubtract.FilterRadius,...
                                      do.BackgroundSubtract.filterMethod, do.BackgroundSubtract.Shifter, do.BackgroundSubtract.SmoothDim, do.BackgroundSubtract.SmoothSe, do.BackgroundSubtract.plotter)
            end
            timedDatai(1).backSubT = toc(tics.backSub);
            % Replacing the normcorr names! So next analysis takes correct file
            files.nameNormcorr = files.nameBackgroundSub;
        else % No background subtraction
            timedDatai(1).backSubT = NaN;
        end
        
        
        %% Eye file file conversion % % % % 
        if do.doEyeConvert
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
        if do.doGetROIs
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
        if do.doTimedFile
            timedDatai(1).error = ME;
            timedDatai(1).memory = memory;
            FillTimedData(timedDatai, do.timedFile, pn, fn, nRois, nSlices, do.doNoRMCorre.alignMethod)
            fprintf('Saved timed data for file %d\n', i)
        end
        rethrow(ME)
    end
    %% fill in the timedData tracker investigation file
    if do.doTimedFile
        FillTimedData(timedDatai, do.timedFile, pn, fn, nRois, nSlices, do.doNoRMCorre.alignMethod)
        fprintf('Saved timed data for file %d\n', i)
    end
end


