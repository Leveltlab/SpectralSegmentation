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
        files.name = filesTxtName;
        files.path = filesTxtFolder;
    elseif ismember()
        
    end
end
if isstruct(files)
    
end

if isempty(do) || isempty(spar) % if either spar or do is missing, start app
    SetAnalysisSettings(files, do, spar)
elseif (isstruct(do) && do.saveLocation) && (isstruct(files) && isfield(files.saveFolder))
    % If it should save in different folders, but those folders are not defined
    SetAnalysisSettings(files, do, spar)
elseif isempty(files) % Select files via pop up if files are not given
    [files, nfiles] = SelectFiles();
end

global info


%% 
% Set the i if you want to run a specific section of the analysis for a
% specific file (after things crash or are unsatisfactory for example)
i = 1;
[files, nSlices] = FilenamesAutomatedAnalysis(files, do);
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
        if do.LineShift && do.LineShiftValue~=0
            tics.lineShift = tic;
            fprintf('Shifting lines & replacing high value of file %s...\n', fn)
            ShiftLinesSbx({pn fn}, do.LineShiftTrans, do.LineShiftValue, pno)
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
        
        [info, nSlices] = UpdateSbxInfo(info);
        
        
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
            clim(cLimits)
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
            
            fnNormcorr = cell(nSlices,1);
            if info.bsplit
                for j = 1:nSlices
                    pnNormcorr = pno;
                    fnNormcorr{j} = [fnUpdate sprintf('_Split%d_normcorr', j)];
                end
            else
                pnNormcorr = pno;
                fnNormcorr{j}= [fnUpdate '_normcorr'];
            end
        else % No NoRMCorre
            nSlices = 1;
            pnNormcorr = pn;
            fnNormcorr = {fnUpdate};
            timedDatai(1).normcorrT = NaN;
        end
        
        %% BACKGROUND SUBTRACTION % % % %
        % Basically only necessary for 1P & miniscope data, to remove extreme
        % vignetting, out of focus fluorescence, blur
        if do.doBackgroundSubtract
            fnBackgroundSub = cell(nSlices, 1);
            fprintf('Doing background subtraction...\n')
            tics.backSub = tic;
            
            for j = 1:nSlices
                fnBackgroundSub{j} = [fnNormcorr{j} '_Sub'];
                BackgroundSubtractSbx([pnNormcorr fnNormcorr{j}],...
                                      fnBackgroundSub{j}, do.BackgroundSubtract.FilterRadius,...
                                      do.BackgroundSubtract.filterMethod, do.BackgroundSubtract.Shifter, do.BackgroundSubtract.SmoothDim, do.BackgroundSubtract.SmoothSe, do.BackgroundSubtract.plotter)
            end
            timedDatai(1).backSubT = toc(tics.backSub);
            % Replacing the normcorr names! So next analysis takes correct file
            fnNormcorr = fnBackgroundSub;
        else % No background subtraction
            timedDatai(1).backSubT = NaN;
        end
        
        
        %% Eye file file conversion % % % % 
        if do.doEyeConvert
            fnEye = strsplit(fn, {'_normcorr','_Split','_split','_copy','.sbx'});
            fnEye = [fnEye{1}, '_eye.mat'];
            % Check presence eye files.
            if exist([pn fnEye], 'files.')
                EyeData2Avi(pn, fnEye, do.EyeDeFlicker, do.EyeInvert)
            else
                fprintf('No Eye file found!\n')
            end
        end
        
        
        %% TRANSPOSE DATASETS % % % % %
        fnTrans = cell(nSlices,1);
        fprintf('Starting Transposing the %d splits!\n', nSlices)
        
        tics.transpose = tic;
        for j = 1:nSlices
            StackTranspose([pnNormcorr fnNormcorr{j} '.sbx'], pno)
            fnTrans{j} = [pno fnNormcorr{j} '_Trans.dat'];
        end
        timedDatai(1).transposeT = toc(tics.transpose);
        
        
        %% DECIMATE DATASETS % % % % % %  
        freqDec = 1;
        fprintf('Decimating the %d datasets to %.1fHz!\n', nSlices, freqDec)
        fnDecTrans = cell(nSlices,1);
        tics.decimate = tic;
        for j = 1:nSlices
            DecimateTrans(fnTrans{j}, freqDec)
            fnDecTrans{j} = [pno fnNormcorr{j} '_DecTrans.dat'];
        end
        timedDatai(1).decimatT = toc(tics.decimate);
        
        
        %% SPECTRAL % % % % % % % %
        fprintf('Doing spectral analysis on the transposed files\n')
        fnSpectral = cell(nSlices, 1);
        
        tics.spectral = tic;
        for j = 1:nSlices
            spectral(fnDecTrans{j});
            fnSpectral{j} = [pno fnNormcorr{j} '_SPSIG.mat'];
        end
        timedDatai(1).spectralT = toc(tics.spectral);
        
        
        %% FLUORESCENCE PROJECTION IMAGES % % % % % 
        clearvars info
        clearvars -global info
        
        tics.fluorescenceImg = tic;
        for j = 1:nSlices
            FluorescenceImgSbx([pno fnNormcorr{j}]);
        end
        timedDatai(1).fluorescenceT = toc(tics.fluorescenceImg);
        
        
        %% Automatic ROI creation
        if do.doGetROIs
            tics.getRois = tic;
            for j = 1:nSlices
                getSpectrois(fnSpectral{j}, spar)
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
        warning('\nCrashed in iteration %d\n%s\n\n', i, ME.message)
        continue
        % rethrow(ME)
    end
    %% fill in the timedData tracker investigation file
    if do.doTimedFile
        FillTimedData(timedDatai, do.timedFile, pn, fn, nRois, nSlices, do.doNoRMCorre.alignMethod)
        fprintf('Saved timed data for file %d\n', i)
    end
end


