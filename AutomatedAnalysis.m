% Doing the normcorr motion correction, transposing, spectral analysis and
% automated ROI creation for all splits, for multiple files in one script
% Optional steps are motion correction, background subtraction, ROI getting
%
% Leander de Kraker

% CROPPING SETTINGS % % % % % (for motion correction!)
% Most important part of cropping is removing the lines (on the left
% and right side) which have value 653357 and will thus overpower all 
% relevant data during motion correction.
% x = 102:774; % HORIZONTAL CROP
% y = 25:512;  % VERTICAL CROP

% Run motion normcorre correction? % % % % % % % % % % % % %
doNoRMCorre = false; % (if false, cropping settings aren't used)

% BACKGROUND SUBTRACTION parameters % % % % % Necessary for 1 photon data
doBackgroundSubtract = true;
filterRadius = 55; % The radius of the filter for background subtraction
filterMethod = 'Circular Average'; % Options Gaussian Average | Circular Average
shifter = 9000; % The data values has to be increased to prevent underexposure/ over correction
% 1D gaussian smoothing to slightly decrease vertical or horizontal banding noise. 
smoothDim = 0; % options: 0 (no smoothing), 1 (horizontal smoothing), 2 (vertical smoothing)
smoothSe  = []; % size for the smoothing
plotter = true; % plot last frame


% Also run getSpectrois to detect the ROIs? % % % % % % % %
getROIs = true;
if getROIs
    % Arm the spar: Spectral roi finder PARameters
    global spar
    spar = Spectroiparm();
end

% Save timing info? % % % % % % % % % % % % % % % % % % % % 
timed = false;
if timed
    % File in which timing info is
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

% Load as many files as user wants
filenames = {};
filepaths = {};
selecting = true;
i = 0;
while selecting
    i = i + 1;
    [filenames{i}, filepaths{i}] = uigetfile('*sbx', sprintf('file %d',i));
    
    if filenames{i} == 0 % Cancel is pressed probably: stop with selecting
        filenames(i) = [];
        filepaths(i) = [];
        selecting = false;
    end
end
filenames = filenames';
filepaths = filepaths';
nfiles = length(filenames);
clearvars selecting

i = 1;
%% process the files
for i = 1:nfiles
    fn = filenames{i};
    pn = filepaths{i};
    
    pnfn = [pn fn(1:end-4)]; % pathname with filename  
    % prepare to load a different sbx file: thoroughly delete info variable
    clearvars info 
    clearvars -global info
    
    % An example frame of the recording (creates new global info variable)
    Img = sbxread(pnfn, 0,1);
    
    global info
    
    % set the name of the recording
    if ~isfield(info, 'strfp')
        info.strfp = pnfn;
    elseif ~strcmp(info.strfp, pnfn) % inconsitent names, maybe it'll all go wrong
        info.strfp
        pnfn
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
    % Place normcorre folder in Matlab path BELOW SpectralSegmentation
    
    % Show the data before committing the analysis
    Img = sbxread(pnfn, 0,100*nSlices);
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
    colormap(cmapL([1     1    1;...
                    1     1    0.75;...
                    0.75  1    0.5;...
                    0.25  0.75 0.5;...
                    0     0.25 0.5;...
                    0     0    0], 256))
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
        
        % START NORMCORR REGISTRATION % % % % %
        normcorrTic = tic;
        simonalign3;
        normcorrtoc = toc(normcorrTic);
        
        filenameNormcorr = cell(nSlices,1);
        if info.bsplit
            for j = 1:nSlices
                filenameNormcorr{j} = [pnfn sprintf('_Split%d_normcorr', j)];
            end
        else
            filenameNormcorr{j}= [pnfn sprintf('_normcorr', j)];
        end
        
    else % No NoRMCorre
        nSlices = 1;
        filenameNormcorr = {pnfn};
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
            BackgroundSubtractSbx(filenameNormcorr{j},...
                                  filenameBackgroundSub{j}, filterRadius,...
                                  filterMethod, shifter, smoothDim, smoothSe, plotter)
        end
        backSubtoc = toc;
        % Replacing the normcorr names! So next analysis takes correct file
        filenameNormcorr = filenameBackgroundSub;
    else % No background subtraction
        backSubtoc = NaN;
    end
    
    
    %% TRANSPOSE DATASETS % % % % %
    filenameTrans = cell(nSlices,1);
    fprintf('Starting Transposing the %d splits!\n', nSlices)
    
    transposeTic = tic;
    for j = 1:nSlices
        StackTranspose([filenameNormcorr{j} '.sbx'], pn)
        filenameTrans{j} = [filenameNormcorr{j} '_Trans.dat'];
    end
    transposetoc = toc(transposeTic);
    
    
    %% DECIMATE DATASETS % % % % % %  
    fprintf('Decimating the %d datasets to 1Hz!\n', nSlices)
    filenameDecTrans = cell(nSlices,1);
    decTic = tic;
    for j = 1:nSlices
        DecimateTrans(filenameTrans{j})
        filenameDecTrans{j} = [filenameNormcorr{j} '_DecTrans.dat'];
    end
    dectoc = toc(decTic);
    
    
    %% SPECTRAL % % % % % % % %
    fprintf('Doing spectral analysis on the transposed files\n')
    filenameSpectral = cell(nSlices, 1);
    
    spectralTic = tic;
    for j = 1:nSlices
        spectral(filenameDecTrans{j});
        filenameSpectral{j} = [filenameNormcorr{j} '_SPSIG.mat'];
    end
    spectraltoc = toc(spectralTic);
    
    
    %% FLUORESCENCE PROJECTION IMAGES % % % % % 
    clearvars info
    clearvars -global info
    
    fluorescenceTic = tic;
    for j = 1:nSlices
        FluorescenceImgSbx(filenameNormcorr{j});
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
        fileInfo = dir([pn fn]);
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
        save(timedFile, 'timedData')
        fprintf('saved timeddata for file %d\n', i)
    end
end


