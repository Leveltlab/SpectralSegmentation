% Doing the normcorr motion correction, transposing, spectral analysis and
% automated ROI creation for all splits, for multiple files in one script
% 

% CROPPING SETTINGS % % % % %  % % %% % 
% Most important part of cropping is removing the lines (on the left
% and right side) which have value 653357 and will thus overpower all 
% relevant data during motion correction.
x = 102:774; % HORIZONTAL CROP
y = 25:512;  % VERTICAL CROP
% x = 50:720;
% y = 10:512;
% for GPU recording
% x = 9:792;
% y = 9:507;

% Also run getspectrois to retrive the ROIs?
getROIs = true;
if getROIs
    % Arm the spar
    global spar
    spar = Spectroiparm();
end

% Save timing info?
timed = true;
if timed
    % File in which timing info is
    timedFile = '\\mvp2nin\Projects\PLInterneurons\spectralTimed2.mat';
    if ~isfile(timedFile)
        warning('%s does not exist!! please say which file to add the timing data to!')
        return
% %         or run following code to create a new timedFile there: 
        timedData = struct('n',0,'computerName',{},'fileSize',[],'fileName',{},...
                     'nSplits', [], 'normcorrT',[],'transposeT',[],'decimatT',[],...
                     'spectralT',[],'backgrndT',[], 'getRoisT', [], 'nRois', []);
        timedData(1).n = 0;
        save(timedFile, 'timedData')
    end
end

filenames = {};
filepaths = {};
selecting = true; % true as long as files are being selected
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
        try
            Splitnm = info.otparam(3);
        catch
%             % Maybe the info file is super bad, ask how many slices there are
%             if ~isfield(info, 'Slices')
%                 answ = inputdlg('How many Slices?', 'Slice info!!...', [1 40],  {'1'});
%                 Splitnm = str2double(answ{1});
%             end

            % The info file could have no slices because only 1 depth is imaged
            Splitnm = 1;
        end
        info.Slices = Splitnm;
    end
    if size(Img,1) < 3
        bVers = 1;
        Img = permute(Img, [2 3 1]);

        if ~isfield(info, 'bVers')
            info.bVers = bVers;
        end
    end

    nSlices = info.Slices;
    info.bsplit = 1; % MULTIPLE SPLITS OR NOT? true = yes multiple splits


    
    %% Normcorr
    
    info.crop.x = x; 
    info.crop.y = y;
    info.d1 = length(info.crop.x);
    info.d2 = length(info.crop.y);
    % info.skipFrame = 65537; % Skip frame 65537 (for old recordings which
    % have a bug that repeats a frame after 65537 (2^16) frames)
    info.Skipframe = -1; % Don't skip any frames
    
    
    % Show the data before committing the analysis
    Img = sbxread(pnfn, 0,100*nSlices);
    %determine pixel value range
    cLimits = prctile(Img(:), [0.05 94]);
    % Activate much tighter subplots
    % [subplot margin top&side],[figure bottomspace,topspace],[leftspace,rightspace]
    subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.01], [0 0.04], [0.1 0]);
    figure('units','normalized','position',[0.1 0.05 0.25 0.85]);
    handles = gobjects(1,nSlices);
    for j = 1:nSlices
        handles(j) = subplot(nSlices, 1, j);
        imagesc(mean(squeeze(Img(1,:,:,j:nSlices:end)),3))
        caxis(cLimits)
        hold on
        rectangle('position',[x(1), y(1), x(end)-x(1), y(end)-y(1)],'edgecolor',[1 1 0.75],'linewidth',2.5)
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

    % START NORMCORR REGISTRATION % % % % %
    normcorrTic = tic;
    simonalign3;
    normcorrtoc = toc(normcorrTic);
    
    
    filenameNormcorr = cell(nSlices,1);
    for j = 1:nSlices
        filenameNormcorr{j} = [pnfn sprintf('_Split%d_normcorr', j)];
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

    
    %% BACKGROUND IMAGES % % % % % 
    backgrndTic = tic;
    for j = 1:nSlices
        BackgroundImgSbx(filenameNormcorr{j});
    end
    backgrndtoc = toc(backgrndTic);
    
    
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
        getRoistoc = nan;
    end
    
    %% fill in the timedData tracker investigation file
    if timed
        load(timedFile)
        fileInfo = dir([pn fn]);
        timedData.fileSize  (end+1) = fileInfo.bytes * 10^-9;
        timedData.normcorrT (end+1) = normcorrtoc;
        timedData.transposeT(end+1) = transposetoc;
        timedData.decimatT  (end+1) = dectoc;
        timedData.spectralT (end+1) = spectraltoc;
        timedData.backgrndT (end+1) = backgrndtoc;
        timedData.getRoisT  (end+1) = getRoistoc;
        timedData.nRois     (end+1) = nRois;
        timedData.nSplits   (end+1) = nSlices;
        timedData.computerName{end+1} = getenv('COMPUTERNAME');
        timedData.fileName  {end+1} = fn;
        timedData.n = timedData.n + 1;
        save(timedFile, 'timedData')
        fprintf('saved timeddata for file %d\n', i)
    end
end


