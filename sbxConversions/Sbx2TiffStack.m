% Select as many sbx files as you want to convert to tiff stack
% tif stack will consist of many grayscale images. 
% Tif can't support files larger than 4 GB,  sometimes not even larger 
% than 2 GB. So large sbx get split up into  multiple different files
% smaller than 3 GB...
% 
% Sbx files currently supported: 1 split
% Only 1 color will be written to tiff file
%
% 
% Leander de Kraker 
% 2021-2-24
% 

maxSize = 3; % max tiff file size in GB

% Load as many files as user wants
% Get the tiff file and where to put the sbx file
fileNameTif = {};
filePathTif = {};
fileNameSbx = {};
filePathSbx = {};
selecting = true; % true as long as files are being selected
i = 0;
while selecting
    i = i + 1;
    [fileNameSbx{i}, filePathSbx{i}] = uigetfile('*.sbx',...
            sprintf('get sbx file %d, press cancel when done selecting',i));

    if fileNameSbx{i} == 0 % Cancel is pressed probably: stop with selecting
        fileNameSbx(i) = [];
        filePathSbx(i) = [];
        selecting = false;
    else % Ask where to save the sbx file
        fileNameSbx{i} = strsplit(fileNameSbx{i}, '.');
        fileNameSbx{i} = fileNameSbx{i}{1};
        [fileNameTif{i}, filePathTif{i}] = uiputfile([filePathSbx{i} fileNameSbx{i} '.tif'],...
                'Where to save the tif file');
    end
end
nfiles = length(fileNameTif);
clearvars selecting
i = 1;
%% Convert the files to tif

for i = 1:nfiles
    
    nametifBase = [filePathTif{i}, fileNameTif{i}];
    namesbx = [filePathSbx{i}, fileNameSbx{i}];
    
    % Get info of the sbx file
    clearvars -global info
    clearvars info
    load(namesbx)
    frame = sbxread(namesbx, 0, 1);
    
    sbxfileProperties = dir([namesbx '.sbx']);
    sbxfileSize = sbxfileProperties.bytes / 10^9; %   GB
    if sbxfileSize > maxSize % if file size is over max size..
        tn = ceil(sbxfileSize/ maxSize); % number of tif files to create
        framesForTifs = round(linspace(0, info.max_idx, tn+1)); % frames for each tif file
        nametif = cell(tn, 1);
        for tifi = 1:tn
            nametif{tifi} = strrep([nametifBase(1:end-4), sprintf('_%3d', tifi), '.tif'], ' ', '0');
        end
        fprintf('sbx file is so large it requires %d tif files...\n', tn)
        wbar = waitbar(0, sprintf('converting sbx file %d to tif file (0/%d)', i, tn));
    else % sbx file can be stored in 1 tif file
        framesForTifs = [0 info.max_idx];
        tn = 1;
        nametif{1} = nametifBase;
        wbar = waitbar(0, sprintf('converting sbx file %d to tif', i));
    end
    
    % tag struct for tif info
    tagstruct.ImageLength = size(frame,1);
    tagstruct.ImageWidth = size(frame,2);
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    if isa(frame, 'uint18')
        tagstruct.BitsPerSample = 8;
    elseif isa(frame, 'uint16')
        tagstruct.BitsPerSample = 16;
    elseif isa(frame, 'uint32')
        tagstruct.BitsPerSample = 32;
    else
        warning('Not an expected format, quiting\n\n')
        return
    end
    tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.RowsPerStrip = 256;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Software = 'MATLAB';
    
    % Write sbx file to to tn tif files
    tic
    for tifi = 1:tn
        waitbar(framesForTifs(tifi)/info.max_idx, wbar,...
                sprintf('converting sbx file %d to tif file (%d/%d)', i, tifi, tn))
        
        % Get which frames to put into this tif file
        if tn == 1 % All frames
            framesToDo = 0:info.max_idx-1;
        else % frames between file limit i and i+1
            framesToDo = framesForTifs(tifi):framesForTifs(tifi+1)-1;
        end
        
        % open tiff file
        t = Tiff(nametif{tifi},'w');
        % write the frames
        for framej = framesToDo
            frame = sbxread(namesbx, framej, 1);
            t.setTag(tagstruct);
            t.write(squeeze(frame(:,:,1,:)));
            t.writeDirectory();
            if mod(framej, 50)==0
                waitbar(framej/info.max_idx, wbar)
            end
        end
        t.close();
    end
    toc
    delete(wbar)
end

