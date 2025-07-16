% Convert MPEG-4 (MP4) files to sbx.
% The three color channels will be averaged into one grayscale video
% 
% 
% Step 1: Select mp4 file.
% Step 2: Select filename and location to save sbx file.
% Step 3: Type the data sampling rate in Hz
% Step 4: Repeat step 1-3 for as many files as you want to convert. Press
%         cancel when done selecting files.
% 
% In my experience mp4 video suffers from compression artifacts which do
% degrade data quality. The 8 bit integers also do not offer much quality.
% 
% Leander de Kraker
% 2022-4-19
% 

% Load as many files as user wants
% Get the tiff file and where to put the sbx file
fileNameMP4 = {};
filePathMP4 = {};
fileNameSbx = {};
filePathSbx = {};
defInput = cell(5, 1); % default input for information request prompt
defInput{4} = 1;
defInput{5} = 'yes';

selecting = true; % true as long as files are being selected
i = 0;
while selecting
    i = i + 1;
    [fileNameMP4{i}, filePathMP4{i}] = uigetfile('*.mp4',...
        sprintf('get mp4 file %d, press cancel when done selecting',i));


    if fileNameMP4{i} == 0 % Cancel is pressed probably: stop with selecting
        fileNameMP4(i) = [];
        filePathMP4(i) = [];
        selecting = false;
    else % Ask where to save the sbx file
        [fileNameSbx{i}, filePathSbx{i}] = uiputfile([filePathMP4{i} fileNameMP4{i}(1:end-4) '.sbx'],...
            'Where to save the sbx file');
        reader = VideoReader([filePathMP4{i}, fileNameMP4{i}]);
        defInput{1} = reader.FrameRate; % hz info present in file
        [hz{i}, scaleUm{i}, FOVum{i}, pixelAspectRatio{i}, squareFOV{i}, defInput] = RequestRecInfo(defInput);

    end
end
nfiles = length(fileNameSbx);
clearvars selecting

i = 1;
%% process the files
for i = 1:nfiles
    
    % global info variable may cause problems, thoroughly delete it
    clearvars -global info
    clearvars info
    
    fileMP4 = [filePathMP4{i}, fileNameMP4{i}];
    fileSbx = [filePathSbx{i}, fileNameSbx{i}];
    
    reader = VideoReader(fileMP4);
    
    
    dims = [reader.Width, reader.Height, reader.NumFrames];
    nframes = dims(3);
    
    fileID = fopen(fileSbx, 'w'); % open sbx file
    
    
    % Read all frames of the tif and write them to the sbx file
    hwb = waitbar(0, 'writing sbx file');
    countEl =  [dims(1:end-1),1]; % How many elements to read along each dimension
    
    for f = 1:nframes
        im = uint16(mean(read(reader, f),3) * 4)';
        fwrite(fileID, im, 'uint16');
        if mod(f, 500)==0
            waitbar(f/nframes, hwb, 'writing sbx file')
            imagesc(im); colormap(cmapL('inferno', 256)); colorbar; title(sprintf('%d',f))
        end
    end
    fclose(fileID);
    close(hwb)
    
    % Write the info (.mat) file that belongs to the sbx
    info = struct();
    info.sz = dims(1:2); %image dimension
    info.channels = 2; %one channel
    info.Slices = 1; % One depth
    info.nchan = 1;
    info.Perm = [2 1 3 4];
    info.Shape = [dims(1:2), 1];
    info.nsamples = info.sz(2) * info.sz(1) * 2 * info.nchan;
    info.simon = 1; %don't invert with intmax
    info.scanbox_version = 2.5;
    info.max_idx = nframes;
    info.recordsPerBuffer = 0;
    info = RequestRecInfoProcess(info, hz{i}, scaleUm{i}, FOVum{i}, pixelAspectRatio{i}, squareFOV{i});
    info.strfp = [fileSbx(1:end-4)];
    
    save([fileSbx(1:end-4) '.mat'], 'info');
    fprintf('done saving %s\n', fileNameSbx{i})
end
