% Convert hdf5 files to sbx.
% Step 1: Select h5 file.
% Step 2: Select filename and location to save sbx file.
% Step 3: Type the data sampling rate in Hz
% Step 4: Repeat step 1-3 for as many files as you want to convert. Press
%         cancel when done selecting files.
% 
% Script is only tested on grayscale data.
% 
% Leander de Kraker
% 2021-4-6
% 

% Load as many files as user wants
% Get the tiff file and where to put the sbx file
fileNameH5 = {};
filePathH5 = {};
fileNameSbx = {};
filePathSbx = {};
defInput = {};
selecting = true; % true as long as files are being selected
i = 0;
while selecting
    i = i + 1;
    [fileNameH5{i}, filePathH5{i}] = uigetfile('*.h5',...
            sprintf('get h5 file %d, press cancel when done selecting',i));


    if fileNameH5{i} == 0 % Cancel is pressed probably: stop with selecting
        fileNameH5(i) = [];
        filePathH5(i) = [];
        selecting = false;
    else % Ask where to save the sbx file and other some info
        [fileNameSbx{i}, filePathSbx{i}] = uiputfile([filePathH5{i} fileNameH5{i}(1:end-3) '.sbx'],...
                'Where to save the sbx file');
        [hz{i}, scaleUm{i}, FOVum{i}, pixelAspectRatio{i}, squareFOV{i}, defInput] = RequestRecInfo(defInput);
    end
end
nfiles = length(fileNameH5);
clearvars selecting

i = 1;
%% process the files
for i = 1:nfiles
    
    % global info variable may cause problems, thoroughly delete it
    clearvars -global info
    clearvars info

    fileH5 = [filePathH5{i}, fileNameH5{i}];
    fileSbx = [filePathSbx{i}, fileNameSbx{i}];

    h5Info = h5info(fileH5);
    
    dims = h5Info.Datasets.Dataspace.Size;
    dataName = ['/', h5Info.Datasets.Name];
    ndims = length(dims);
    nframes = dims(3);    
    
    if isfield(h5Info.Datasets.Datatype, 'Type') && strcmp(h5Info.Datasets.Datatype.Type, 'H5T_STD_U8LE')
        % Increase values by multiplying with 256 to allow increase in data
        % detail in later analysis
        multiplier = 2^8;
        fprintf("Increasing data's max values from 256 to 65536 to allow for more detail in 16 bit file\n")
    else
        multiplier = 1;
    end
    
    
    fileID = fopen(fileSbx, 'w'); % open sbx file
    
    
    % Read all frames of the tif and write them to the sbx file
    hwb = waitbar(0, 'writing sbx file');
    countEl =  [dims(1:end-1),1]; % How many elements to read along each dimension

    for f = 1:nframes
        startLoc = [ones(1,ndims-1),f]; % from which sample to start reading the data
        im = uint16(h5read(fileH5, dataName, startLoc, countEl)) .* multiplier;
        fwrite(fileID, im, 'uint16');
        if mod(f, 500)==0
            waitbar(f/nframes, hwb, 'writing sbx file')
            imagesc(im); colormap(cmapL('inferno', 256)); title(sprintf('%d',f))
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
