% Takes multiple tiff file(s) which consists of multiple grayscale images
%    and converts all of them to one 16 bit unsigned integer sbx format, 
%    with accompanying .mat info file
% 
% This is to combine multiple tiff files of one long recording into one sbx
% file.
%
% Step 1. Select a tiff file and press ok. Repeat this.
% Step 2. Press cancel when done with selecting all the tiff files. 
% Step 3. Set filename and folder for the one sbx file to output. Press ok.
%
%
% Script is only tested on grayscale data.
% 
% Leander de Kraker
% 2021-6-9
%

% Load as many files as user wants
% Get the tiff file and where to put the sbx file
fileNameTif = {};
filePathTif = {};
selecting = true; % true as long as files are being selected
i = 0;
while selecting
    i = i + 1;
    [fileNameTif{i}, filePathTif{i}] = uigetfile('*.tif',...
        sprintf('get tif file %d, press cancel when done selecting',i));
    
    if fileNameTif{i} == 0 % Cancel is pressed probably: stop with selecting
        fileNameTif(i) = [];
        filePathTif(i) = [];
        selecting = false;
    end
end

nfiles = length(fileNameTif);
clearvars selecting

if fileNameTif{end} == 0
    warning('Sbx name not selected, select sbx file name. (Do not press cancel when asked where to save sbx file)')
else
    [fileNameSbx, filePathSbx] = uiputfile([filePathTif{1} fileNameTif{1}(1:end-4) '.sbx'],...
        'Where to save the sbx file');
    
    fileSbx = [filePathSbx, fileNameSbx];
    [hz, scaleUm, FOVum, pixelAspectRatio, squareFOV] = RequestRecInfo();
end

i = 1;
%% process the files
tic


fileID = fopen(fileSbx, 'w'); % open sbx file
hwb = waitbar(0, 'writing sbx file');

for i = 1:nfiles
    
    % global info variable may cause problems, thoroughly delete it
    clearvars -global info
    clearvars info
    
    fileTif = [filePathTif{i}, fileNameTif{i}];
    
    tiffInfo = imfinfo(fileTif);
    if isfield(tiffInfo(1), 'BitDepth') && tiffInfo(1).BitDepth == 8
        % Increase values by multiplying with 256 to allow increase in data
        % detail in later analysis
        multiplier = 2^8;
        fprintf("Increasing data's max values from 256 to 65536 to allow for more detail in 16 bit file\n")
    else
        multiplier = 1;
    end
    
    nframes = length(tiffInfo);
    for f = 1:nframes
        % Read all frames of the tif and write them to the sbx file
        im = uint16(imread(fileTif,f)) .* multiplier;
        fwrite(fileID, im, 'uint16');
    end
    waitbar(i/nfiles, hwb, 'writing sbx file')
end

fclose(fileID);
close(hwb)
% Write the info (.mat) file that belongs to the sbx
info = struct();
info.sz = size(im); %image dimension
info.channels = 2; %one channel
info.Slices = 1; % One depth
info.nchan = 1;
info.Perm = [2 1 3 4];
info.Shape = [size(im), 1];
info.nsamples = info.sz(2) * info.sz(1) * 2 * info.nchan;
info.simon = 1; %don't invert with intmax
info.scanbox_version = 2.5;
info.max_idx = nframes;
info.recordsPerBuffer = 0;
info.strfp = fileSbx(1:end-4);
info = RequestRecInfoProcess(info, hz, scaleUm, FOVum, pixelAspectRatio, squareFOV);

save([fileSbx(1:end-4) '.mat'], 'info');
fprintf('done saving %s\n', fileNameSbx)

toc
