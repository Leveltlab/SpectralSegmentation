% Takes tiff file(s) which consists of multiple grayscale images
%    and converts it to 16 bit unsigned integer sbx format, with
%    accompanying .mat info file
% Pop up enables loading of as many tif files as you want. Press cancel
% when done with selecting
%
% Leander de Kraker
% 2020-11-5
%

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
    [fileNameTif{i}, filePathTif{i}] = uigetfile('*.tif',...
            sprintf('get tif file %d, press cancel when done selecting',i));
    
    if fileNameTif{i} == 0 % Cancel is pressed probably: stop with selecting
        fileNameTif(i) = [];
        filePathTif(i) = [];
        selecting = false;
    else % Ask where to save the sbx file
        [fileNameSbx{i}, filePathSbx{i}] = uiputfile([filePathTif{i} fileNameTif{i}(1:end-4) '.sbx'],...
                'Where to save the sbx file');
         Hz(i) = str2double(inputdlg('What is the framerate for this image sequence? (Hz)'));
    end
end
nfiles = length(fileNameTif);
clearvars selecting

i = 1;
%% process the files
for i = 1:nfiles
    
    % global info variable may cause problems, thoroughly delete it
    clearvars -global info
    clearvars info

    fileTif = [filePathTif{i}, fileNameTif{i}];
    fileSbx = [filePathSbx{i}, fileNameSbx{i}];

    fileID = fopen(fileSbx, 'w'); % open sbx file

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
    
    % Read all frames of the tif and write them to the sbx file
    hwb = waitbar(0, 'writing sbx file');
    for f = 1:nframes
        im = uint16(imread(fileTif, f)) .* multiplier;
        fwrite(fileID, im, 'uint16');
        waitbar(f/nframes, hwb, 'writing sbx file')
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
    info.Freq = Hz(i);
    info.recordsPerBuffer = 0;
    info.strfp = [fileSbx(1:end-4)];

    save([fileSbx(1:end-4) '.mat'], 'info');
    fprintf('done saving %s\n', fileNameSbx{i})
end

