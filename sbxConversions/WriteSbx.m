function WriteSbx(fileSbx, data, varargin)
% WriteSbx(fileSbx, data)
% WriteSbx(fileSbx, data, Hz)
% Write sbx file from workspace variable
% 
% Input: 
%   fileSbx (string): filename, with path possibly, 
%   data (3D double | uint16): video data to write
%   Hz (double): sampling rate of data
% 
% As script:
% Step 1: Call the variable that has the 2p video data: data
% Step 2: Select filename and location to save sbx file.
% Step 3: Set the correct sampling rate
% Step 4: run the code
% 
% 
% 
% Leander de Kraker
% 2021-9-17. 
% 2021-12-13 Turned into function
% 

if exist('varargin', 'var') && nargin == 3
    fileNameSbx = strsplit(fileSbx, '\');
    fileNameSbx = fileNameSbx{end};
    Hz = varargin{1};
else
    Hz = 30;
    filePathSbx = 'D:\NAOMi data\MyData\';
    fileNameSbx = 'typicalVolume';
    fileSbx = [filePathSbx, fileNameSbx];
end

%% process the files

% global info variable may cause problems, thoroughly delete it
clearvars -global info
clearvars info

dims = size(data);
nframes = dims(3);

fileID = fopen([fileSbx, '.sbx'], 'w'); % open sbx file


% Read all frames of the tif and write them to the sbx file
hwb = waitbar(0, 'writing sbx file');

for f = 1:nframes
    fwrite(fileID, uint16(data(:, :, f))', 'uint16');
    if mod(f, 1000)==0
        waitbar(f/nframes, hwb, 'writing sbx file')
        imagesc(uint16(data(:, : , f))); 
        colormap(cmapL('inferno', 256)); title(sprintf('%d',f))
    end
end
fclose(fileID);
close(hwb)

% Write the info (.mat) file that belongs to the sbx
info = struct();
info.sz = dims([2 1]); %image dimension
info.channels = 2; %one channel
info.Slices = 1; % One depth
info.nchan = 1;
info.Perm = [2 1 3 4];
info.Shape = [dims([2 1]), 1];
info.nsamples = info.sz(1) * info.sz(2) * 2 * info.nchan;
info.simon = 1; %don't invert with intmax
info.scanbox_version = 2.5;
info.max_idx = nframes;
info.Freq = Hz;
info.recordsPerBuffer = 0;
info.strfp = fileSbx;

save([fileSbx '.mat'], 'info');
fprintf('done saving %s\n', fileNameSbx)


