% Converts one folder full of tif files which each consist of one frame
% and converts those to one sbx file
% 
% Step 1. Select the folder that contains all the tiff frames. press ok.
% Step 2. Select the folder and filename for the output sbx file
% 
%
% Script is only tested on grayscale data.
% 
% Chris v.d. Togt
% 2022-8-15. Added info and better prompt names - Leander
% 

fp = uigetdir('', 'Select folder that contains the tiff frames');
files = dir([fp '/*.tif*']);
[~, fn] = fileparts(files(1).name);  
[Fn, strSavepath] = uiputfile([fn '.sbx'], 'Select folder and filename of output sbx file');
fileID = fopen([strSavepath Fn], 'w');
[hz, scaleUm, FOVum, pixelAspectRatio, squareFOV] = RequestRecInfo();

hwb = waitbar(0, 'writing sbx file');
for i = 1:length(files)
    fname = fullfile(fp, files(i).name);
    im = imread(fname);
    fwrite(fileID,im,'uint16');
    waitbar(i/length(files), hwb, 'writing sbx file')
end       
fclose(fileID);


info.sz = size(im); %image dimension
info.channels =  2; %one channel
info.nchan = 1;
info.Perm = [2 1 3 4];
info.Shape = [size(im), 1];
info.nsamples = info.sz(2) * info.sz(1) * 2 * info.nchan;
info.simon = 1; %don't invert with intmax
info.scanbox_version = 2.5;
info.max_idx = length(files);
info.recordsPerBuffer = 0;
info = RequestRecInfoProcess(info, hz, scaleUm, FOVum, pixelAspectRatio, squareFOV);

save([strSavepath Fn(1:end-4) '.mat'], 'info');
close(hwb)  
