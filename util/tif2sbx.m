%function tif2sbx()
%take a tifstack and converts to sbx

fp = uigetdir();
files = dir([fp '/*.tif*']);
[~, fn] = fileparts(files(1).name);  
[Fn, strSavepath] = uiputfile([fn '.sbx']);
fileID = fopen([strSavepath Fn], 'w');
        
for i = 1:length(files)
    fname = fullfile(fp, files(i).name);
    im = imread(fname);
    fwrite(fileID,im,'uint16');
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
Strfreq = inputdlg('What is the framerate for this image sequence?');
info.Freq = str2double(Strfreq);
info.recordsPerBuffer = 0;

save([strSavepath Fn(1:end-4)], 'info');
   
