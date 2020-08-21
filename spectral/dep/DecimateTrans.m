function DecimateTrans(varargin)
% decimate a trans.dat file, creating a new transposed file which has the
% signal decimated to 1Hz, with the data format changed to double, instead
% of uint16
% 
% Can be executed as script or as function
% input:
%   string, folder and filename of the trans.dat file that will be
%           decimated
% 
% Chris v.d. Togt
% 2020-7-25
% Edited by Leander to enable call as function
% 2020-7-29
% 

%% Get the filename
if exist('varargin', 'var') && nargin == 1
    strfp = varargin{1};
    [filePath,filenameTrans] = fileparts(strfp);
    filePath = [filePath '\'];
else
    [filenameTrans, filePath] = uigetfile('*_Trans.dat');   
    strfp = [filePath filenameTrans];
end

filenameBase = split(filenameTrans, '_Trans');

gcp

[dim, frequency] = Zgetinfo(strfp);

sbxt = memmapfile(strfp, 'Format', 'uint16', 'Offset', 500);

fileDecimated = fopen([filePath filenameBase{1} '_DecTrans.dat'],'w');

Lng = dim(1); % length of image stack
width = dim(2) * dim(3); % length * height

steps = ceil(Lng/5000)*10;
chunk = floor(width/steps);
rem = width - steps*chunk;

Df1 = ceil(frequency); % decimation factor 
freq = frequency/Df1;  % downsampled frequency (about 1 Hz)
LngD1 = ceil(Lng/Df1); % number of samples after decimating
nbytes = fprintf(fileDecimated, "%d %d %d %d %s", LngD1, dim(2), dim(3), freq, datetime);
fwrite(fileDecimated,zeros(500-nbytes,1));
D1 = zeros(LngD1, chunk);

hw = waitbar(0, 'processed: 0%', 'Name', 'Decimating transposed dataset');
for j = 1:steps
    indices = (1 + (j-1)*chunk):(j*chunk);
    D = double(XYgetZ(indices, sbxt, Lng));
    D = reshape(D, Lng, [] ); 
    parfor i = 1:length(indices)
        D1(:,i) = decimate(D(:,i), Df1);
    end
    fwrite(fileDecimated, D1,'double');
    waitbar(j/steps, hw, ['processed: ' num2str(round(j/steps*1000)*0.1) '%'])
end

D1 = zeros(LngD1, rem);
indices = (1 + steps*chunk):(steps*chunk + rem);
D = double(XYgetZ(indices, sbxt, Lng));
D = reshape(D, Lng, [] );
parfor i = 1:length(indices)
	D1(:,i) = decimate(D(:,i), Df1);
end
fwrite(fileDecimated, D1,'double');
close(hw)
fclose(fileDecimated);

%% Reading the Dec file  %%
% strfp = [fp mfn{1} '_D1.dat']
% fileID = fopen(strfp,'r');
% tline1 = fgetl(fileID);
% D = split(tline1, ' ');
% lng = str2double(D{1});
% width =  str2double(D{2});
% height = str2double(D{3});
% freq = str2double(D{4});
% fclose(fileID)
%
% sbxt = memmapfile(strfp, 'Format', 'double', 'Offset', 500);
% D = sbxt.Data(1:lng*width*height);
% D = reshape(D, lng, []);
% Dm = mean(D);
% Dm = reshape(Dm, width, height);
% figure, surf(Dm), shading interp, view(2), axis tight


