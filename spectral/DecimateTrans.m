function DecimateTrans(varargin)
% decimate a trans.dat file, creating a new transposed file which has the
% signal decimated to 1Hz, with the data format changed to double, instead
% of uint16
% 
% DecimateTrans(strfp, freqDec)
% DecimateTrans([], freqDec)
% DecimateTrans
% 
% Can be executed as script or as function
% input (both are optional):
%   1: strfp: (string): folder and filename of the trans.dat file that will be
%           decimated
%   2: freqDec: (double 1x1): frequency in Hz to decimate data to. default
%                                                                1.0 Hz
% 
% Chris v.d. Togt
% 2020-7-25
% 2020-7-29: Edited by Leander to enable call as function
% 2022-8-25: Edited by Leander to aliviate memory usage, as memory usage of
%               memorymapped file grows as more data is retrieved from it.
% 

%% Get the filename

freqDec = 1.0;
strfp = [];
if exist('varargin', 'var') && nargin == 1
    strfp = varargin{1};
elseif exist('varargin', 'var') && nargin == 2
    strfp = varargin{1};
    freqDec = varargin{2};
end

if ~exist(strfp, 'file')
    [filenameTrans, filePath] = uigetfile('*_Trans.dat', 'select transposed data');   
    strfp = [filePath filenameTrans];
end
strfpBase = split(strfp, '_Trans');

gcp

[dim, frequency] = Zgetinfo(strfp);

sbxt = memmapfile(strfp, 'Format', 'uint16', 'Offset', 500);

fileDecimated = fopen([strfpBase{1} '_DecTrans.dat'],'w');

Lng = dim(1); % length of image stack
width = dim(2) * dim(3); % length * height

steps = ceil(Lng/5000)*10;
chunk = floor(width/steps);
rem = width - steps*chunk;

Decf = round(frequency/freqDec); % decimation factor 
freqDec = frequency/Decf;  % downsampled frequency
LngD1 = ceil(Lng/Decf); % number of samples after decimating
nbytes = fprintf(fileDecimated, "%d %d %d %d %s", LngD1, dim(2), dim(3), freqDec, datetime);
fwrite(fileDecimated,zeros(500-nbytes,1));
D1 = zeros(LngD1, chunk);

% Calculate how many times to reload the memorymapped data
switch sbxt.Format
    case 'uint8'
        multiplier = 1;
    case 'uint16'
        multiplier = 2;
    case 'single'
        multiplier = 4;
    case 'double'
        multiplier = 8;
    otherwise % whatever
        multiplier = 8;
end
fileSize = numel(sbxt.Data)*multiplier ./ (1*10^9);
if ispc
    [~, memSize] = memory;
    memSize = memSize.PhysicalMemory.Available ./ (1*10^9);
    memSize = max(2, memSize); % Do not accept less than 2GB memory
else % Hardoded value for non-windows computers. Pablo B
    memSize = 8;
end
timesToReload = ceil(fileSize/(memSize/3));
timeToReload = round(linspace(0, steps+1, timesToReload+1));

% Load, decimate and write data
hw = waitbar(0, 'processed: 0%', 'Name', 'Decimating transposed dataset');
% timer = zeros(steps, 1); memSystem = zeros(steps, 1); memPhysical = zeros(steps, 1);
for j = 1:steps
%     tic
    indices = (1 + (j-1)*chunk):(j*chunk);
    D = double(XYgetZ(indices, sbxt, Lng));
    D = reshape(D, Lng, [] ); 
    parfor i = 1:length(indices)
        D1(:,i) = decimate(D(:,i), Decf);
    end
    fwrite(fileDecimated, D1,'double');
    waitbar(j/steps, hw, ['processed: ' num2str(round(j/steps*1000)*0.1) '%'])
    
    if ismember(j, timeToReload) % Reload sbxt to aliveate memory usage
        clearvars sbxt
        pause(0.02)
        sbxt = memmapfile(strfp, 'Format', 'uint16', 'Offset', 500);
    end
%     timer(j) = toc; 
%     [~, memmi] = memory; 
%     memSystem(j) = memmi.SystemMemory.Available;
%     memPhysical(j)= memmi.PhysicalMemory.Available;
end

% Write remaining data
D1 = zeros(LngD1, rem);
indices = (1 + steps*chunk):(steps*chunk + rem);
D = double(XYgetZ(indices, sbxt, Lng));
D = reshape(D, Lng, [] );
parfor i = 1:length(indices)
	D1(:,i) = decimate(D(:,i), Decf);
end
fwrite(fileDecimated, D1,'double');
close(hw)
fclose(fileDecimated);
sbxt = [];
fprintf('done with decimating data\n')




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


