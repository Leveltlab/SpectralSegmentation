function spectral(varargin)
%  
% calculate cross-spectral density for each pixel with it's surrounding 8 pixels
% When input is not given a pop-up will ask for a filename
% The file should contain binary data with time in the first dimension and decimated to
% approximately 1Hz.
%
% Chris van der togt, 01-07-2020 
% Netherlands Institute for Neuroscience, 
% Amsterdam, the Netherlands
% INPUT: filename_DecTrans.dat; transposed decimated data of an aligned calcium imaging movie.
% OUTPUT: Saves data to SPSIG.mat file
% UPDATE: 01-05-2021 Number of steps estimated from working memory and worker number 

%% Get the filename
if exist('varargin', 'var') && nargin == 1
    strfp = varargin{1};
    [fp,fn] = fileparts(strfp);
    fp = [fp '\'];
else
    [fn, fp] = uigetfile('*_DecTrans.dat');  
    strfp = [fp fn];
end
mfn = strsplit(fn, '_DecTrans');

[dim, freq] = Zgetinfo(strfp);

P = gcp;
NumWorkers = P.NumWorkers;

Segm = 128; %fft window length for 1Hz sampling rate
Lng = dim(1); %length of image stack
Wdth = dim(2);%horizontal orientation
Lines = dim(3)-2; %number of lines to process = height - 2
sampd2 = Segm/2;
sampd4 = Segm/4;

BytesA = 8 * 4 * Wdth * Lng;%(width of line, length of trace * 4 * 8 bytes)
BytperLine = BytesA * (NumWorkers + 1);

%number of lines to process given amount(0.5) of available memory and estimated memory
%usage per line.
usr = memory();
SysMem = usr.MemAvailableAllArrays;
BW = floor(SysMem/2/BytperLine);
W = BW-2;
if W < 1  %set minimum of W and BW
    W = 1;
    BW = 3;
end
steps = floor(Lines/W);
if steps < 1 %set maximum of W and BW
    steps = 1;
    W = Lines;
    BW = Lines+2;
end

lastW = 0;
%remaining step if necessary
if rem(Lines, W) > 0
    steps = steps + 1;
    lastW = rem(Lines, W);
end
nwLines = steps*W;

%% SPECTRAL POWER OF SURROUNDING 8 PIXELS

%now open as memory mapped file
sbxt = memmapfile(strfp, 'Format', 'double', 'Offset', 500);

Sax = (0:sampd4)/(Segm)/2*freq;

cnt0 = length(1:Segm:Lng)-1;
cnt1 = length(sampd2+1:Segm:Lng)-1;
cnt = cnt0 + cnt1;
SPic = zeros(sampd4+1, Wdth-2, W, steps);
Win = hamming(Segm); 

tic
%%
disp(['Runs in ' num2str(steps) ' steps, please wait.'])

hwb = waitbar(0, 'Estimating Cross-spectral Power');
q = parallel.pool.DataQueue;
afterEach(q, @updatewb);

for r = 1:steps
    Skip = W;
    if (r == steps && lastW > 0) 
        BW = lastW+2;
        W = lastW;
    end
    indices = (1:Wdth*BW) + Skip * Wdth *(r-1);
     
    Dec = XYgetZ(indices, sbxt, Lng);
    Dec = reshape(Dec, Lng, Wdth, BW); 

    %reshape to block of Width x BW with Lng length
    %to obtain 50% overlapping data segments and slicable for parallel computing   
    D0 = Dec(1:cnt0*Segm,:,:);
    D0 = reshape(D0, Segm, cnt0, Wdth, BW);
    D1 = Dec(sampd2+1:sampd2+cnt1*Segm,:,:);  
    D1 = reshape(D1, Segm, cnt1, Wdth, BW); 
    Dc = cat(2, D0, D1);
    
    CSpect = zeros(sampd4+1, Wdth-2, W, 8);
    ASpect = zeros(sampd4+1, Wdth, BW);
    p = 0;  
    parfor j = 1:cnt
        Tmp = squeeze(Dc(:,j,:,:));
        Tmp = reshape(Tmp, Segm, Wdth * BW);
        Tmp = detrend(Tmp); 
        Tmp = reshape(Tmp, Segm, Wdth, BW);
        [C, A] = getcsdf(Tmp, Segm, Wdth, BW, Win);
        CSpect = CSpect + C;
        ASpect = ASpect + A;
        send(q, j)
    end

    CSpect = CSpect./cnt;
    ASpect = ASpect./cnt;

    Powsum = getcorpow(CSpect, ASpect);
    SPic(:,:,1:W,r) = Powsum;

end
close(hwb)

SPic = reshape(SPic, sampd4+1, Wdth-2, nwLines);
SPic = SPic(:,:,1:Lines);
SPic = shiftdim(SPic,1);
SPic = padarray(SPic, [1 1 0]); %make original size by padding borders

save([fp mfn{1} '_SPSIG.mat'], 'SPic', 'Sax', 'freq')
toc

%%
    function updatewb(~)
        p = p + 1;
        fraction = p/cnt; 
        waitbar(fraction, hwb, ['Estimating Cross-spectral Power, Step ' num2str(r) ': ' num2str(round(fraction * 1000)/10) '%' ], '%3.1f')
    end
end
