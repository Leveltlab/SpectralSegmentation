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

%% Get the filename
if exist('varargin', 'var') && nargin == 1
    strfp = varargin{1};
    [fp,fn] = fileparts(strfp);
    fp = [fp '\'];
else
    [fn, fp] = uigetfile('*_DecTrans.dat');  
    strfp = [fp fn];
end

[dim, freq] = Zgetinfo(strfp);

mfn = strsplit(fn, '_DecTrans');

gcp


%% SPECTRAL POWER OF SURROUNDING 8 PIXELS

%now open as memory mapped file
sbxt = memmapfile(strfp, 'Format', 'double', 'Offset', 500);

Wdth = dim(2);%horizontal orientation
Lines = dim(3)-2; %number of lines to process = height - 2

%find optimal number of lines to process simultaneusly(more lines - more memory) 
f = factor(Lines); %factorize line dimension
c = combnk(1:length(f),2); %all combinations of 2 factors
p = [prod(f(c')) f]; %the products of these combinations + individual factors
[~, i] = min(abs(p-20)); %The closest to 10
W = p(i);
BW = W+2;
steps = Lines/W;

%total number of lines should be divisable by W
if rem(Lines, W) > 0
    disp(['Error ' num2str(dim(3)) 'not divisible by ' num2str(W)])
end

Segm = 128;
Lng = dim(1); %length of image stack


sampd2 = Segm/2;
sampd4 = Segm/4;
Sax = (0:sampd4)/(sampd4)/2*freq;


cnt0 = length(1:Segm:Lng-Segm+1);
cnt1 = length(sampd2+1:Segm:Lng-Segm+1);
cnt = cnt0 + cnt1;
SPic = zeros(sampd4+1, Wdth-2, W, Lines/W);
Win = hamming(Segm); 

tic
%%
hw = waitbar(0, 'Estimating Cross-spectral Power');
for r = 1:steps
    %slice = (1:W)+(W-2)*(double(r)-1); %W vertical locations
    %D = ZT.Read(slice,b); %read horizontal lines of z ordered data
    indices = (1:Wdth*BW) + W * Wdth *(r-1);
    D = XYgetZ(indices, sbxt, Lng);
    Dec = reshape(D, Lng, Wdth, BW); 

    %reshape to block of Width x BW with Lng length
    %to obtain 50% overlapping data segments and slicable for parallel computing   
    D0 = Dec(1:cnt0*Segm,:,:);
    D0 = reshape(D0, Segm, cnt0, Wdth, BW);
    D1 = Dec(sampd2+1:sampd2+cnt1*Segm,:,:);  
    D1 = reshape(D1, Segm, cnt1, Wdth, BW); 

    Dc = cat(2, D0, D1);

    CSpect = zeros(sampd4+1, Wdth-2, W, 8);
    ASpect = zeros(sampd4+1, Wdth, BW);

    parfor j = 1:cnt
        Tmp = squeeze(Dc(:,j,:,:));
        Tmp = reshape(Tmp, Segm, Wdth * BW);
        Tmp = detrend(Tmp); 
        Tmp = reshape(Tmp, Segm, Wdth, BW);
        [C, A] = getcsdf(Tmp, Segm, Wdth, BW, Win);
        CSpect = CSpect + C;
        ASpect = ASpect + A;
    end

    CSpect = CSpect./cnt;
    ASpect = ASpect./cnt;

    Powsum = getcorpow(CSpect, ASpect);
    SPic(:, :, :, r) = Powsum;

    waitbar(r/steps, hw, ['Estimating Cross-spectral Power: ' num2str( round(r/steps * 1000)/10) '%' ], '%3.1f') 
end
close(hw)

SPic = reshape(SPic, sampd4+1, Wdth-2, Lines); 
SPic = shiftdim(SPic,1);
SPic = padarray(SPic, [1 1 0]); %make original size by padding borders

save([fp mfn{1} '_SPSIG.mat'], 'SPic', 'Sax', 'freq')
toc

