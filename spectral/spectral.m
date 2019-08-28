%  
% calculate crossspectral density for each pixel with it's surrounding 8 pixels
%javaaddpath( [gitpath '\2Pimage\java\Jimutil\Jimutil.jar'])

%license checkout Distrib_Computing_Toolbox
%distcomp.feature( 'LocalUseMpiexec', false )
%start parallel pool
gcp

%% READ SLICES _________________________________
%Hor oriented slice (b=0; y value from Image) or 
%Vert oriented slice (b=1; x value from image)
[fn, fp] = uigetfile('*_Trans.dat');
% ZT = Spect.Zslice([fp fn] );
% dim = ZT.getDimensions();
% freq = ZT.getFrequency();
[dim, frequency] = Zgetinfo([fp fn]);

%split off '_Trans.dat'
mfn = strsplit(fn, '_Trans');
%load([fp mfn{1}(1:end-6) '.mat']);
 
%% SPECTRAL POWER OF SURROUNDING 8 PIXELS

%now open as memory mapped file
sbxt = memmapfile([fp fn], 'Format', 'uint16', 'Offset', 500);

b = 0; Wdth = dim(2);%horizontal orientation
%b = 1; Wdth = dim(3);%vertical orientation
Lines = dim(3)-2; %number of lines to process = height - 2

%find optimal number of lines to process simultaneusly(more lines - more memory) 
f = factor(Lines); %factorize line dimension
c = combnk(1:length(f),2); %all combinations of 2 factors
p = [prod(f(c')) f]; %the products of these combinations + individual factors
[~, i] = min(abs(p-20)); %The closest to 10
W = p(i);
BW = W+2;

%total number of lines should be divisable by W
if rem(Lines, W) > 0
    disp(['Error ' num2str(dim(3)) 'not divisible by ' num2str(W)])
end

%since we decimate the signal to approximately 2Hx, 
Df = floor(frequency/2); %decimation factor 
freq = frequency/Df;     %downsampled frequency

f = 2.^(4:12); %possible segment lengths
Segm = f(f>60*freq); % power of 2 %60 seconds = lowest wavelength
Segm = Segm(1);

Lng = dim(1); %length of image stack
LngDec = floor(Lng/Df);
Lngf = LngDec*Df; %length floored by decimated sampling rate
Dec = zeros(LngDec, Wdth*BW);

sampd2 = Segm/2;
sampd4 = Segm/4;
Sax = (0:sampd4)/(sampd4)/2*freq;


cnt0 = length(1:Segm:LngDec-Segm+1);
cnt1 = length(sampd2+1:Segm:LngDec-Segm+1);
cnt = cnt0 + cnt1;
SPic = zeros(sampd4+1, Wdth-2, W, Lines/W);
Win = hamming(Segm); 

%%
 for r = 1:Lines/W
     %slice = (1:W)+(W-2)*(double(r)-1); %W vertical locations
     %D = ZT.Read(slice,b); %read horizontal lines of z ordered data
     indices = (1:Wdth*BW) + W * Wdth *(r-1);
     D = XYgetZ(indices, sbxt, Lng);
     D = reshape(D, Lng, [] ); 
     D = D(1:Lngf,:);
     parfor i = 1:length(indices)
         Dec(:,i) = decimate(double(D(:,i)), Df);
     end
     
     Dec = reshape(Dec, LngDec, Wdth, BW); 
    %reshape to block of Width x BW with Lng length
     %to obtain 50% overlapping data segments and slicable for parallel computing   
     D0 = Dec(1:cnt0*Segm,:,:);
     D0 = reshape(D0, Segm, cnt0, Wdth, BW);
     D1 = Dec(sampd2+1:sampd2+cnt1*Segm,:,:);  
     D1 = reshape(D1, Segm, cnt1, Wdth, BW); 
     
     D = cat(2, D0, D1);
     %D = reshape(D, Segm, cnt, Wdth, BW);
     
     CSpect = zeros(sampd4+1, Wdth-2, W, 8);
     ASpect = zeros(sampd4+1, Wdth, BW);
         
     parfor j = 1:cnt
         Tmp = double(squeeze(D(:,j,:,:)));
         Tmp = reshape(Tmp, Segm, Wdth * BW);
         Tmp = detrend(Tmp); 
         Tmp = reshape(Tmp, Segm, Wdth, BW);
         [C, A] = getcsdf(Tmp, Segm, Wdth, BW, Win);
         CSpect = CSpect + C;
         ASpect = ASpect + A;
     end
     
     CSpect = CSpect./cnt;
     ASpect = ASpect./cnt;
     
    % Powsum = getcoher(CSpect, ASpect);
     Powsum = getcorpow(CSpect, ASpect);
     SPic(:, :, :, r) = Powsum;
     
     disp([num2str(round(r*1000*W/Lines)/10) '%'])
 end
 
 
  
 SPic = reshape(SPic, sampd4+1, Wdth-2, Lines); 
 SPic = shiftdim(SPic,1);
 SPic = padarray(SPic, [1 1 0]); %make original size by padding borders

 save([fp mfn{1} '_SPSIG.mat'], 'SPic', 'Sax', 'freq')
 
 %delete(gcp)
 %clear all
 
 