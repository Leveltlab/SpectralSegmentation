function POW = getcorpow(CSpect, ASpect) 

[Lng, Width, Height] = size(ASpect);

%CCOR = real(ifft(CSpect));
ACOR = real(ifft(ASpect));  %inverse fft, correlation function, 0 delay is variance
VAR = squeeze(ACOR(1,:,:)); %estimate of total variance for each pixel

TVAR = VAR(2:Width-1,2:Height-1);

%surrounding 8 pixels
winx = 1:Width-2;
winy = 1:Height-2;
SVAR = VAR(winx,winy);
VAR2 = shiftdim(repmat(TVAR.*SVAR, 1, 1, Lng), 2);
POW = CSpect(:,:,:,1).*conj(CSpect(:,:,:,1))./VAR2;


winx = 2:Width-1;
winy = 1:Height-2;
SVAR = VAR(winx,winy);
VAR2 = shiftdim(repmat(TVAR.*SVAR, 1, 1, Lng), 2);
POW = POW + CSpect(:,:,:,2).*conj(CSpect(:,:,:,2))./VAR2;

winx = 3:Width;
winy = 1:Height-2;
SVAR = VAR(winx,winy);
VAR2 = shiftdim(repmat(TVAR.*SVAR, 1, 1, Lng), 2);
POW = POW + CSpect(:,:,:,3).*conj(CSpect(:,:,:,3))./VAR2;

winx = 1:Width-2;
winy = 2:Height-1;
SVAR = VAR(winx,winy);
VAR2 = shiftdim(repmat(TVAR.*SVAR, 1, 1, Lng), 2);
POW = POW + CSpect(:,:,:,4).*conj(CSpect(:,:,:,4))./VAR2;

winx = 3:Width;
winy = 2:Height-1;
SVAR = VAR(winx,winy);
VAR2 = shiftdim(repmat(TVAR.*SVAR, 1, 1, Lng), 2);
POW = POW + CSpect(:,:,:,5).*conj(CSpect(:,:,:,5))./VAR2;

winx = 1:Width-2;
winy = 3:Height;
SVAR = VAR(winx,winy);
VAR2 = shiftdim(repmat(TVAR.*SVAR, 1, 1, Lng), 2);
POW = POW + CSpect(:,:,:,6).*conj(CSpect(:,:,:,6))./VAR2;

winx = 2:Width-1;
winy = 3:Height;
SVAR = VAR(winx,winy);
VAR2 = shiftdim(repmat(TVAR.*SVAR, 1, 1, Lng), 2);
POW = POW + CSpect(:,:,:,7).*conj(CSpect(:,:,:,7))./VAR2;

winx = 3:Width;
winy = 3:Height;
SVAR = VAR(winx,winy);
VAR2 = shiftdim(repmat(TVAR.*SVAR, 1, 1, Lng), 2);
POW = POW + CSpect(:,:,:,8).*conj(CSpect(:,:,:,8))./VAR2;

POW = POW/8.0;

