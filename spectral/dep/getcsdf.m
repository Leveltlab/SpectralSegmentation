function [SDF, AUTO] = getcsdf(data, Segment, Width, Height, Win)
% on parallel pool
%
%
%


Seg2 = round(Segment/4)+1; %reduce output
SDF = zeros(Seg2, Width-2, Height-2,  8);

data = bsxfun(@times, data, Win);
Fdata = fft(data);  %Fdata = fft(data, [], 2);
%Fdata(1,:,:) = 0;

%Fdata = Fdata.'; %transposed
Fdatat = conj(Fdata); %complex conjugated transposed
wins = 1:Seg2;

TARG = Fdata(wins,2:Width-1,2:Height-1);
AUTO = Fdata(wins,:,:).*Fdatat(wins,:,:);

%surrounding 8 pixels
winx = 1:Width-2;
winy = 1:Height-2;
SRC = Fdatat(wins,winx,winy);
SDF(:,:,:,1) = TARG.*SRC;

winx = 2:Width-1;
winy = 1:Height-2;
SRC = Fdatat(wins,winx,winy);
SDF(:,:,:,2) = TARG.*SRC;

winx = 3:Width;
winy = 1:Height-2;
SRC = Fdatat(wins,winx,winy);
SDF(:,:,:,3) = TARG.*SRC;

winx = 1:Width-2;
winy = 2:Height-1;
SRC = Fdatat(wins,winx,winy);
SDF(:,:,:,4) = TARG.*SRC;

winx = 3:Width;
winy = 2:Height-1;
SRC = Fdatat(wins,winx,winy);
SDF(:,:,:,5) = TARG.*SRC;

winx = 1:Width-2;
winy = 3:Height;
SRC = Fdatat(wins,winx,winy);
SDF(:,:,:,6) = TARG.*SRC;

winx = 2:Width-1;
winy = 3:Height;
SRC = Fdatat(wins,winx,winy);
SDF(:,:,:,7) = TARG.*SRC;

winx = 3:Width;
winy = 3:Height;
SRC = Fdatat(wins,winx,winy);
SDF(:,:,:,8) = TARG.*SRC;

