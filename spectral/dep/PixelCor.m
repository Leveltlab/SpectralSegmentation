function [Con, A, F, V, roundedness, Rvar] = PixelCor(...
                Dim, F, sbxt, py, px, Iy, Ix, yrange, xrange, freq, ThC)
% Correlates all pixels in an ROI selected as valid in roisfromlocalmax, with 
% the median trace associated with the pixels surrounding a max in the spectral image
%
%takes as input:
%Dim: for dimension purposes (Image size)
%F :  Mask of selected ROI in Voxel
%sbxt : memory mapped transposed image stack
%px : x position of Max in F
%py : y position of Max in F
%Ix : x indices of mask F
%Iy : y indices of mask F
%xrange : x indices of Voxel in Spectral Image
%yrange : y indices of Voxel in Spectral Image
%freq   : sampling frequency
%ThC    : Correlation threshold cutoff for contourc height
%Output
%Con : new contour
%A   : new area
%F   : new voxel mask
%V   : spatial map of pixel correlations in Voxel
%Ro  : Roundedness
%Rvar : %variance of pixel correlations
%
%first extract the traces from the Trans file associated with
%this ROI
%Chris van der Togt, 2018
%added correlation cutoff and median filtering of correlations
%Chris van der Togt, 2019
% Leander de Kraker 2020-7-3 changed decimate to movmean filtering,

global DISPLAY
Con = []; %found contour
A = 0; %area of found contour
roundedness = 0; %Roundedness of contour
Rvar = 0; %variance of pixel correlations
Mt =  zeros(Dim);
Mt(yrange, xrange) = F;

Mt = logical(Mt');%transposed
sig = sbxt.Data.y(:,Mt(:));
R = floor(freq)*2; % mean filter and subsample to 0.5Hz
sig = movmean(sig, R, 1);
sig = sig(1:R:end, :);

%some index accounting to find the traces associated with the ROI
%max
Ft = F'; %transpose indices
linix = find(Ft(:));
Seed = zeros(size(Ft));
Seed(px-1:px+1, py-1:py+1) = 1;
Seedx = find(Seed);
[~, Inix] = intersect(linix, Seedx);
Msig = median(sig(:,Inix), 2); %take median of signal of max and surrounding 8 pixels,

Cor = corr(Msig, sig, 'rows', 'pairwise');

%  figure(4), plot(Cor) 

% In voxel
V = zeros(size(Ft));
V(linix) = Cor;

%determine correlation strength at which to set contourc level
thr = ThC * max(abs(Cor)); 

V = V'; %transpose back
if thr < 0.5
 %smooth with median filter, since with lower cutoff more variable
 %correlations are included leading to noisy contours
    m = 3;
    if thr < 0.25
        m = 5;
    end
    J = medfilt2(V, [m m]);
else
    J = V;
end

%get contour
c = contourc(J,[thr thr]);
s = getcontourlines(c);
v = arrayfun(@(x) eq(x.c,1),s); %which contours are closed
iv = find(v); %select these only

for j = 1:length(iv)
    vx = s(iv(j)).x;
    vy = s(iv(j)).y;
    %is my spectral max (px, py) in this contour?
    In = find(inpolygon(px,py,vx,vy),1);
    if ~isempty(In)
        %does this one have a valid area
       Atmp = polyarea(vx,vy);
       if Atmp > A
            Con.x = vx;
            Con.y = vy;
            A = Atmp;
            F(:) = 0; %surround pixels
            Inix = inpolygon(Ix,Iy,vx,vy);
            F(Inix) = 1; %pixels in contour
            %roundedness
            roundedness = perimarea(vx, vy);
            %mean pixel variance in new contour
            Rvar = mean(V(Inix).^2);
       end
    end
end
if A > 0 && DISPLAY
    Vimg = V;
    Vimg(Vimg==0) = nan;
    figure(3), hold off, imagesc(Vimg), hold on
    plot(Con.x, Con.y, 'w')
    colormap([0 0 0 ; parula(256)]); colorbar;
    drawnow
end
    
