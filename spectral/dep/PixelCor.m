function [Con, A, F, V, roundedness, Rvar] = PixelCor(...
                Dim, F, sbxt, py, px, Iy, Ix, yrange, xrange, ThC)
% Correlates all pixels in an ROI selected as valid in roisfromlocalmax, with 
% the median trace associated with the pixels surrounding a localmax in the spectral image
% 2018
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
%Rvar : %Mean covariance of pixels
%
% added correlation cutoff and median filtering of correlations
% Added median correlation of surround to obtain a better cutoff 
% Chris van der Togt, 2021
% Netherlands Institute for Neuroscience 

global DISPLAY
Con = []; %found contour
A = 0; %area of found contour
roundedness = 0; %Roundedness of contour
Rvar = 0; %variance of pixel correlations

se = strel('disk', 6);
Fd = imdilate(F,se); %dilated roi in mask

Fb = Fd - F; %border pixels
se = strel('disk', 2);
Fb = imerode(Fb, se);

Mt =  zeros(Dim);
Mt(yrange, xrange) = F;

Mb =  zeros(Dim);
Mb(yrange, xrange) = Fb;

Mt = logical(Mt');%transposed
Sigpix = double(sbxt.Data.y(:,Mt(:)));
%PreRoiSz = sum(Mt(:));

Mb = logical(Mb');
Sigborder = double(sbxt.Data.y(:,Mb(:)));

if DISPLAY == 1
        figure(4)
        plot(median(Sigpix,2), 'b'), hold on
        plot(median(Sigborder,2), 'g'), hold off
end

Ft = F'; %transpose indices
linix = find(Ft(:));
Seed = zeros(size(Ft));
Seed(px-1:px+1, py-1:py+1) = 1;
Seedx = find(Seed);
[~, ISix] = intersect(linix, Seedx);
Msig = median(Sigpix(:,ISix), 2); %take median of signal of max and surrounding 8 pixels,

Cor = corr(Msig, Sigpix, 'rows', 'pairwise');
Corb = corr(Msig, Sigborder, 'rows', 'pairwise');
Cmx = max(Cor);
Cmd = median(Corb);


% In voxel
V = zeros(size(Ft));
V(linix) = Cor;

%determine correlation strength at which to set contourc level
thr = ThC * (Cmx - Cmd) + Cmd; %difference between center and border

V = V'; %transpose back
if thr < 0.5
    % smooth with median filter, since with lower cutoff more variable
    % correlations are included leading to noisy contours
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
        figure(3), hold off, imagesc(Vimg), caxis([0 1]), hold on
        plot(Con.x, Con.y, 'w')
        colormap([0 0 0 ; parula(256)]); colorbar;
        drawnow

% For debugging
  %  if PreRoiSz >  1.5 * A
        figure(4)
        Rix = find(F);
        [~, Ix] = intersect(linix, Rix);
        [~, Inx] = intersect(linix, find(~F));
        plot(median(Sigpix(:,Ix), 2), 'b')
        hold on
        plot(median(Sigpix(:,Inx), 2), 'g')
 
        Ft = Fd';
        lx = find(Ft(:));
        M =  zeros(Dim);
        M(yrange, xrange) = Fd;
        Mt = logical(M');%transposed
        Sig = sbxt.Data.y(:,Mt(:));
        [W, H] = runnmf(Sig, lx, size(V));
        figure(4), plot(W(:,1)/7 + 4500, 'm'), hold off 
%     end
     pause
end
    
