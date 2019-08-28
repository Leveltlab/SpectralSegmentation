function [Con, A, F, V] = PixelCor(Dim, F, sbxt, py, px, Iy, Ix, yrange, xrange, freq)
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
%Output
%Con : new contour
%A   : new area
%F   : new voxel mask
%V   : spatial map of pixel correlations in Voxel
%
%Chris van der Togt, 2018

%first extract the traces from the Trans file associated with
%this ROI
    Con = []; %found contour
    A = []; %area of found contour
    Mt =  zeros(Dim);
    Mt(yrange, xrange) = F;

    Mt = logical(Mt');%transposed
    % Tix = find(Mt); %indices of ROI in transposed mask
    Sigpix = sbxt.Data.y(:,Mt(:));
    Sp = [];
    R = floor(freq);      
    for q = 1:size(Sigpix,2)
        Sp(:,q) = decimate(double(Sigpix(:,q)), R); %signals for all pixels in ROI decimated by freq
    end
%some index accounting to find the traces associated with the ROI
%max
    Ft = F'; %transpose indices
    linix = find(Ft(:));
    Seed = zeros(size(Ft));
    Seed(px-1:px+1, py-1:py+1) = 1;
    Seedx = find(Seed);
    [~, Inix] = intersect(linix, Seedx);
    Msig = median(Sp(:,Inix), 2); %take median of signal of max and surrounding 8 pixels,

    Cor = corr(Msig, Sp, 'rows', 'pairwise');
    
%   figure(4), plot(Cor)

% In voxel
    V = zeros(size(Ft));
    V(linix) = Cor;
    V = V';

    mxm = max(abs(Cor));
    mnm = min(abs(Cor));
    thr = min(0.5, (mxm - mnm)/2 + mnm);
    
 %get contour
    c = contourc(V,[thr thr]);
    s = getcontourlines(c);
    v = arrayfun(@(x) eq(x.c,1),s); %closed
    iv = find(v);
    for j = 1:length(iv)
        vx = s(iv(j)).x;
        vy = s(iv(j)).y;
        In = find(inpolygon(px,py,vx,vy),1);
        if ~isempty(In)
           Con.x = vx;
           Con.y = vy;
           A = polyarea(vx,vy);
           F(:) = 0;
           F(inpolygon(Ix,Iy,vx,vy)) = 1;
           
 %  figure(3), imagesc(V), hold on
 %  plot(Con.x, Con.y, 'w')
%   drawnow
%   pause(0.1)
        end
     end
  
    
