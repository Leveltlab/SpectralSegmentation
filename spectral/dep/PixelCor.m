function [Con, A, F, V, Ro, Rvar] = PixelCor(Dim, F, sbxt, py, px, Iy, Ix, yrange, xrange, freq, ThC)
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

global DISPLAY
    Con = []; %found contour
    A = 0; %area of found contour
    Ro = 0; %Roundedness of contour
    Rvar = 0; %variance of pixel correlations
    Mt =  zeros(Dim);
    Mt(yrange, xrange) = F;

    Mt = logical(Mt');%transposed
    % Tix = find(Mt); %indices of ROI in transposed mask
    Sigpix = sbxt.Data.y(:,Mt(:));
    R = floor(freq);
    
    Sp = zeros(ceil(size(Sigpix,1)/R), size(Sigpix,2));      
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
        m = floor((1-thr)*3)* 2 + 1; %let median filter depend on threshold height
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
                F(:) = 0;
                Inix = inpolygon(Ix,Iy,vx,vy);
                F(Inix) = 1;
                %roundedness
                Ro = perimarea(vx, vy);
                %mean pixel variance in new contour
                Rvar = mean(V(Inix).^2);
           end
        end
    end
    if A > 0 && DISPLAY == 1        
        figure(3), imagesc(V), hold on
        plot(Con.x, Con.y, 'w')
        drawnow
   end
end
    
