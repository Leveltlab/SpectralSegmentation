function [Savedrois, Mask, SC, rejlog] = roisfromlocalmax(IM, Savedrois, Mask, varargin)
% roisfromlocalmax(IM, threshold, distance = 10, area = 20, border = 10, )
% IM : image
% Savedrois : array of saved Roi parameters
% Mask : rois as image regions with index value
% threshold : 0 - 1  minimum amount of correlation 
% distance: minimal distance between rois
% area: minimal area of an roi
% border: edge of image to ignore
% 
% this function first determines the local maxima in an image, then
% for a certain fraction of higest points finds a contour around each maximum, 
% if not already contained in a contour.
% then discards overlapping contours, or separates contours in different
% rois.
% Chris van der Togt, 2017, Netherlands Institute for Neuroscience 

global DISPLAY
% Parameters for ROI search
distance = 10;     % Minimum distance between roi maxima
fractionmax = 0.2; % Higest fraction of maxima (pnt) to use
threshold = 0.80;  % Contour threshold greater than percentile of pixel range
Significance = 0.90; % Maximum should be greater than percentile of pixel range
% Parameters which can be chosen by spar
area = [60 120];   % Minimum/max size of roi
border = 10;       % Edge of image to ignore
VoxelSz = 50;      % Determines area to find contours
                   % All rois should be within the size of an area of VoxelSz x VoxelSz
PAf = 0.70;        % Minimal ratio perimeter to squared area; roundedness factor (circle = 1.0)
cutOffCorr = 0.5;  % Cutoff threshold for pixel correlations


if ~isempty(varargin)
    % Minimum distance between two correlation maxima
    p = varargin{1};
        
    if isfield(p, 'areasz')
        area = p.areasz;
    end
    if isfield(p, 'border')
        border = p.border;
    end
    if isfield(p, 'voxel')
        VoxelSz = p.voxel;
        if mod(VoxelSz, 2) == 1 % Voxel size should be an even number
            VoxelSz = VoxelSz + 1;
        end
    end
    if isfield(p, 'roundedness')
        PAf = p.roundedness;
    end
    if isfield(p, 'cutOffCorr')
        cutOffCorr = p.cutOffCorr;
    end
end

if length(varargin) > 1
    sbxt = varargin{2}; %memory mapped file to transposed image data
    freq = varargin{3}; %frequency of transposed images
    SC = varargin{4}; % SpatialCorr
end

% Get reasonable maxima sorted on height
%pnt = mylocalmax(IM, border);
pnt = localmax(IM, border, VoxelSz);

% Reduce number of points to 1/4 with largest peak
Nmp = round(length(pnt)*fractionmax);
dim = size(IM);
%hold on, plot(pnt(1:Nmp, 2), pnt(1:Nmp, 1), 'r+')

% Size of image sample in pixels
SzImg = VoxelSz*VoxelSz;
thx = round(threshold*SzImg); %threshold number
SigTh = round(Significance*SzImg);

[Iy, Ix] = find(ones(VoxelSz,VoxelSz));

rejlog = nan(Nmp, 4);
Cnt = Savedrois.Cnt;
for i = 1:Nmp        
    % Detect if this point is on a previously selected roi
    if ~(Mask(pnt(i,1),pnt(i,2)) > 0)
        pVal = pnt(i,3); %pixelvalue of maximum
        lyw = pnt(i,1)-VoxelSz/2;
        hyw = lyw + VoxelSz-1;
        if lyw < 1, lyw = 1; hyw = VoxelSz; end
        if hyw > dim(1), hyw = dim(1); lyw = hyw - VoxelSz+1; end

        lxw = pnt(i,2)-VoxelSz/2;
        hxw = lxw+VoxelSz-1;
        if lxw < 1, lxw = 1; hxw = VoxelSz; end
        if hxw > dim(2), hxw = dim(2); lxw = hxw - VoxelSz+1;end

        M = Mask(lyw:hyw, lxw:hxw); 
        I = IM(lyw:hyw, lxw:hxw);   % Voxel for contours

        % Cut out previous rois setting to zero
        I(M>0) = 0;            % Set previous ROI locations to zero
        If = imgaussfilt(I,1); % Smooth sample

        Ws = sort(If(:));            % Range of pixel values              
        th = Ws(thx);                % Pixel value threshold
        thxd = 0.1*(pVal-th); %5 % Difference between max and threshold as step

        py = pnt(i,1)-lyw+1;
        px = pnt(i,2)-lxw+1;
        RoiMx = If(py, px);

        if RoiMx > Ws(SigTh) % Pixel value maximum should be significantly higher than surrounding pixels
            it = 0;
            bval = false;
            NwA = [];
            Ro = [];
            while ~bval && it < 10
                [Con, A, F, Pin, Ro] = getCon(If, th, area, PAf*0.9, py, px, Iy, Ix);
                if ~isempty(Con) %valid contour: contains point and has valid roundedness
                    if DISPLAY == 1
                        figure(2), hold off, imagesc(I), colormap gray, hold on
                        plot(px, py, '+r')
                        plot(Con.x, Con.y, 'r')
                        title(sprintf('seedpoint x: %d, y: %d', pnt(i,1), pnt(i,2)))
                    end
                    Con.y = Con.y + lyw - 1;
                    Con.x = Con.x + lxw - 1;
                    bval = true; %valid maximum not in previously selected contour
                    if Cnt > 0
                        %if other points of previous rois are in this contour
                        %determine if it is simply overlapping or much larger
                        In = find(inpolygon(Savedrois.P(2,:),Savedrois.P(1,:), Con.x, Con.y),1);
%                         figure(3), imagesc(Mask(lyw:hyw, lxw:hxw) + F.*Cnt)
                        if ~isempty(In)
                            %the overlapping contours should be at least area(2) larger,
                            %and the containing maximum further than
                            %<b>distance<\b> pixels away, otherwize
                            %ignore this contour and pnt
                            if (A - Savedrois.A(In) > area(2)) && ( sqrt( (Savedrois.P(2,In)- pnt(i,2))^2 + (Savedrois.P(1,In)- pnt(i,1))^2) > distance )
                                %remove values of previous roi from
                                %sample using the mask and smooth sample
                                I(M==In) = 0;
                                I = imgaussfilt(I,2);
                                %try again with higher threshold
                                th = th + thxd;
                            else
                                it = 10; %skip this point 
                            end
                            bval = false; 
                        end
                    end

                    if bval
                        yrange = lyw:hyw;
                        xrange = lxw:hxw;                             
                        [NwCon, NwA, NwF, NwV, Ro, Rvar] = PixelCor(size(Mask), F, sbxt, py, px, Iy, Ix, yrange, xrange, cutOffCorr);
                        if ~isempty(NwCon) &&  NwA > area(1) && Ro >= PAf % Exists and area greater than minimum, and roundedness > roundedness factor
                            Con = NwCon;
                            Con.y = Con.y + lyw - 1;
                            Con.x = Con.x + lxw - 1;
                            F = NwF; %new voxel mask (roi pixels:1, surround: 0)
                        
                            SCM = SC(lyw:hyw, lxw:hxw);
                            SCM(F>0) = NwV(F>0); % Substitute
                            SC(lyw:hyw, lxw:hxw) = SCM; % Update correlation map R^2

                            Cnt = Cnt + 1;
                            Savedrois.Con(Cnt) = Con;
                            Savedrois.A(Cnt) = NwA;
                            Savedrois.Rvar(Cnt) = Rvar;
                            Savedrois.Roundedness(Cnt) = Ro;
                            Savedrois.P(1:3,Cnt) = pnt(i,:)';
                            
                             M(F>0) = Cnt;
                            Mask(lyw:hyw, lxw:hxw) = M;
                            
                            rejlog(i,:) = [2 NwA, Ro, Rvar];
                            if DISPLAY == 1
                                figure(3), 
                                title(sprintf('good!: rounded %.2f, rvar %.2f,', Ro, Rvar)) 
                                %pause
                            end
                        elseif NwA > area(1) && Ro < PAf  %not round enough but large area, => increase threshold
                            bval = false;
                            th = th + thxd; 
                            
                        elseif DISPLAY == 1
                                figure(3), 
                                title([ 'rejected: A ' num2str(NwA) ', Ro '  num2str(Ro) ', Rvar ' num2str(Rvar)] )
                                %pause
                        else 
                            rejlog(i,:) = [isempty(NwCon) NwA, Ro, 0];
                        end
                    end

                elseif A <= area(1) && Pin   % The area was too small, no valid contours
                    % but there were some closed contours with the peak inside, 
                        th = th - thxd; % Lower the threshold  
                else
                    th = th + thxd; % Increase the threshold, none closed or too large
                end
                it = it + 1;
            end
            if bval == false || it >= 10
                rejlog(i,:) = [isempty(Con), A, Ro, 0];
                if DISPLAY == 1
                        figure(3), 
                        title([ 'rejected: A ' num2str(NwA) ', Ro '  num2str(Ro) ] )
                        %pause
                end 

            end
        else
               rejlog(i,:) = [-1, 0, 0, 0];
        end
    else % This point was in an existing contour
        rejlog(i,:) = [nan, nan, pnt(i,2), pnt(i,1)];        
    end
    if(isnan(rejlog(i,1)) && isnan(rejlog(i,2)) && isnan(rejlog(i,3)))
        disp(i)
    end
end
Savedrois.Cnt = Cnt;


