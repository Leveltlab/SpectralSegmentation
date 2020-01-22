function [Savedrois, Mask, SC] = roisfromlocalmax(IM, Savedrois, Mask, varargin)
% roisfromlocalmax(IM, threshold, distance = 10, area = 20, border = 15)
% IM : image
% Mask : rois as image regions with index value
% threshold : 0 - 1  minimum amount of correlation 
% distance: minimal distance between rois
% area: minimal area of an roi
% border: edge of image to ignore
% 
% this function first determines the local maxima in an image, then
% for the 1/4 higest points finds a contour around each maximum, 
% if not already contained in a contour.
% then discards overlapping contours, or separates contours in different
% rois.
% Chris van der Togt, 2017, Netherlands Institute for Neuroscience 
global DISPLAY  
distance = 20;     % Distance between rois (maxima)
area = [60 120];   % Minimum/max size of roi
border = 15;       % Edge of image to ignore
fractionmax = 0.25;% Higest fraction of maxima (pnt) to use
VoxelSz = 60;      % Determines area to find contours
% All rois should be within the size of an area of VoxelSz x VoxelSz    
threshold = 0.90;  % Contour threshold greater than percentile of pixel range
Significance = 0.90; % Maximum shoud be greater than percentile of pixel range
PAf = 0.70;        % Minimal ratio perimeter to squared area; roundedness factor (circle = 1.0)
cutoffcorr = 0.5;  % Cutoff threshold for pixel correlations


if ~isempty(varargin)
    % Minimum distance between two correlation maxima
    p = varargin{1};
    if isfield(p, 'pixelthreshold')
        threshold = p.pixelthreshold;
    end
    if isfield(p, 'significance')
        Significance = p.significance;
    end
    if isfield(p, 'fractionmax')
        fractionmax = p.fractionmax;
    end
    if isfield(p, 'minimaldistance')
        distance = p.minimaldistance;
    end         
    if isfield(p, 'areasz')
        area = p.areasz;
    end
    if isfield(p, 'border')
        border = p.border;
    end
    if isfield(p, 'roundedness')
        PAf = p.roundedness;
    end
    if isfield(p, 'voxel')
        VoxelSz = p.voxel;
    end
    if isfield(p, 'cutoffcorr')
        cutoffcorr = p.cutoffcorr;
    end
end

if length(varargin) > 1
    sbxt = varargin{2};
    freq = varargin{3};
    SC = varargin{4}; % SpatialCorr
end
% lowest value to zero
% IM = setminlevel(IM);

% Get reasonable maxima sorted on height
pnt = mylocalmax(IM, border);

% Reduce number of points to 1/4 with largest peak
Nmp = round(length(pnt)*fractionmax);
dim = size(IM);
% hold on, plot(pnt(1:Nmp, 2), pnt(1:Nmp, 1), 'r+')

% Size of image sample in pixels
SzImg = VoxelSz*VoxelSz;
thx = round(threshold*SzImg); %threshold number
SigTh = round(Significance*SzImg);

[Iy, Ix] = find(ones(VoxelSz,VoxelSz));

%find a particular max
%find(pnt(:,1) == 387 & pnt(:,2) == 135)

Cnt = Savedrois.Cnt;
for i = 1:Nmp
% 	if debug
%           if i == 181
%               DISPLAY = 1;
%           end

    % Detect if this point is on a previously selected roi
    if ~(Mask(pnt(i,1),pnt(i,2)) > 0)
%         figure(3), imagesc(Mask), hold on
%         plot(pnt(i,1), pnt(i,2), '+r')
%         pause
%         else  % not already in a contour
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
            bval = 0;
            while ~bval && it < 10
                [Con, A, F] = getCon(If, th, area, PAf, py, px, Iy, Ix);
                if ~isempty(Con) %valid contour: contains point and has valid roundedness
                    if DISPLAY == 1
                        figure(2), imagesc(I), colormap gray, hold on
                        plot(px, py, '+r')
                        plot(Con.x, Con.y, 'r')
                    end
                    Con.y = Con.y + lyw - 1;
                    Con.x = Con.x + lxw - 1;
                    bval = 1; %valid maximum not in previously selected contour
                    if Cnt > 0
                        %if other points of previous rois are in this contour
                        %determine if it is simply overlapping or much
                        %larger
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
                            bval = 0; 
                        end
                    end
                    if bval % It might still be a bad roi; too many negative values below the threshold
                        pixv = I(F>0);
                        if length(find(pixv<th))/length(pixv) > 0.25 % More than 20% of the samples are below threshold
                            bval = 0;
                            th = th + thxd; % Step up contour threshold
                        end
                    end
                    if bval

                        yrange = lyw:hyw;
                        xrange = lxw:hxw;                             
                        [NwCon, NwA, NwF, NwV, Ro, Rvar] = PixelCor(size(Mask), F, sbxt, py, px, Iy, Ix, yrange, xrange, freq, cutoffcorr);
                        if ~isempty(NwCon) &&  NwA > area(1) && Ro > PAf % Exists and area greater than minimum, and roundedness > roundedness factor
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
                            
                        elseif NwA > area(1) && Ro < PAf  %not round enough but large area, => increase threshold
                            bval = 0;
                            th = th + thxd; 
                        end
                    end

                elseif A <= area(1)    % The area was too small, no valid contours
                    % but there were some closed contours,              
                    th = th - thxd; % Lower the threshold                        
                else
                    th = th + thxd; % Increase the threshold, none closed or too large
                end
                it = it + 1;
            end
%             pause(0.5); %hold off
        else
            % Ignore this maximum, below pixel significance!
%             figure(2), imagesc(I), colormap gray, hold on
%             px = pnt(i,1)-lxw+1;
%             py = pnt(i,2)-lyw+1;
%             plot(py, px, '+r')
%             notS = notS + 1;
%             disp(['Not significant : ' num2str(notS)])
        end
    else % This point was in an existing contour
%         figure(2), plot(pnt(i,2), pnt(i,1), 'or')
        In = Mask(pnt(i,1),pnt(i,2));
        % Is this point distance away from earlier point, could be a
        % separate roi (rare cases)
        if ( sqrt( (Savedrois.P(1,In)- pnt(i,1))^2 + (Savedrois.P(2,In)- pnt(i,2))^2) > distance )
            % hold on, plot(pnt(i,2), pnt(i,1), '+w')
            vx = Savedrois.Con(In).x;
            vy = Savedrois.Con(In).y;
            wx = floor(min(vx)-1):ceil(max(vx)+1);
            wy = floor(min(vy)-1):ceil(max(vy)+1);
            wx = wx(wx>0 & wx < size(IM,2));
            wy =  wy(wy>0 & wy < size(IM,1));

%             I = IM(wy,wx);
%             figure(2), imagesc(I), hold on 
            M = Mask(wy,wx);
            Mz = M; %copy
            Mz(M~=In) = 0; %set pixels of other rois to zero


            py = pnt(i,1) - wy(1)+1; %y not x
            px = pnt(i,2) - wx(1)+1;
            % plot(px, py, 'r+')
            % plot(vx-wx(1)+1, vy-wy(1)+1, 'r')

            [iy, ix] = find(ones(size(Mz))); % image indices
            [NwCon, NwA, NwF, NwV, Ro, Rvar] = PixelCor(size(Mask), Mz, sbxt, py, px, iy, ix, wy, wx, freq, cutoffcorr);
            if ~isempty(NwCon) &&  NwA > area(1) && Ro > PAf % exists and area greater than minimum
                Mz(NwF==1) = 0;
                % First original pnt
                ox = Savedrois.P(1,In) - wy(1) + 1;
                oy = Savedrois.P(2,In) - wx(1) + 1;
                % plot(ox, oy, '+r')
                [Con1, A1, F1] = getCon(Mz, 0, area, PAf, ox, oy, iy, ix);
                if  ~isempty(Con1) &&  A1 > area(1)
                    % At this point we have two valid new rois
                    % First ROI:
                    Con1.y = Con1.y + wy(1) - 1;
                    Con1.x = Con1.x + wx(1) - 1;
                    Savedrois.Con(In) = Con1;
                    Savedrois.A(In) = A1;
                    
                    Mask(Mask == In) = 0; %remove original roi
                    M(M==In) = 0;         %also from voxel mask
                    M(F1>0) = In;         %set replaced first (M cut from global Mask)

                    % Second ROI
                    Cnt = Cnt + 1;
                    NwCon.y = NwCon.y + wy(1) - 1;
                    NwCon.x = NwCon.x + wx(1) - 1;
                    Savedrois.Con(Cnt) = NwCon;
                    Savedrois.A(Cnt) = NwA;
                    Savedrois.Rvar(Cnt) = Rvar;
                    Savedrois.Roundedness(Cnt) = Ro;
                    Savedrois.P(1:3,Cnt) = pnt(i,:)';
                    M(NwF>0) = Cnt;
                    % figure(4), imagesc(M)
                    Mask(wy,wx) = M;      %add to global mask

                    SCM = SC(wy, wx);
                    SCM(NwF>0) = NwV(NwF>0); %substitute
                    SC(wy,wx) = SCM; %update correlation map
                end
            end
        end
    end
end
Savedrois.Cnt = Cnt;


