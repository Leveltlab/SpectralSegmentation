function PPfromMaskFile(varargin)
% Find contours around 2D ROI mask which is inside a file. Save those
% contours !and potentially edited ROI mask! into the file
% 
% This script has the ability to de-island ROIs:
% In case the mask has island ROIs (non- touching parts of ROI), the ROI
% get split into multiple parts/ROIs. !This edits the Mask!
% 
% ROIs that are too small are rejected. This also edits the Mask!
% 
% 
% PPfromMask; % Get pop up to select file
% PPfromMask('folderspath\', 'mySPSIGfile.mat')
% PPfromMask(filepath, filename, deIsland, minSize)
% 
% Input (all option, can be left empty):
%   - filepath (string).
%   - filename (string).
%   - deIsland (boolean). default = true. De-island unconnected ROI.
%   - minSize  (digid). default = 10. Below which size to reject ROI.
% 
% 
% If you already have fluorescence signal for the ROIs, and ROIs get
% deleted or split during this process, do not forget update the sig by
% running the last section!!
% 
% See also: PPfromMask, to create contours from workspace variables
%           PPModernize, to calculate the other properties of PP!
%           Mask3D_Flatten, to create a 2D mask from a 3D mask.
% 
% Leander de Kraker
% 2020-6-3, 2020-6-10, 2022-10-18
% 
% 

if exist('varargin', 'var') && nargin >= 2 && ~isempty(varargin{2})
    filepath = varargin{1};
    if ~isempty(filepath)
        if ~strcmp(filepath(end), '\')
            filepath = [filepath '\']; 
        end
    end
    filename = varargin{2};
else
    [filename, filepath] = uigetfile('*.mat', 'Select the file which has a Mask but no good PP');
end

load([filepath filename], 'PP', 'Mask', 'BImg')
    
if exist('varargin', 'var') && nargin >= 3 && ~isempty(varargin{3})
    deIsland = varargin{3};
else
    deIsland = true; % turn de-islanding of ROIs on or off. true | false
end

if exist('varargin', 'var') && nargin >= 4 && ~isempty(varargin{4})
    minSize = varargin{4};
else
    minSize = 10; % ROIs aren't allowed to be smaller then this many pixels
end

nROIs = max(Mask(:));
nROIs2 = length(unique(Mask(:)))-1;
if nROIs ~= nROIs2
    warning('The mask is the problem!! missing ROIs at some places\n')
    find(diff(unique(Mask(:)))>1)
    return
end

if exist('PP', 'var')
    if nROIs ~= PP.Cnt
        fprintf('PP has different amount of ROIs than mask: %d vs %d\n',...
                    nROIs, PP.Cnt)
        fprintf('replacing PP')
        clear PP
    end
end

PP.Con = struct('x', [], 'y', []);
PP.P = [];
PP.A = [];
minval = min(BImg(:));

maskX = repmat(1:size(Mask,2),[size(Mask,1), 1]);
maskY = repmat((1:size(Mask,1))', [1, size(Mask,2)]);
roisUpdatedIdx_ByPPfromMask = 1:max(Mask(:));
splitted = 0;
i = 1;
while i <= max(Mask(:)) % mask nr ROIs can grow because of de-islanding ROIs
    maski = double(Mask==i);

    if sum(maski(:)) < minSize % if ROI is too small, delete it and try again (next)
        
        fprintf('ROI %d too small, deleting it\n', i)
%         figure; imagesc(maski);
%         hold on; plot(PP.P(1,1:i), PP.P(2,1:i),'xr');
%         title(sprintf('ROI %d too small, %d pixels', i, sum(maski(:))))
        Mask(Mask==i) = 0;
        Mask(Mask>i) = Mask(Mask>i)-1;
        roisUpdatedIdx_ByPPfromMask(end) = [];
        roisUpdatedIdx_ByPPfromMask(i:end) = roisUpdatedIdx_ByPPfromMask(i:end)+1;
        
    else % surely a contour will be found
        coni = contourc(maski, [0.5 0.5]);
        
        
        if mod(i, 10)== 1
            fprintf('%d/%d\n', i, max(Mask(:))); 
            hold off
            imagesc(maski)
            hold on
            plot(coni(1,:), coni(2,:),'r')
            drawnow
            pause(0.01)
        end
        
        n = coni(2,1); % number of coordinates in this contour
        if deIsland && (size(coni,2) > (n+1)) % if the size is longer, there are multiple contours in this contour
            % DE-ISLANDING
            fprintf('ROI %d has multiple contours. Splitting ROI\n', i)
            PP.Con(i).x = coni(1,2:n); % get coordinates of this one contour
            PP.Con(i).y = coni(2,2:n);

            inCon = inpolygon(maskX, maskY, PP.Con(i).x, PP.Con(i).y);
            if sum(inCon(:))~=0
                % Because the Mask gets edited, the next iteration will be what was the next part of this ROI
                Mask(Mask>=i) = Mask(Mask>=i) + 1; % Add a gap between previous and this & all other ROIs
                Mask(inCon) = i; % That's this contour! 
                roisUpdatedIdx_ByPPfromMask(i+1:end+1) = roisUpdatedIdx_ByPPfromMask([i:end]);
                splitted = splitted + 1;
            else
                 % pixels inside polygon detection went wrong, abort islanding
                 fprintf('in polygon went wrong, keeping original ROI, not de-islanding')
            end
            
            
        elseif ~deIsland && (size(coni,2) > (n+1))
            % Retrieve at which indexes is the information about coordinates stored
            pos = 1;
            posi = 1;
            while pos < size(coni, 2) 
                posi = posi + 1;
                pos(posi) = coni(2,pos(end))+1+pos(end);
            end
            pos(end) = [];
            % Remove those info points of the contour, and nan them to make them invisible during plotting
            coni(:, pos) = nan;
            coni(:, 1) = [];
            PP.Con(i).x = coni(1,:);
            PP.Con(i).y = coni(2,:);
        else
            coni = getcontourlines(coni);
            PP.Con(i).x = coni.x;
            PP.Con(i).y = coni.y;
        end
        
        % Find the maximum in the ROI, to find 'seedpoint' location
        BImgi = BImg;
        BImgi(maski==0) = minval;
        [~, PP.P(1,i)] = max(max(BImgi));
        [~, PP.P(2,i)] = max(BImgi(:,PP.P(1,i)));
        PP.A(i) = sum(maski(:));
        
        i = i + 1;
    end
end

fprintf('%d ROIs were edited/ de-islanded in the Mask\n', splitted)

nROIs = max(Mask(:));
nCons = length(PP.Con);
if nCons  ~= nROIs
    fprintf('Somehow the number of ROIs got fucked up between Mask and PP\n')
end
PP.Cnt = nROIs;

figure;
maski = Mask;
imagesc(maski)
hold on
for i = 1:PP.Cnt
    plot(PP.Con(i).x, PP.Con(i).y, 'r','Linewidth',2)
end
plot(PP.P(1,:),PP.P(2,:), 'xg')


%% Save the change

save([filepath, filename], 'PP', 'Mask', 'roisUpdatedIdx_ByPPfromMask', '-append')
fprintf('\nsaved PP and potentially updated Mask too %s\n\n',filename)


%% If you had signals for the Mask and ROIs have been deleted, the signal 
% needs to be updated.
% Call the signal variable sig and run this section
% sig = [nroisOld x time] -> [nroisUpdated x time]

% sig = sig(roisUpdatedIdx_ByPPfromMask, :);
% feel free to rename this variable and save it to the file

