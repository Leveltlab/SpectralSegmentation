% Find contours around 2D ROI mask
%
% This script has the ability to de-island ROIs:
% In case the mask has island ROIs (non- touching parts of ROI), the ROI
% get split into multiple parts/ROIs!
% This de-islanding can be turned on or off
% 
% If you already have fluorescence signal for the ROIs, and ROIs get
% deleted or split during this process, do not forget update the sig by
% running the last section!!
% 
% Leander de Kraker
% 2020-6-3, 2020-6-10
% 
% 

[fn, pn] = uigetfile('*.mat', 'Select the problematic file');
load([pn fn], 'PP', 'Mask', 'BImg')


deIsland = false; % turn de-islanding of ROIs on or off. true | false
minSize = 10; % ROIs aren't allowed to be smaller then this many pixels


nROIs = max(Mask(:));
nROIs2 = length(unique(Mask(:)))-1;
if nROIs ~= nROIs2
    warning('The mask is the problem!! missing ROIs at some places\n')
    find(diff(unique(Mask(:)))>1)
    return
end

if exist('PP', 'var')
    PPold = PP;
    PPexisted = true;
    if nROIs ~= PP.Cnt
        fprintf('PP has different amount of ROIs then mask: %d vs %d\n',...
                    nROIs, PP.Cnt)
    end
else
    PPexisted = false;
end

clear PP
PP.Con.x = [];
PP.Con.y = [];
PP.P = [];
minval = min(BImg(:));

maskX = repmat(1:size(Mask,2),[size(Mask,1), 1]);
maskY = repmat((1:size(Mask,1))', [1, size(Mask,2)]);
roisUpdatedIdx_ByPPfromMask = 1:max(Mask(:));
splitted = 0;
i = 1;
while i <= max(Mask(:)) % mask nr ROIs can grow because of de-islanding ROIs
    maski = double(Mask==i);
%     [~, ~, maskib] = BufferMask(maski, 1);
%     maski = maski+maskib;
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
        
        % Find the maximum in the ROI
        BImgi = BImg;
        BImgi(maski==0) = minval;
        [~, PP.P(1,i)] = max(max(BImgi));
        [PP.P(3,i), PP.P(2,i)] = max(BImgi(:,PP.P(1,i)));

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
% imagesc(CreateRGB2({BImg, maski},[0 0 1; 1 1 0], true, true))
hold on
% colors = flipud(parula(nROIs));
for i = 1:PP.Cnt
    plot(PP.Con(i).x, PP.Con(i).y, 'r','Linewidth',2)
end
% colormap()
plot(PP.P(1,:),PP.P(2,:), 'xg')


%% Save the change

save([pn, fn], 'PP', 'Mask', 'roisUpdatedIdx_ByPPfromMask', '-append')
fprintf('\nsaved PP and potentially updated Mask to %s\n\n',fn)
clearvars maski i BImgi coni maskX maskY minval n nROIs2 inCon

%% Also edit the rest of PP!!
PP = PPModernize(pn, fn);

%% If you had signals for the Mask and ROIs have been deleted, the signal 
% needs to be updated.
% Call the signal variable sig and run this section
% sig = [nroisOld x time] -> [nroisUpdated x time]

% sig = sig(roisUpdatedIdx_ByPPfromMask, :);
% feel free to rename this variable and save it to the file

