function varargout = CorrRoiCorners(corner, check, Mask, sbxt, dim)
% This function correlates signals of each pixel from the ROIs in check
% with the signals of the 'corner' coordinates which are part of the
% same ROI.
%
% output: - corrImg: cells with 2D arrays, each cell containing the 
%               correlation coefficients resulting from a corner signal with
%               all pixels of the ROIs
%         - pieceSigs: the fluoresence signals of the corners
%         - corrs: the correlation coefficients for each ROI, with each
%               corner (columns) correlated to all pixels in the ROI (rows)
%         - corrsSub: The correlation coefficients with the correlation  
%               coefficients from other corners subtracted, so similar
%               correlations result in values near 0, instead of positive
%               correlation. (important for clustering)
%         - pos: struct with the pixel coordinates of the ROI (1D, x and y)
%           
 

nCheck = length(check);
if nCheck > 10
    f = waitbar(0,'Initializing...');
    tic
end

Mt = Mask'; %the mask needs to be transposed
indices = find(ismember(Mt(:), check));
% indices = indices + double(dim(2)); %add a line?! OR NO LINE ?!?!!

Sigrois = sbxt.Data.y(:,indices);
Sigrois = double(Sigrois);
% sig2 = cell(nCheck,1); % This sig keeps the signal of every pixel in all ROI
pos = struct('d1', [], 'x', [], 'y', []);
pos(nCheck).x = [];

ncorners = size(corner,2);

pieceSigs = cell(nCheck,1);
corrs = cell(nCheck,1);
corrsSub = cell(nCheck,1);
corrImg = zeros(numel(Mt), ncorners);

for i = 1:nCheck
    % Retrieve signal and get 1D indexes to the ROIs
    pos(i).d1 = find(Mt(:)== check(i));
    [~, ia] = intersect(indices, pos(i).d1); % get indices of to this roi
    sig2 = Sigrois(:,ia); % get data from this selection
    
    % Get the mean signal from 9 pixels around the corner coordinates
    pieceSigs{i} = zeros(ncorners,dim(1));
    for j = 1:ncorners
        xmesh = reshape(meshgrid(corner(1,j,i)-1:corner(1,j,i)+1), [9,1]);
        ymesh = reshape(meshgrid(corner(2,j,i)-1:corner(2,j,i)+1)', [9,1]);
        corner1D = (ymesh-1).*dim(2)+xmesh; % convert to 1D
        corner1D = ismember(pos(i).d1, corner1D); % find the actual ROI pixels
        
        pieceSigs{i}(j,:) = mean(sig2(:, corner1D),2);
    end
    
    % Correlate the four 'start-correlation' with rest of the ROI
    corrs{i} = corr(double(sig2), double(pieceSigs{i})',...
        'rows','pairwise');
    
    % Subtract correlation of others from each other, so positive
    % correlations are cancelled out.
    corrsSub{i} = zeros(size(corrs{i}));
    for j = 1:ncorners
        others = 1:ncorners~=j;
        corrsSub{i}(:,j) = corrs{i}(:,j) - mean(corrs{i}(:,others),2);
    end
    
    corrImg(pos(i).d1,:) = corrs{i};
    
    % Convert 1D indexes to 2D
    pos(i).x = mod(pos(i).d1, dim(2));
    pos(i).y = floor(pos(i).d1 / dim(2)) + 1;
    
    if mod(i,11)==0
        waitbar(i/(nCheck*1.1),f,'Calculating...');
    end
end

corrImg = reshape(corrImg, [size(Mt), ncorners]);


varargout = {corrImg, pieceSigs, corrs, corrsSub, pos};

if nCheck>10
    close(f)
    toc
end
end

