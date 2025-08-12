function linkMat = SquishLinkedRois(linked)
% Go from linked ROIs between every recording pair to linked ROIs matrix
% 
% output: 
%   linkMat ([nMatches x nrec] double): ROI ID of matched ROIs
% 

nfiles = size(linked, 2);
nmasks = 1:nfiles;

% all links in a matrix
fprintf('Compressing link data into match matrix...\n')
linkMatV1 = cell(nfiles-1, 1);
for m = 1:nfiles-1
    otherMasks = m+1:nmasks(end);
    linkMatV1{m} = zeros(1,nfiles);
    for c = otherMasks
        if ~isempty(linked{m,c})
            own = linked{m,c}(:,1);
            compared = linked{m,c}(:,2);
            overlap = linked{m,c}(:,5);
            linkMatV1{m}(end+1:end+length(own), [m, c]) = [own, compared];
        end
    end
    linkMatV1{m}(1,:) = [];
    % Sort to the neuron of current m
    [~, idx] = sort(linkMatV1{m}(:,m));
    linkMatV1{m} = linkMatV1{m}(idx,:);
end

% Compress the matrix m so the neuron of mask m only appears once
linkMat = cell(nfiles-1, 1);
for m = 1:nfiles-1
    [rois, limits] = unique(linkMatV1{m}(:,m));
    limits(end+1) = size(linkMatV1{m},1)+1;
    linkMat{m} = zeros(length(rois), nfiles);
    for j = 1:length(rois) % for each roi. put others in one row
        linkMat{m}(j,:) = max(linkMatV1{m}(limits(j):limits(j+1)-1,:),[],1);
    end
end
linkMat = cell2mat(linkMat);

% Sort & squish the matrix
for i = 2:nfiles
    [~,sorted] = sort(linkMat(:,i));
    linkMat = linkMat(sorted,:);
    [rois, limits] = unique(linkMat(:,i));
    if length(limits) > 1
        limits(end+1) = size(linkMat,1)+1;
        dif = diff(limits); % number of elements with ROI j
        % If one of the recordings matched another ROI keep it. So even if
        % the matched ROI is not found via some of the recordings, it can still
        % get linked with those recordings.
        for j = 2:length(rois)
            if dif(j) > 1
                rows = limits(j):(limits(j)+dif(j)-1);
                linkMat(limits(j),:) = max(linkMat(rows,:));
            end
        end
        % Now we made sure we don't lose any links, squish the matrix
        linkMat = linkMat([1:limits(2), limits(3:end-1)'],:);
    end
end