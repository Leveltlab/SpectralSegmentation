function [linked, linked1] = LinkOverlappingRois(inRoi, thres)
% Linking ROIs between recording pairs based on minimum overlap threshold
% 
% output:
% linked ([nfiles x nfiles] cell, with [nlinked x 5] double) in upper
% triangle
%         column 1: ROI ID of recording m, corresponds to row nr. 
%         column 2: linked ROI of other recording, corresponds to column nr.
%         column 3: percentage overlap of this ROI with other ROI
%         column 4: percentage overlap of other ROI with this ROI
%         column 5: percentage overlap for both the ROIs
% 


% linked1: column 1: ROI of mask i, 
%          column 2: linked ROI of other mask
%          column 3: number of pixels overlap
%          column 4: percentage overlap of this ROI with other ROI
nfiles = length(inRoi);
nrois = cellfun(@size, inRoi, repmat({1}, size(inRoi)));
linked1 = cell(nfiles);
nmasks = 1:nfiles;
fprintf('linking ROIs overlapping above threshold...\n')

for m = nmasks % m = current mask
    otherMasks = nmasks(nmasks~=m); % Skip own overlap
    for c = otherMasks % c = compared to
        for r = 1:nrois(m) % r = each ROI
            link = inRoi{m}{r,c}(:,3)>thres;
            if any(link)
                thisLink = [ones(sum(link),1).*r, inRoi{m}{r,c}(link,:)];
                linked1{m, c} = [linked1{m, c}; thisLink];
            end
        end
    end
end

%__________________________________________________________________________
% Creating consistent links by looking at how much % the ROIs overlap,
% viewing from both the masks, instead of 1 towards the other
linked = cell(nfiles);
fprintf('Checking whether ROI overlap is high enough between the two recordings...\n')

for m = nmasks
    otherMasks = nmasks(nmasks~=m);
    for c = otherMasks
        for i = 1:size(linked1{m,c},1)
            r = linked1{m,c}(i,1); %  roi of mask
            rc = linked1{m,c}(i,2); % r was linked with roi from other mask
            overlap1 = linked1{m,c}(i,4);
            % Overlap looking from other recording
            found = inRoi{c}{rc, m}(inRoi{c}{rc, m}(:,1)==r, :); 
            
            overlap2 = found(3);
            overlap = (overlap1 + overlap2) ./ 2;
            if overlap > thres % true if mutual overlap is large enough!
                linked{m,c} = [linked{m,c}; [r, rc, overlap1, overlap2, overlap]];
            end
        end
    end
end

% Now combine the linked ROIs from recording 1-> recording 2 and 2->1 etc
% to one array
for c = nmasks
    otherMasks = c+1:nmasks(end);
    for m = otherMasks
        % Switch columns to match the order from other recording comparison
        if ~isempty(linked{m,c})
            linked{m,c} = [linked{m,c}(:,2),linked{m,c}(:,1),linked{m,c}(:,4),linked{m,c}(:,3),linked{m,c}(:,5)];
            linked{c,m} = [linked{c,m}; linked{m,c}];
            linked{m,c} = [];
            [~, sorted]  = sort(linked{c,m}(:,1));
            linked{c,m} = linked{c,m}(sorted,:);
            linked{c,m} = unique(linked{c,m},'rows');
        end
    end
end
