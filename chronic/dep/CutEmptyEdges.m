function [BImg2, Masks2, PP] = CutEmptyEdges(BImg2, Masks2, PP, nfiles)
%  [BImg2, Masks2, PP] = CutEmptyEdges(BImg2, Masks2, PP)
% Cuts away edges which only contain zeros in the BImg2 from the BImg and
% corresponding Masks, and also edits the contour information accordingly.
% 
% This is desireable because consecutive registration rounds can cause
% empty edges around the images.
% 
% Leander de Kraker
% 2021

BImgMat = cat(3, BImg2{:});
BImgMat = max(BImgMat,[],3);
zeroCols = find(all(BImgMat'==0));
zeroRows = find(all(BImgMat==0));

if ~isempty(zeroCols)
    % Cut at the end (bottom)
    if ismember(size(BImgMat,1), zeroCols) % Does the end have columns with only zeros?
        skipper = find(diff(zeroCols)>1);
        if ~isempty(skipper)
            % Black columns at multiple different places in the BImgMat?!
            if length(skipper)>2
                warning('contact Leander'); 
                skipper = skipper(end); 
            end
            skipper = skipper + 1;
        else
            skipper = 1;
        end
        for i = 1:nfiles
            BImg2{i}(zeroCols(skipper):end,:) = [];
            Masks2{i}(zeroCols(skipper):end,:) = [];
        end
    end
    % Cut at the beginning (top) (also change PP)
    if ismember(1, zeroCols)
        skipper = find(diff(zeroCols)>1);
        if ~isempty(skipper)
            % Black columns at multiple different places in the BImgMat?!
            if length(skipper)>2
                warning('contact Leander');
                skipper = skipper(1); % Only delete until the first non-zero row
            end
        else
            skipper = length(zeroCols); 
        end
        for i = 1:nfiles
            BImg2{i}(1:zeroCols(skipper),:) = [];
            Masks2{i}(1:zeroCols(skipper),:) = [];
            PP(i).P(1,:) = PP(i).P(1,:) - skipper;
            for j = 1:PP(i).Cnt
                PP(i).Con(j).x = PP(i).Con(j).x - skipper;
            end
        end
    end
end
if ~isempty(zeroRows)
    % Cut at the end (right)
    if ismember(size(BImgMat,2), zeroRows) % Does the end have rows with only zeros?
        skipper = find(diff(zeroRows)>1);
        if ~isempty(skipper)
            % Black rows at multiple different places in the BImgMat?!
            if length(skipper)>2
                warning('contact Leander'); 
                skipper = skipper(end); 
            end
            skipper = skipper + 1;
        else
            skipper = 1;
        end
        for i = 1:nfiles
            BImg2{i}(:,zeroRows(skipper):end) = [];
            Masks2{i}(:,zeroRows(skipper):end) = [];
        end
    end
    % Cut at the beginning (left) (also change PP)
    if ismember(1, zeroRows)
        skipper = find(diff(zeroRows)>1);
        if ~isempty(skipper)
            % Black columns at multiple different places in the BImgMat?!
            if length(skipper)>2
                warning('contact Leander'); 
                skipper = skipper(1); % Only delete until the first non-zero row
            end
        else
            skipper = length(zeroRows); 
        end
        for i = 1:nfiles
            BImg2{i}(:,1:zeroRows(skipper)) = [];
            Masks2{i}(:,1:zeroRows(skipper)) = [];
            PP(i).P(2,:) = PP(i).P(2,:) - skipper;
            for j = 1:PP(i).Cnt
                PP(i).Con(j).y = PP(i).Con(j).y - skipper;
            end
        end
    end
end

