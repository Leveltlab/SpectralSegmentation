function img = CutEmptyEdgesImg(img, cutting)
%  img = CutEmptyEdgesImg(img)
% Cuts away edges as denoted by the chronic registration
% 
% input:
%   img (cell array with images of same size)
%   cutting ([4 x 1] double): 
% 
% output:
%   img (cell array with images of same size, cropped)
% 
%
% Leander de Kraker
% 2022-7-11
%

nfiles = length(img);

if ~isnan(cutting(1)) 
    % Cut at the end (bottom)
    for i = 1:nfiles
        img{i}(cutting(1):end,:,:) = [];
    end
end
if  ~isnan(cutting(2))
    % Cut at the beginning (top)
    for i = 1:nfiles
        img{i}(1:cutting(2),:,:) = [];
    end
end
if ~isnan(cutting(3)) 
    % Cut at the end (right)
    for i = 1:nfiles
        img{i}(:, cutting(3):end,:) = [];
    end
end
if ~isnan(cutting(4))
    % Cut at the beginning (left)
    for i = 1:nfiles
        img{i}(:, 1:cutting(4),:) = [];
    end
end