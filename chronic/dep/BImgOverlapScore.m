function [corrScore, calculatedBool] = BImgOverlapScore(imgs, toCorr)
% Calculate correlation between pixelvalues of multiple images
% If the 2D double matrixes are not the same size, the larger ones will be
% cut.
% 
% [corrScore, calculatedBool] = BImgOverlapScore(imgs, toCorr)
% corrScore = BImgOverlapScore(imgs)
%
% Input: 
%   - imgs ([n x 1] cell array with doubles inside)
%   - (optional) toCorr ([m x 1] double): For which images to calculate
%               correlation to all the other images to
%
% Output:
%   - corrScore ([n x n] double): correlation between each image pair
%   - calculatedBool ([n x n] logical): which value got calculated, based
%                                       on toCorr
% 
% Leander de Kraker
% 2018-8-3
% 2025-11-20: Changed to be able to ignore nans, and added ability to
%             select which images to calculate correlations for
% 
arguments (Input)
    imgs
    toCorr = 1:length(imgs)
end

if size(toCorr, 1)>1 && size(toCorr,2)==1
    toCorr = toCorr';
elseif all(size(toCorr)>1) 
    warning('toCorr variable has unexpected size')
end

nfiles = length(imgs);
% Calculate once if there are nans in the images
nanMasks = cell(nfiles, 1);
for i = 1:nfiles
    nanMasks{i} = isnan(imgs{i});
end

corrScore = zeros(nfiles);
for i = toCorr
    others = find((1:nfiles) > i);
    for j = others
        if all(size(imgs{i}) == size(imgs{j}))% The images are the same size :)
            nanPresent = nanMasks{i}(:) | nanMasks{j}(:);
            if any(nanPresent)
                corrScore(i, j) = corr(imgs{i}(~nanPresent), imgs{j}(~nanPresent));
            else
                corrScore(i, j) = corr(imgs{i}(:), imgs{j}(:)); % Get linear correlation value
            end
        else
            % The images need to be cut
            minsize1 = 1:min(size(imgs{i}, 1), size(imgs{j},1));
            minsize2 = 1:min(size(imgs{i}, 2), size(imgs{j},2));
            BImgi = imgs{i}(minsize1, minsize2);
            BImgj = imgs{j}(minsize1, minsize2);
            nanPresent = nanMasks{i}(minsize1, minsize2) | nanMasks{j}(minsize1, minsize2);
            if any(nanPresent(:))
                corrScore(i, j) = corr(BImgi(~nanPresent(:)), BImgj(~nanPresent(:)));
            else
                corrScore(i, j) = corr(BImgi(:), BImgj(:));
            end
        end
    end
end
% To save time, only the lower triangle has been calculated by the for loop
% fill upper triangle of matrix from lower triangle (output should be
% symmetrical over diagonal)
corrScore = corrScore + triu(corrScore,-1).';

calculatedBool = false(nfiles);
calculatedBool(toCorr, :) = true;
calculatedBool(:, toCorr) = true;

