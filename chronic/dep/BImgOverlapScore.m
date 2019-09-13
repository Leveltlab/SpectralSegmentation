function [corrScore] = BImgOverlapScore(BImg)
% Takes a cell array of 2D doubles as input.
% outputs the correlation values
% of the comparisons between the different 2D doubles.
%
% If the 2D double matrixes are not the same size, the larger ones will be
% cut.
%
% Leander de Kraker
% 2018-8-3

nfiles = length(BImg);
corrScore = zeros(nfiles);
for i = 1:nfiles
    others = find((1:nfiles) < i);
    for j = others
        if all(size(BImg{i}) == size(BImg{j}))% The images are the same size :)
            corrScore(i, j) = corr(BImg{i}(:), BImg{j}(:)); % Get linear correlation value
        else
            % The images need to be cut
            minsize1 = 1:min(size(BImg{i}, 1), size(BImg{j},1));
            minsize2 = 1:min(size(BImg{i}, 2), size(BImg{j},2));
            BImgi = BImg{i}(minsize1, minsize2);
            BImgj = BImg{j}(minsize1, minsize2);
            corrScore(i, j) = corr(BImgi(:), BImgj(:));
        end
    end
end
% To save time, only the lower triangle has been calculated by the for loop
% fill upper triangle of matrix from lower triangle (output should be
% symmetrical over diagonal)
corrScore = corrScore + tril(corrScore,-1).';

end

