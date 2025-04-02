function score = ExtendScore(score, nmatches)
% extend scores with 0: there is no overlap for single neurons
score(end+1:nmatches) = 0;