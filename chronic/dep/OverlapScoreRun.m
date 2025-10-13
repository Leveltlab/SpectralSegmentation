function [score, h] = OverlapScoreRun(linkMat,inRoi, nLinks, thres, doPlot)
% The 'OVERLAP SCORE' is based on the pixel overlap between the linked ROIs. 
% This also takes into account the ROIs which are not officially linked with
% each other, but were linked together because of mutual links.
arguments
    linkMat
    inRoi
    nLinks double = sum(linkMat~=0,2)
    thres (1,1) double = NaN
    doPlot (1,1) boolean = true
end

nMatches = size(linkMat, 1);
nfiles   = size(linkMat, 2);
score = zeros(nMatches,1);
for i = 1:nMatches
    score(i) = OverlapScore(inRoi, linkMat, i);
end

if doPlot % Plot scores of the ROIs
    h = figure; hold on
    plot(nLinks+(rand(length(nLinks),1)./5-0.1), score, 'o')
    % plot(nLinks+(score2-0.5), score2, 'o')
    h.Children.XTick = 2:nfiles;
    ylim([0 1])
    title(sprintf('overlap scores for chronic matched ROIs. thr=%.2f',thres))
    xlabel('ROI found in n recordings')
    ylabel('score (% pixel overlap average)')
    % Average score per number of links in the match
    scoreMean = zeros(1,nfiles-1);
    for i = 2:nfiles
        scoreMean(i-1) = mean(score(nLinks==i));
    end
    plot(2:nfiles, scoreMean,'ro-')
    line([1.5 nfiles+0.5],[thres, thres],'color',[0 0 0 0.4],'linestyle','--')
    legend({'match score','average score','threshold'},...
        'location','southeast')
end