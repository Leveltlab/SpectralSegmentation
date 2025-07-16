function h = PlotChronicConnectedROIs(linkMat, PP, colors)
% Plot Connection lines between matched ROIs
% 
% Input: 
%   linkMat ([nmatches x nfiles] 2D double)
%   PP (struct): contour & ROI information 
% 
% Output:
%   h (handles to every line)
% 
% Leander de Kraker
% 2023-10-26
% 

if isempty(colors)
    colors = flipud(cmapL([0 0 1; 0 1 1; 0 1 0; 1 1 0; 1 0 0; 1 0 1], size(linkMat, 2)));
end

nLinks = sum(linkMat~=0, 2);
n = nLinks-1;
n = (n.*(n+1))./2;

h = gobjects(sum(n), 1);
counter = 0;
for i = 1:size(linkMat, 1)
    idx = linkMat(i,:);
    present = find(idx);
    for j = present
        presentRest = present(present>j);
        for k = presentRest
            counter = counter + 1;
            cordinatesX = [PP(j).P(1, idx(j)), PP(k).P(1, idx(k))];
            cordinatesY = [PP(j).P(2, idx(j)), PP(k).P(2, idx(k))];
            colorkj = mean(colors([j k], :), 1);
            h(counter) = line(cordinatesX, cordinatesY, 'color', colorkj);
%             hCoMLines(counter) = annotation('arrow', cordinatesX, cordinatesY, 'color', colorkj);
        end
    end
end