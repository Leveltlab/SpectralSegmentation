function imgs = MoveLikeChronicFun(imgs, tr)
% Apply transformations to a cell array of pictures
% 
%
% Leander de Kraker
% 2022-7-11


%% Load the background images

nfiles = length(imgs);

for i = 1:length(tr)
    tri = tr(i);
    
    % Apply offsets
    regiVert = tri.vert;
    regiHori = tri.hori;
    imgNew = cell(1, nfiles);
    for f = 1:nfiles
        imgNew{f} = zeros(max(regiVert(:,2)), max(regiHori(:,2)), size(imgs{f},3));
        imgNew{f}(regiVert(f,1):regiVert(f,2), regiHori(f,1):regiHori(f,2), :) = imgs{f};   
    end
    
    % Cut empty edges away
    imgs = CutEmptyEdgesImg(imgNew, tri.cutting);

    for f = find(abs(tri.rotation)>0.02)' % If rotation angle < 0.02, don't rotate
        imgs{f} = imrotate(imgs{f}, tri.rotation(f), 'bilinear', 'crop');
    end
    
end

