function pnt = localmax(IM, edge, vox)
% retrieve local maxima in image after substracting median filtered image
%
% IM; the spectral image
% edge; border of image to ignore
% vox; voxel size to process
% Pxd; number of pixels of separation between maxima
% Chris van der Togt, 2021
% Netherlands Institute for Neuroscience 

    Pxd = 7;
    
    % median filter
    J = ordfilt2(IM, vox^2/2, ones(vox));
    I = IM - J + mean(J(:), 'omitnan');
    
    %remove border
    dsz = size(I);
    d = I(edge:dsz(1)-edge,edge:dsz(2)-edge);
    
    % max filter
    B = ordfilt2(d, Pxd^2, ones(Pxd));
    peaks = (d >= B);

    %get maxima
    v = d(peaks); %values
    [i, j] = find(peaks); %coordinates
    x = i + edge-1;
    y = j + edge-1;
    pnt = [x, y, v];
       
    %order by magnitude
    [~, ix] = sort(pnt(:,3), 'descend');
    pnt = pnt(ix,:);
    
    

            
            