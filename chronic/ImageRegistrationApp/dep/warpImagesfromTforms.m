function warpedImg = warpImagesfromTforms(transforms, rawImages)
    n_imgs = max(size(rawImages));
    outView = @(x) affineOutputView(size(rawImages{x}),transforms(x),"BoundsStyle","FollowOutput");
    outviews = arrayfun(outView, 1:n_imgs)';
    
    finalView = getLargestImgRef(outviews);
    warpedImg = cell(size(rawImages));
    for i = 1:n_imgs
        warpedImg{i} = imwarp(rawImages{i},transforms(i), "OutputView", finalView, "interp", "nearest");
    end
end