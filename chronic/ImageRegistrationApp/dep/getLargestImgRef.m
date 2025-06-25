function finalImgRef = getLargestImgRef(ImgRefs)
    arguments (Input)
        ImgRefs (:, 1) imref2d
    end
    arguments (Output)
        finalImgRef imref2d
    end
    %first retrieve the world limit property into a 2*n matirx, then
    %reshape to an array of 1*2n. 
    XWorldLimitsArray = reshape(vertcat(ImgRefs.XWorldLimits), 1, []);
    YWorldLimitsArray = reshape(vertcat(ImgRefs.YWorldLimits), 1, []);
    %retrieve the range in coords to find the smallest and largest coord
    %needed to fit everything, round it to remove small fraction errors
    [xWorldmin, xWorldmax] = bounds(round(XWorldLimitsArray));
    [yWorldmin, yWorldmax] = bounds(round(YWorldLimitsArray));
    %calculate size
    xSize = xWorldmax - xWorldmin;
    ySize = yWorldmax - yWorldmin;
    finalImgRef = imref2d([ySize xSize], [xWorldmin, xWorldmax], [yWorldmin, yWorldmax]);
    if finalImgRef.PixelExtentInWorldX ~= 1 || finalImgRef.PixelExtentInWorldY ~= 1
        error('the calculated Pixel extent is wrong, this is a bug');
    end
end