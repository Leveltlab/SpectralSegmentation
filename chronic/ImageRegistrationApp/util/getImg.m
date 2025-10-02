function images = getImg(filenames, filepaths)
    arguments (Input)
        filenames (:, 1) string
        filepaths (:, 1) string
    end
    arguments (Output)
        images (:, 1) cell
    end
    fetchStr = 'BImg';
    
    images = cell(length(filenames), 1);
    rawImagesTable = fetchData(filenames, filepaths, {fetchStr});
    
    for i = 1:length(filenames)
        images{i} = cleanImage(rawImagesTable.(fetchStr){i});
    end
    
    function cleanImg = cleanImage(Img)
        isClean = false;
        cleanImg = Img;
        
        while not (isClean) %trim outer image border
            isClean = true;
            if isequal(cleanImg(1,:),cleanImg(end,:))
                cleanImg = cleanImg(2:end-1,:);
                isClean = false;
            end
            if isequal(cleanImg(:,1),cleanImg(:,end))
                cleanImg = cleanImg(:,2:end-1);
                isClean = false;
            end
        end
        infMask = isinf(cleanImg); % create mask, 1 if value at index is inf, 0 if not
        [minVal, maxVal] = bounds(cleanImg(~infMask)); % check the bounds of the values not at inf
        cleanImg(infMask) = minVal; % clean inf values with true min vals
        cleanImg = (cleanImg - minVal) ./ range([minVal, maxVal]); % normalise image to have range [0, .., 1]
    end
end


