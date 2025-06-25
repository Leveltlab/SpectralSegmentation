function masks = getMask(filenames, filepaths)
    arguments (Input)
        filenames (:, 1) string
        filepaths (:, 1) string
    end
    arguments (Output)
        masks (:, 1) cell
    end
    fetchStr = 'Mask';

    masks = cell(length(filenames), 1);
    rawMaskTable = fetchData(filenames, filepaths, {fetchStr});
    
    for i = 1:length(filenames)
        mask = rawMaskTable.(fetchStr){i};
        mask(rawMaskTable.(fetchStr){i} > 0) = 1;
        masks{i} = mask;
    end
end