function [filenames, filepaths] = fetchFile(filter)
    arguments (Input)
        filter (:, 1) string = '';
    end
    arguments (Output)
        filenames (:, 1) string
        filepaths (:, 1) string
    end
    
    [rawNames, rawPaths] = uigetfile(filter, "MultiSelect","on");
    
    if ischar(rawNames)
        nFiles = 1;
    elseif iscell(rawNames)
        nFiles = length(rawNames);
    elseif rawNames == 0
        nFiles = 0;
    else
        error("error fetching file, returned filename is not 0, a cell or a char")

    end

    if ~ischar(rawPaths) && rawPaths ~= 0
        error("error fetching file, returned filepath is not a char as expected")
    end

    filenames = strings(nFiles, 1);
    filepaths = strings(nFiles, 1);

    for i=1:nFiles
        if ischar(rawNames)
            filenames(i) = rawNames;
        else
            filenames(i) = rawNames{i};
        end
        filepaths(i) = rawPaths;
    end
end