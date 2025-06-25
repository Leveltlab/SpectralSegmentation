function [filenames, filepaths] = fetchDir()
    arguments (Output)
        filenames (:, 1) string
        filepaths (:, 1) string
    end
    
    [rawNames, rawPaths] = selectFolder();
    nFiles = length(rawNames);
    filenames = strings(nFiles, 1);
    filepaths = strings(nFiles, 1);
    
    for i=1:nFiles
        filenames(i) = rawNames{i};
        filepaths(i) = rawPaths{i};
    end

    function [rawNames, rawPaths] = selectFolder()
        folderPath = uigetdir();
        if isequal(folderPath, 0)
            rawNames = {};
            rawPaths = {};
            
        else
            tempFolder = dir(folderPath);
            tempFolder([tempFolder.isdir]) = []; 
            
            rawNames = {tempFolder.name}';
            rawPaths = {tempFolder.folder}';
        end
    end
end