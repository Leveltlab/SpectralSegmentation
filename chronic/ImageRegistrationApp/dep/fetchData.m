function outputTable = fetchData(filenames, filepaths, vars)
    arguments (Input)
        filenames (:,1) string
        filepaths (:,1) string
        vars (1, :) cell
    end
    arguments (Output)
        outputTable (:, :) table
    end
    
    nFiles = length(filenames);
    for i = 1:nFiles
        filename = fullfile(filepaths(i), filenames(i));
        loadedStruct(i) = load(filename, vars{:});
    end
    fields = fieldnames(loadedStruct);
    
    if length(fields) == 1 && nFiles == 1
        loadedStruct.(fields{1}) = {loadedStruct.(fields{1})};
    end
    
    outputTable = struct2table(loadedStruct);
end
