%{ 
Use on startup to add relevant SpectralSeg files to the matlab Path

In the below excludedPaths variable can add folder extensions to ignore
%}
excludedPaths = {'\.git', '\docs'};
currentDir = pwd;

excludePattern = cellfun(@(str) fullfile(currentDir, str), excludedPaths, 'UniformOutput', false);


allPaths = genpath(currentDir);
allPaths = strsplit(allPaths, pathsep);
allPaths = allPaths(~cellfun('isempty', allPaths));
keepPaths = true(size(allPaths));

for i = 1: length(excludedPaths)
    keepPaths = keepPaths & ~contains(allPaths, excludedPaths{i});
end

filteredPaths = allPaths(keepPaths);

for i = 1: length(filteredPaths)
    addpath(filteredPaths{i});
end

disp('Paths added to MATLAB:');
disp(filteredPaths');

clear