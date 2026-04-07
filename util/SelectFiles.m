function [files, nfiles]= SelectFiles()
% [files, nfiles] = SelectFiles();
% 
% select as many files as the user wants via pop-up dialog. Multiselect is
% possible per individual selection. Press cancel to stop selecting
% 
% output:
%   - files (table). column 1=path, column 2=name 
%   - nfiles(scalar integer). Number of files selected.
% Leander de Kraker
% 2026-3-13
% 
 
%%
varNames = {'path', 'name'};
files = table('Size', [0, 2], 'VariableTypes', {'string', 'string'}, ...
                              'VariableNames', varNames);
selecting = true;
i = 1;
while selecting
    if i > 1
        prompt = sprintf('file %d. Previous file: %s', i, files.name{end});
    else
        prompt = sprintf('file %d. Press cancel when done.', i);
    end
    
    [filenames, filepaths] = uigetfile('*sbx', prompt, 'MultiSelect', 'on');    
    
    if ~iscell(filenames) & filenames == 0 % Cancel is pressed probably: stop with selecting
        selecting = false;
    elseif iscell(filenames) % multiple files have been selected at once
        nselected = length(filenames);
        
        files = [files; table(repmat({filepaths}, [nselected, 1]), filenames', 'VariableNames', varNames)];
        i = i + nselected;
    else
        files = [files; table({filepaths}, {filenames}, 'VariableNames', varNames)];
        i = i + 1;
    end
end
nfiles = height(files);

files.name = erase(files.name, '.sbx');