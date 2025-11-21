function [filepaths, filenames, fileList] = CollectFiles(searcher, strFolder, depth, removeDuplicates)
% [filepaths, filenames, fileList] = CollectFiles('SPSIG_Res.mat', 'E:\mydata', 1, true)
% 
% [filepaths, filenames] = CollectFiles; Returns all Res files of a mouse, 
% by looking in the folders inside the current folder (so a search of 
% 'depths' deep). default depth=1
% 
% Input (all optional):  
% - searcher (string): string to use to search for files. 
%                       (by default: 'SPSIG_Res.mat')
% - strFolder (string): which folder to start search from
% - depth (integer): how deep, folder-wise, to search for files. 0=search
% in current folder. 1 (default), look one folder deep, 2, look inside the
%   folders that are inside those folders, etc.
% - removeDuplicates (true | false): true deletes the different splits of a 
%       recording, leaving only unique recordings. (by default: false)
% 
% Output:
% - filepaths ([n x 1] cell with strings)
% - filenames ([n x 1] cell with strings)
% - fileList  ([n x 1] struct with all info on the files)
% 
% Leander de Kraker
% 2025-11-21
%

arguments (Input)
    searcher = 'SPSIG_Res.mat'
    strFolder = cd
    depth (1, 1) = 1
    removeDuplicates (1, 1) logical = false
end

% remove \ if present in folderStr
if ~strcmp(strFolder(end), '\')
    strFolder = [strFolder '\'];
end

searcher = strcat(repmat('*', [1, min(depth, 2)]),...
                  repmat('\', [1, depth>0]),...
                  '*',...
                  searcher);
fileList = dir([strFolder, searcher]);
filenames = {fileList.name}';
filepaths = {fileList.folder}';
filepaths = strcat(filepaths, '\');

% Remove duplicate files of same recording, recorded via different splits
if removeDuplicates
    duplicateKeyStr = 'split'; % If this string is present in the filename, 
                               % allow removal from list.
    uniqueFileNamePart = 3; % This part, denoted by '_' of the filename 
                            % should be unique. Use it to determine if the
                            % file has been seen or not
    [foldersCount, foldersUnique] = groupcounts(filepaths);
    foldersToCheck = foldersUnique(foldersCount>1);
    nfiles = length(filenames);
    toDel = false(nfiles, 1);
    for i = 1:length(foldersToCheck)
        idx = find(strcmp(filepaths, foldersToCheck{i}));
        n = length(idx);
        % the 3rd part of the filename should be an identifier for the
        % recording (within each folder)
        recNr = ShortenFileNames(filenames(idx), uniqueFileNamePart);
        recNrPresent = unique(recNr);
        recNrSeen = false(n, 1);
        % Also check if the word 'split' is present in the file.
        splitStrPresent = contains(lower(filenames(idx)), duplicateKeyStr);
        for j = 1:n
            seenIdx = find(strcmp(recNrPresent, recNr{j}));
            if recNrSeen(seenIdx) && splitStrPresent(j)
                % if the recNr has been seen before and 'split' is in the
                % name, prepare to delete file.
                toDel(idx(j)) = true;
            else
                % We now saw that recNr, tell recNrSeen
                recNrSeen(seenIdx) = true;
            end
        end
    end
    
    if sum(toDel)>0
        fprintf('Deleting %d of %d files because they seem duplicates\n', sum(toDel), nfiles)
        fileList(toDel) = [];
        filenames(toDel) = [];
        filepaths(toDel) = [];
    else
        fprintf('Only unique files seem present :)\n')
    end
end



