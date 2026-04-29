% When a chronic file contains recordings, which share a mask with other
% recordings, which are not in the chronic file, you can add those
% recordings in this script
% 
% Leander de Kraker
% 2022-12-14
% 

clear; clc
% Load chronic file
[chronicName, chronicPath] = uigetfile('*chronic.mat','get Chronic file');
load([chronicPath, chronicName], 'transformed', 'linkMat', 'score',...
            'PPs', 'Masks', 'nLinksMask', '', 'thres',...
            'filenames','filepaths','filenamesShort','filedates', 'nfiles')


%% Select compatible files AUTOMATICALLY (from same folder as respective chronic rec)

files = cell(nfiles, 1);
nfilesSame = zeros(1, nfiles);
easyPass = false; % Check for identical ROI size (may be too stringent)
for i = 1:nfiles
    % Find all SPSIG files in the same folder
    files{i} = dir([filepaths{i} '*SPSIG.mat']);
    toDel = strcmp({files{i}.name}, filenames{i});
    files{i}(toDel) = [];

    % Check if it has the same mask/contours
    nfilesSame(i) = length(files{i});
    toDel = true(nfilesSame(i), 1);
    for j = 1:nfilesSame(i)
        clearvars PP
        load([filepaths{i}, files{i}(j).name], 'PP')
        if PP.Cnt == PPs(i).Cnt % Same amount of ROIs?
            if all(PPs(i).A == PP.A) || easyPass % Really same ROIs more torough.
                toDel(j) = false;
            else
                fprintf('file %s has same amount of ROIs as %s\n',...
                        file{i}(j).name, filenames{i})
                fprintf('However, the ROI sizes were different\n')

            end
        end
    end
    files{i}(toDel) = [];
    nfilesSame(i) = length(files{i});
    fprintf('%d files with same mask found for file %d %s\n', nfilesSame(i), i, filenames{i})
end
filepathsAdd = cell(nfiles, 1);
filenamesAdd = cell(nfiles, 1);
for i = 1:nfiles
    for j = 1:nfilesSame(i)
        filepathsAdd{i}{j} = files{i}(j).folder;
        if ~strcmp(filepathsAdd{i}{j}(end), '\')
            filepathsAdd{i}{j} = [filepathsAdd{i}{j}, '\'];
        end
        filenamesAdd{i}{j} = files{i}(j).name;
    end
    filepathsAdd{i} = filepathsAdd{i}';
    filenamesAdd{i} = filenamesAdd{i}';
end
    


%% Select files MANUALLY

easyPass = false;  % Check for identical ROI size (may be too stringent)
filepathsAdd = cell(nfiles, 1);
filenamesAdd = cell(nfiles, 1);
currentDir = cd;
for i = 1:nfiles
    selecting = true;
    j = 0;
    toDel = true(1);

    cd(filepaths{i})
    while selecting
        j = j + 1;
        toDel(j) = true;
        if j == 1
            str = sprintf('load additional file for file %d: %s. Press cancel when done', i, filenames{i});
        else
            str = sprintf('load additional file %d for file %d: %s. Press cancel when done', j, i, filenames{i});
        end
        [filenamesAdd{i}{j}, filepathsAdd{i}{j}] = uigetfile('*SPSIG.mat', str);
        
        if filenamesAdd{i}{j} == 0 % Cancel is pressed probably: stop with selecting
            selecting = false;
        else % Get the date
            % double check for same mask
            clearvars PP
            load([filepathsAdd{i}{j}, filenamesAdd{i}{j}], 'PP')
            if PP.Cnt == PPs(i).Cnt % Same amount of ROIs?
                if all(PPs(i).A == PP.A) || easyPass % Really same ROIs more torough.
                    toDel(j) = false; % accept this file            
                else
                    fprintf('file %s has same amount of ROIs as %s\n',...
                            filenamesAdd{i}{j}, filenames{i})
                    fprintf('However, the ROI sizes were different\n')
                end
            end
        end
    end
    
    filenamesAdd{i}(toDel) = [];
    filepathsAdd{i}(toDel) = [];
    
    nfilesSame(i) = length(filenamesAdd{i});
    filepathsAdd{i} = filepathsAdd{i}';
    filenamesAdd{i} = filenamesAdd{i}';
end
cd(currentDir)


%% Make new filenames list

% Ask for cool short filename?
askShortName = true;

nfilesNew = (sum(nfilesSame)+nfiles);

% Figure out sorting of the new and original files in the new file list
idx = 1:nfiles;
idx(2:end) = idx(2:end)+cumsum(nfilesSame(1:end-1));
idxa = 1:nfilesNew;
idxa(idx) = [];

filenamesNew = cell(nfilesNew, 1);
filepathsNew = cell(nfilesNew, 1);
filenamesShortNew = cell(nfilesNew, 1);
filedatesNew = NaT(nfilesNew, 1);
filenamesNew(idx) = filenames;
filepathsNew(idx) = filepaths;
filenamesShortNew(idx) = filenamesShort;
filedatesNew(idx) = filedates;

count = 1;
for i = 1:nfiles
    
    % Propose new name for file which had other files with same mask
    if nfilesSame(i)~=0    
        strTit = sprintf('short original name file %d/%d', i, nfiles);
        filenamesShortNew(idx(i)) = ProposeGoodName(filenames{i},...
                                                    filepaths{i},...
                                                    filenamesShort{i},...
                                                    askShortName, strTit);
    end
    
    % Set filenames for newly added files
    for j = 1:nfilesSame(i)
        filenamesNew{idxa(count)} = filenamesAdd{i}{j};
        filepathsNew{idxa(count)} = filepathsAdd{i}{j};
        filedatesNew(idxa(count)) = FindDate(filepathsAdd{i}{j}, filenamesAdd{i}{j});
         
        strTit = sprintf('short name file %d/%d', count, sum(nfilesSame));
        strProposal = ShortenFileNames(filenamesNew(idxa(count)), 1, filedatesNew(idxa(count)));

        filenamesShortNew(idxa(count)) = ProposeGoodName(filenamesNew{idxa(count)},...
                                                         filepathsNew{idxa(count)},...
                                                         strProposal{:},...
                                                         askShortName, strTit);
        count = count + 1;
    end
end

filenames = filenamesNew;
filepaths = filepathsNew; 
filedates = filedatesNew;
filenamesShort = filenamesShortNew;
nfiles = nfilesNew;

clearvars filenamesShortNew filenamesNew filepathsNew filedatesNew nfilesNew
clearvars files j i count strTit PP easyPass askShortName


%% Register images & Masks & PPs

% Maximum registration shift



%% Add in linkMat

toPlace = zeros(nfiles, 1);
toPlace(idx) = 1;
toPlace = cumsum(toPlace);

linkMat = linkMat(:,toPlace);


%% Add info about which vari


%% A function which gets called in this script.
% DO NOT EXECUTE THE CODE WHERE THIS GETS CALLED VIA F9, (line by line). 
% then MATLAB will not understand it can find the function in this script
% EXECUTE THE CODE WHERE THIS GETS CALLED via run section (ctrl+enter)
function str = ProposeGoodName(fname, fpath, str, askShortName, strTit)
    baseName = strsplit(fname, {'_normcorr', '_SPSIG', '.mat'});
    baseName = baseName{1};
    
    % use JSON to get an idea of the stimulus type  
    found = false;
    if askShortName
        jsonName = [fpath baseName, '_session.json'];
        if exist(jsonName, 'file')
            fid = fopen(jsonName);
            rawjson = fread(fid,inf);
            fclose(fid);
            contentjson =  jsondecode(char(rawjson'));
            if isfield(contentjson, 'stimulus')
                str = [str '_' contentjson.stimulus];
                found = true;
            end
        end
        if ~found
            str = baseName;
        end
        
        strQuest = sprintf('How to call %s', fname);
        str = inputdlg(strQuest, strTit, 1, {str});
    else
        str = baseName;
    end
end

