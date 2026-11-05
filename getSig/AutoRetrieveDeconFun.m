function AutoRetrieveDeconFun(filepaths, filenames, toRun, doParamEstimation, sigToLoad, decToSave)
% Retrieve signal and deconvolve for as many files as you want. 
% Can also apply Zscore or SEAL (Simple Estimate of Activity by Leander)
% 
% All input is optional. If run without input files, the function will ask
% how you want to select files, giving two options for file selection
% Without the other input, they will default as explained next:
% Input:
%  filepaths (string or [n x 1] cell array with strings): if empty, asks
%                                                         for selection
%  filenames (string or [n x 1] cell array with strings)
%  toRun (string or [1-4 x 1] cell array with strings): Options=  
%      'retrieveSignals', 'SealSignals' 'ZscoreSignals', 'DeconvolveSignals'
%       capitalisation not necessary
%  doParamEstimation = false; % estimate parameters for deconvolution?
%  sigToLoad ([m x 1] cell array with strings): Which signals  
%           DeconvolveSignals should analyse 
%           (defaults to 'sigCorrected' if empty)
%  decToSave ([m x 1] cell array with strings): What should the output 
%           names be of the loaded signals in DeconvolveSignals 
%           default: 'deconCorrected';
% 
% Output: new variables in the requested SPSIGs
% 
% Example:
% filepaths = 'E:\data\mymouse\';
% filenames = 'recordingA.SPSIG';
% toRun = 'retrieveSignals';
% doParamEstimation = false;
% sigToLoad = 'sigCorrected';
% decToSave = 'deconCorrected';
%
% Leander de Kraker
% NOT FULLY FINISHED TO LIKING: FILE SELECTION POP UPS SUCK:
%  TODO: filepaths only input can be used to pre-select file selection
%        method

arguments (Input)
    filepaths = {};
    filenames = {};
    toRun = {'retrieveSignals', 'DeconvolveSignals'};
    doParamEstimation = false;
    sigToLoad = {'sigCorrected'};
    decToSave = {'deconCorrected'};
end

% convert to cell if needed
if ~iscell(filepaths) && (ischar(filepaths) || isstring(filepaths))
    filepaths = {filepaths};
end
if ~iscell(filenames) && (ischar(filenames) || isstring(filenames))
    filepaths = {filenames};
end
if ~iscell(sigToLoad) && (ischar(sigToLoad) || isstring(sigToLoad))
    filepaths = {sigToLoad};
end
if ~iscell(decToSave) && (ischar(decToSave) || isstring(decToSave))
    filepaths = {decToSave};
end
nfiles = length(filenames);
btns = {'via current folder, automatically collect files', 'select folder and automatically collect files', 'Collect files manually'};

%% Is file selection needed?
if isempty(filenames)
    selectType = questdlg('', 'how to select files', btns{1}, btns{2}, btns{2});
    
    %% SELECTION option 1: 
    switch selectType
        case btns{1} % Collect files automatically via current folder
            %% Collect all files
            removeDuplicates = false;
            depth = 1;
            searchFolder = cd;
            [filepaths, filenames] = CollectFiles('SPSIG.mat', searchFolder, depth, removeDuplicates);
            nfiles = length(filenames);
        
        %% SELECTION option 2: 
        case btns{2} % Collect files manually
            %% pop-up, Press 'cancel' when done selecting.
            
            filenames = {};
            filepaths = {};
            
            selecting = true; % true as long as files are being selected
            i = 0;
            while selecting
                i = i + 1;
                [filenames{i}, filepaths{i}] = uigetfile('*SPSIG.mat', sprintf('file %d. Press cancel when done',i));
                
                if filenames{i} == 0 % Cancel is pressed probably: stop with selecting
                    filenames(i) = [];
                    filepaths(i) = [];
                    selecting = false;
                end
            end
            filenames = filenames';
            filepaths = filepaths';
            nfiles = length(filenames);
            clearvars selecting
        otherwise
            warning('some sort of impossible button pressed, should not happen')
    end
end

%% Do retrieve signals and deconvolution

for i = 1:nfiles
    if any(strcmpi(toRun, 'retrieveSignals'))
        retrievesignals([filepaths{i},  filenames{i}]);
    end
    if any(strcmpi(toRun, 'SealSignals'))
        SealSignals([filepaths{i},  filenames{i}]);
    end
    if any(strcmpi(toRun, 'ZscoreSignals'))
        ZscoreSignals([filepaths{i},  filenames{i}])
    end
    if any(strcmpi(toRun, 'DeconvolveSignals'))
        DeconvolveSignals([filepaths{i},  filenames{i}], doParamEstimation, sigToLoad, decToSave);
    end
end
