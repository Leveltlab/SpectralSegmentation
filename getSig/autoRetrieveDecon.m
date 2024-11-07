% Retrieve signal and deconvolve for as many files as you want. 
%  
% Press 'cancel' when done selecting.
% 
% 

filenames = {};
filepaths = {};

selecting = true; % true as long as files are being selected
i = 0;
while selecting
    i = i + 1;
    [filenames{i}, filepaths{i}] = uigetfile('*SPSIG.mat', sprintf('file %d',i));
    
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

%% Do retrieve signals and deconvolution

doParamEstimation = false; % estimate parameters for deconvolution?
sigToLoad = {'sigCorrected'};
decToSave ={'deconCorrected'};

for i = 1:nfiles
    % retrievesignals([filepaths{i},  filenames{i}]);
    % SealSignals([filepaths{i},  filenames{i}]);
    DeconvolveSignals([filepaths{i},  filenames{i}], doParamEstimation, sigToLoad, decToSave);
end
