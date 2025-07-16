function infoSigPar = DeconvolveSignals(filename, varargin)
% Deconvolves requested, or sigCorrected (default) variables which are in 
% the file and saves the results to that file.
% 
% filename = 'C:\data\mymousy\Mikey_SPSIG.mat';
% estimatePars = true; % let MLspike search parameters
% signalNames = {'sigCorrected', 'sig'};
% deconNames = {'deconCorrected', decon'};
% infoSigPar = DeconvolveSignals(filename, estimatePars, signalNames, deconNames);
% infoSigPar = DeconvolveSignals;
% 
% 
% Input: 
%  - filename (string): SPSIG filename, with sig variable present
%                       can include path
% Optional inputs (any can by left empty):
%  - estimatePars (boolean): Use parameters from file (false) |default|
%                            or find parameters (true)
%  - signalNames (string or cell array with strings): signals to deconvolve
%                                       default = 'sigCorrected'
%  - deconNames (string or cell array with strings): names of the output
%            deconvolved signals, corresponding to the signalNames
%            default = 'deconCorrected', or the input sigNames with decon 
%               in front of name and the original first letter capitalized.
% 
% Output: Updated SPSIG file.
%   - infoSigPar (struct): info on Signal retrieval Parameters
%
%
% Can also be executed as a script/ without input, in that case a filename
% and whether to use parameter estimation will be requested.
% 
% Deconvolving using MLspike
% 
% 
% Chris v.d. Togt & Leander de Kraker
% 2024-10-24: Enabled processing of specific requested signals
%             Enabled estimation of deconvolution parameters.
% 

%% Select file if code is being called as script

if exist('varargin', 'var') && nargin>=1 && ~isempty(filename)
    asFunction = true;
else
    asFunction = false;
    [fn, pn] = uigetfile('*_SPSIG.mat');
    if fn == 0 % if cancelled
        fprintf('no file selected: exiting\n')
        return
    end
    filename = [pn fn];
end
varNames = who('-file', filename);

if asFunction && nargin==2 && ~isempty(varargin{1}) % to do estimation of parameters is given
    doParEstimation = varargin{1};
elseif asFunction % Not given as input but executed as function. Use standard option
    doParEstimation = false;
else
    doParEstimation = questdlg('Estimate deconvolution parameters?', 'estimate MLspike parameters', ...
                            'yes', 'Load from premade file', 'Load from premade file');
    if strcmp(doParEstimation, 'yes')
        doParEstimation = true;
    else 
        doParEstimation = false;
    end
end

% Sigs to deconvolve
if asFunction && nargin==3
    sigNames = varargin{2};
    sigNamesGiven = true;
    if iscell(sigNames) % multiple sig names were given
        nsigs = length(sigNames);
    else
        nsigs = 1;
        sigNames = {sigNames};
    end
else % Use default signal name to load
    sigNamesGiven = false;

    % find valid signal
    if ismember('sigCorrected', varNames)
        sigNames = {'sigCorrected'};
    elseif ismember('sig', varNames)
        sigNames = {'sig'};
    else
        error('No signal found to deconvolve! Please specify in input')
    end
    nsigs = 1;
end

% Names to give to deconvolved signals
if asFunction && nargin==4 % Output names are given in input
    deconNames = varargin{3};
    if ~iscell(deconNames)
        deconNames = {deconNames};
    end
else
    if sigNamesGiven % No output names defined, but input names were given
        for i = 1:nsigs
            deconNames = ['decon', upper(sigNames{i}(1)), sigNames{i}(2:end)];
        end
    else
        if strcmp(sigNames{1}, 'sigCorrected')
            deconNames = {'deconCorrected'};
        elseif strcmp(sigNames{1}, 'sig')
            deconNames = {'decon'};
        else
            warning('Do not know how to name deconvolved signal! Please specify')
            warning('help DeonvolveSignals')
        end
    end
end


% Load the data
data = load(filename, sigNames{:}, 'freq', 'infoSigPar');
if isfield(data, 'infoSigPar')
    infoSigPar = data.infoSigPar;
else % Is the really old name infoSigRetrieval present?
    data2 = load(filename, 'infoSigRetrieval');
    if isfield(data2, 'infoSigRetrieval')
        infoSigPar = data2.infoSigRetrieval;
    else
        infoSigPar = struct([]);
    end
end
% Struct with info of how the retrieved signal is done
infoSigPar(1).deconv = 'MLspike';
infoSigPar.deconParEstimated = doParEstimation;

freq = data.freq;
nrois = size(data.(sigNames{1}), 2);
nt = size(data.(sigNames{1}), 1);
t = (1:nt+1)/freq;

if ~doParEstimation % parameters from file, please check params with G = spk_demoGUI
    parFile = 'MLpest20180927.mat';
    % parFile = 'MLPest20210313.mat';
    load(parFile, 'par'); % load par
else % pre-allocate variables for parameter estimation
    par = tps_mlspikes('par'); 
    par.saturation = 0.1;
    par.drift.parameter = .01;
    aest = nan(nrois, 1);
    tauest = nan(nrois, 1);
    sigmaest = nan(nrois, 1);
end
par.dt = 1/freq;

% MLspike deconvolution
tic
decon = struct();
for s = 1:nsigs
    sig = data.(sigNames{s});

    % signal should be fluctuating with baseline on 0
    if median(sig(:))<0.5
        sig = sig + 1;
        fprintf('Summed one to the signal to get baseline to 1\n')
    end
    
    deconS  = zeros(size(sig));
    
    % Define waitbar
    if nsigs>1
        strSigCounter =  sprintf(' (%d/%d)', s, nsigs);
    else
        strSigCounter = '';
    end
    hWb = waitbar(0,...
                  sprintf('Running MLspike on %s %s (Runtime %1.2f min)',sigNames{s}, strSigCounter, 0.00),...
                  'Name','Deconvolution');
    global h N p
    h = hWb;
    N = nrois;
    p = 1;
    D = parallel.pool.DataQueue;
    afterEach(D, @retrievesignals_parWB);
    parfor roi = 1:N % parallel processing should be possible, only the waitbar won't make much sense.
        if doParEstimation
            [tauest(roi), aest(roi), sigmaest(roi), spikes] = DeconvolveWithParEstimation(sig(:,roi), freq, par);
        else
            [spikes, ~, ~] = spk_est(sig(:,roi), par);
        end

        if ~isempty(spikes)
            out = arrayfun(@(x) x>=t(1:end-1)&x<t(2:end), spikes', 'UniformOutput', false);
            deconS(:,roi) = sum(cell2mat(out))';
        end
        send(D, roi)
    end
    close(hWb)
    fprintf('MLspike took %1.2f minutes to process %d rois. signal %s%s.\n',...
            toc/60, nrois, sigNames{s}, strSigCounter)
    
    % Check deconvolution for dissapointments
    dissapointment0 = find(all(deconS==0)); % no spikes at all for these ROIs
    dissapointment1 = find(all(deconS==1)); % one spike at every timepoint for these ROIs...
    if ~isempty(dissapointment0)
        fprintf('%d ROIs do not have any spikes in %s. ROIs:', length(dissapointment0), deconNames{s})
        for i = 1:length(dissapointment0); fprintf(' %d', dissapointment0(i)); end; fprintf('\n')
    end
    if ~isempty(dissapointment1)
        fprintf('%d ROIs have constant spikes in %s. ROIs:', length(dissapointment1), deconNames{s})
        for i = 1:length(dissapointment1); fprintf(' %d', dissapointment1(i)); end; fprintf('\n')
        % Set all the ROIs that only have 1 to 0
        deconS(:,dissapointment1) = 0;
        fprintf('Their values have been set to 0\n')
    end
    infoSigPar.(deconNames{s}).dissapointment0 = dissapointment0;
    infoSigPar.(deconNames{s}).dissapointment1 = dissapointment1;
    
    % Save the estimated parameters in a more managable format
    if doParEstimation 
        par.a = aest;
        par.tau = tauest;
        par.finetune.sigma = sigmaest;
    end
    infoSigPar.(deconNames{s}).par = par;
    
    decon.(deconNames{s}) = deconS;
end

decon.infoSigPar = infoSigPar;

% Add fields of struct as variables
save(filename, '-struct', 'decon', '-append')



