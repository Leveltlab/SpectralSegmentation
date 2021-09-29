function infoSigPar = DeconvolveSignals(filename)
%
% infoSigRetreival = DeconvolveSignal(pn, fn) deconvolves the sig and
% sigCorrected variables which are present in the file [pn fn]
% 
% Input: fn = (string) spectral filename, with sig variable present
%        pn = (string) pathname
% In case no input is given the script will ask for the filename
%
% Output: infoSigPar = info on Signal retrieval Parameters (struct)
%
% Deconvolving using MLspike
%

%% Select file if code is being called as script
try
    asfunction = nargin == 1;
catch
    asfunction = false;
end
if ~asfunction
    [fn, pn] = uigetfile('*_SPSIG.mat');
    if fn == 0 % if cancelled
        fprintf('no file selected: exiting\n')
        return
    end
    filename = [pn fn];
end

% Activate much tighter subplots
% [subplot margin top&side],[figure bottomspace,topspace],[leftspace,rightspace]
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.06], [0.1 0.1], [0.1 0.02]);

% Load the data
load(filename, 'sig', 'sigCorrected', 'infoSigPar', 'freq');
doSigCorrected = exist('sigCorrected','var');

% parameters, please check params with G = spk_demoGUI
parallel_process = true; % default
%MLspike_params = load('MLpest20180927', 'MLpest');
MLspike_params = load('MLPest20210313.mat', 'par'); %these are better parameters, but please check with G = spk_demoGUI
MLspike_params = MLspike_params.par;
MLspike_params.dt = 1/freq;

nrois = size(sig,2);
t = (1:length(sig)+1)/freq;

drifting_baseline = true; % default
if ~drifting_baseline
    % Prevent drifting baseline estimation during MLspike deconvolution
    MLspike_params = rmfield(MLspike_params,'drift');
    dbsl_process = 'fixed baseline';
else
    dbsl_process = 'drifting baseline';
end

MLpar = repmat(MLspike_params, nrois, 1);

% Data matrix preallocation
decon  = zeros(size(sig));
den    = zeros(size(sig));
spikes = cell(nrois,1);
drift  = zeros(size(sig));
if doSigCorrected
    deconCorrected  = zeros(size(sig));
    denCorrected    = zeros(size(sig));
    spikesCorrected = cell(nrois,1);
end

% MLspike deconvolution
tic
if parallel_process
    if ~verLessThan('matlab','9.3')
        % Define waitbar
        hWb = waitbar(0,sprintf('Running MLspike (Runtime %1.2f min)',0.00),...
                                'Name','Deconvolution');
        global h N p
        D = parallel.pool.DataQueue;
        afterEach(D, @retrievesignals_parWB);
        h = hWb;
        N = nrois;
        p = 1;
        parfor roi = 1:N % parallel processing should be possible, only the waitbar won't make much sense.
            [spikes{roi}, den(:,roi), drift(:,roi)] = spk_est(sig(:,roi), MLpar(roi));
            send(D,roi)
        end
        close(hWb)
        
        if doSigCorrected
            % Start deconvolution on neuropil corrected data
            hWb = waitbar(0,sprintf('Running MLspike on neuropil corrected data (Runtime %1.2f min)',0.00),...
                                    'Name','Deconvolution');
            h = hWb;
            parfor roi = 1:N % parallel processing should be possible, only the waitbar won't make much sense.
                [spikesCorrected{roi}, denCorrected(:,roi)] = spk_est(sigCorrected(:,roi), MLpar(roi));
                send(D,roi)
            end
            close(hWb)
        else
            fprintf('no sigCorrected present. skipping sigCorrected\n')
        end
        
    else % older versions of matlab (<2017b)
        
        N = nrois;
        parfor roi = 1:N % parallel processing should be possible, only the waitbar won't make much sense.           
            [spikes{roi}, den(:,roi), drift(:,roi)] = spk_est(sig(:,roi),MLpar(roi));
        end
        if doSigCorrected
            parfor roi = 1:N % parallel processing should be possible, only the waitbar won't make much sense.           
                [spikesCorrected{roi}, denCorrected(:,roi)] = spk_est(sigCorrected(:,roi),MLpar(roi));
            end
        else
            fprintf('no sigCorrected present. skipping sigCorrected\n')
        end
    end
else % Not parallel process
    
    fprintf('deconvolution at ROI: ')
    for roi = 1:nrois % parallel processing should be possible, only the waitbar won't make much sense.
        [spikes{roi},den(:,roi), drift(:,roi)] = spk_est(sig(:,roi), MLpar(roi));
        fprintf('%d ', roi)
    end
    fprintf('\n')
    if doSigCorrected
        fprintf('neuropil corrected deconvolution at ROI: ')
        for roi = 1:nrois % parallel processing should be possible, only the waitbar won't make much sense.
            [spikesCorrected{roi}, denCorrected(:,roi)] = spk_est(sigCorrected(:,roi), MLpar(roi));
            fprintf('%d ', roi)
        end    
        fprintf('\n')
    else
        fprintf('no sigCorrected present. skipping sigCorrected\n')
    end
end

fprintf('MLspike took %1.2f minutes to process %d rois (%s).\n',toc/60,...
        nrois,dbsl_process)


for roi = 1:nrois % parallel processing should be possible, yields little.
    % decon(:,roi) = sum(spikes{roi}'>=t(1:end-1)&spikes{roi}'<t(2:end))';
    if ~isempty(spikes{roi})
        out = arrayfun(@(x) x>=t(1:end-1)&x<t(2:end), spikes{roi}', 'UniformOutput', false);
        decon(:,roi) = sum(cell2mat(out))';
    end
    if doSigCorrected
        if ~isempty(spikesCorrected{roi})
            out = arrayfun(@(x) x>=t(1:end-1)&x<t(2:end), spikesCorrected{roi}', 'UniformOutput', false);
            deconCorrected(:,roi) = sum(cell2mat(out))';
        end
    end
end

% Check deconvolution for dissapointments
dissapointment0 = find(all(decon==0));
dissapointment1 = find(all(decon==1));

if ~isempty(dissapointment0)
    fprintf('%d ROIs do not have any spikes in decon. ROIs:', length(dissapointment0))
    for i = 1:length(dissapointment0); fprintf(' %d', dissapointment0(i)); end; fprintf('\n')
    infoSigPar.dissapointment0 = dissapointment0;
end
if ~isempty(dissapointment1)
    fprintf('%d ROIs have constant spikes in decon. ROIs:', length(dissapointment1))
    for i = 1:length(dissapointment1); fprintf(' %d', dissapointment1(i)); end; fprintf('\n')
    infoSigPar.dissapointment1 = dissapointment1;
    % Set all the ROIs that only have 1 to 0
    decon(:,dissapointment1) = 0;
    fprintf('Their values have been set to 0\n')
end

if doSigCorrected
    dissapointmentCorrected0 = find(all(deconCorrected==0));
    dissapointmentCorrected1 = find(all(deconCorrected==1));
    if ~isempty(dissapointmentCorrected0)
        fprintf('%d ROIs do not have any spikes in deconCorrected. ROIs:', length(dissapointmentCorrected0))
        for i = 1:length(dissapointmentCorrected0)
            fprintf(' %d', dissapointmentCorrected0(i))
        end; fprintf('\n')
        infoSigPar.dissapointmentCorrected0 = dissapointmentCorrected0;
    end
    if ~isempty(dissapointmentCorrected1)
        fprintf('%d ROIs have constant spikes in deconCorrected. ROIs:', length(dissapointmentCorrected1))
        for i = 1:length(dissapointmentCorrected1)
            fprintf(' %d', dissapointmentCorrected1(i))
        end; fprintf('\n')
        infoSigPar.dissapointmentCorrected1 = dissapointment1;
        % Set all the ROIs that only have value 1 to 0
        deconCorrected(:,dissapointmentCorrected1) = 0;
        fprintf('Their values have been set to 0\n')
    end
end

% Struct with info of how the retrieved signal is done
infoSigPar.deconv = 'MLspike';
infoSigPar.decon  = MLspike_params;

if doSigCorrected
    save(filename, 'den', 'decon', 'infoSigPar', 'denCorrected', 'deconCorrected', 'drift', '-append')
else
    save(filename, 'den', 'decon', 'infoSigPar', 'drift', '-append')
end

% Plot deconvoluted data, yellow for non-neuropil corrected, blue for neuropil corrected 
figure('units', 'normalized', 'position',  [0.75 0.2 0.24 0.5], 'Name', 'deconvolved')
subplot(1,1,1)
xas = (1:length(sig))./freq;
RGB = cat(3,decon',decon',decon');
if doSigCorrected % neuropil corrected decon signal in blue color channel
    RGB(:,:,3) = deconCorrected';
end
imagesc(xas./60, 1:size(decon,2), RGB);
title('spike estimates'); xlabel('time (minutes)')
xlim([0 xas(end)/60]); ylabel('ROI')

% % Plot all spikes using dot markers (slow)
% set(gca,'Color','k')
% hold on
% for i = 1:length(spikes)
%     if ~isempty(spikes{i})
%         plot(spikes{i},i,'.','color',[1 0.5 0],'markersize',3)
%     end
% end
% if doSigCorrected % neuropil corrected decon signal in blue color channel
%     for i = 1:length(spikesCorrected)
%         if ~isempty(spikesCorrected{i})
%             plot(spikesCorrected{i},i,'.w','markersize',3)
%         end
%     end
% end
% set(gca,'YDir','reverse')

