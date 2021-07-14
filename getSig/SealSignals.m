function infoSigPar = SealSignals(filename)
% Simple Estimate of Activity by Leander (seal)
%
% infoSigRetreival = SealSignal(filename) does a simple and fast
%   calculation on sig and sigCorrected in the file 'filename'
% 
% Input: filename = (string) spectral filename, with sig variable present
% In case no input is given the script will ask for the filename
%
% Output: infoSigPar = info on Signal retrieval Parameters (struct)
%
% SEAL smooths the signal a bit with a gaussian filter of windowSize
% SEAL takes the derivative of the smoothed signal, 
% and removes negative values.
% This is very super fast and will result in high signal when the
% fluorescence signal is rising. 
% It will underestimate activity in case of 'sustained high
% fluorescence because of longer activity'. 
% But the signal is not sparse like the spike estimate

% Select file if code is being called as script
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

% Load the data
load(filename, 'sig', 'sigBack', 'infoSigPar', 'freq');
if ~exist('infoSigPar', 'var')
    infoSigPar.alpha = 0.7;
end

% Activate much tighter subplots
% [subplot margin top&side],[figure bottomspace,topspace],[leftspace,rightspace]
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.06], [0.1 0.1], [0.1 0.02]);

noiseFloor = 0;
% filter parameters to reduce noise
windowSize = 1;
% b = (1/windowSize)*ones(1,windowSize);

seal = imgaussfilt(sig, [windowSize, 0.01]);
seal = diff(seal) - noiseFloor;
seal(end+1, :) = 0; % fill at end
seal(seal<0) = 0;
sealBack = imgaussfilt(sigBack, [windowSize, 0.01]);
sealBack = diff(sealBack) - noiseFloor;
sealBack(end+1, :) = 0;
sealBack(sealBack<0) = 0;
sealCorrected = seal - infoSigPar.alpha * sealBack;
sealCorrected(sealCorrected<0) = 0;

% plot seal
figure('units', 'normalized', 'position',  [0.5 0.2 0.24 0.5], 'Name', 'SEAL')
subplot(1,1,1)
RGB = CreateRGB({seal', sealBack'}, 'rg b');
xas = (1:length(sig))./freq;
imagesc(xas./60, 1:size(seal,2), RGB);
title(sprintf('seal (Simple Estimate of Activity by Leander)\n%s',filename))
xlabel('time (minutes)') 
ylabel('ROI (n)')
xlim([0 xas(end)/60]);


%% save seal
infoSigPar.noiseFloorSEAL = noiseFloor;
infoSigPar.filterSizeSEAL = windowSize;

save(filename, 'seal', 'sealBack', 'sealCorrected', 'infoSigPar', '-append')
fprintf('saved the seal\n')


