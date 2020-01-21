% This script runs SpecProfileCalcFun.m to calculate the PP.SpecProfile 
% and PP.peakFreq variable
% 
% PP.SpecProfile: 2D double [ROIs x frequency]: average spectral power
%               values of each ROI for each spectral frequency/ component
% PP.peakFreq: 1D vecor [ROIs]: which frequency has the highest power 
% 
%
% Leander de Kraker
% 2020-1-21
% 

% Load the mask, contour and spectral data
[fnSpsig, pnSpsig] = uigetfile('*SPSIG.mat','load mask data (SPSIG file)');
nameSpsig = [pnSpsig, fnSpsig];
load(nameSpsig, 'PP', 'Mask', 'SPic', 'Sax')

spec = log(SPic); % spectral data
vals = sort(unique(spec(:))); % remove -inf values
if vals(1)==-inf
    spec(spec==-inf) = vals(2);
end

% Check if SpecProfile is already present etc
if isfield(PP, 'SpecProfile')
    fprintf('Specprofile is already present in the datasets PP variable\n')
    if size(PP.SpecProfile,2) ~= PP.Cnt
        fprintf('But there seems to be an incorrect number of ROIs in PP.specProfile\n')
    end
end


% Delete the 0Hz component if it's present
if Sax(1) == 0
    if size(spec, 3) == length(Sax)
        spec = spec(:, :, 2:end);
    end
    Sax = Sax(2:end);
end

% Get the Specprofile & peak frequency
tic
[PP.SpecProfile, PP.peakFreq] = SpecProfileCalcFun(spec, Mask, Sax);
fprintf('done (%.1f seconds)\n',  toc)


%% Save the updated PP variable

fprintf('saving the updated PP to the SPSIG file\n')
save(nameSpsig, 'PP', '-append')
fprintf('saved %s\n', fnSpsig)

