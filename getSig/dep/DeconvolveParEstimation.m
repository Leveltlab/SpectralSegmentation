function [tauest, aest, sigmaest, spikes] = DeconvolveWithParEstimation(sig, freq)
% Apply autocalibration of deconvolution parameters based on signal of one
% ROI
%
% Input:
%   sig ([t x 1] double): signal of one ROI
%   freq (scalar double): frequency of recorded signal
% 
% Output: Important parameters of the deconvolution process
%       the values of parameters A, tau and sigma.
%
% 
% adapted by Huub Terra & Leander de Kraker
% 2024-10-24
% 


% Initialize struct
try
    pax = spk_autocalibration('par');
    
    % (set saturation parameter)
    pax.saturation = 0.1;
    pax.dt = 1/freq;
    % (set limits for A and tau)
    pax.amin = 0.03;
    pax.amax = 0.20;
    pax.taumin = 0.1;
    pax.taumax = 1.5;
    % (when running MLspike from spk_autocalibratio, do not display graph summary)
    pax.mlspikepar.dographsummary = false;
    
    % perform auto-calibration
    [tauest, aest, sigmaest] = spk_autocalibration(sig, pax);
catch
    tauest   = NaN;
    aest     = NaN;
    sigmaest = NaN;
end