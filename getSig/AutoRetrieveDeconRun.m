% This script sets the input variables for AutoRetrieveDeconFun as an
% example.
% AutoRetrieveDeconFun retrieve signal and deconvolve for as many files as you want. 
% Can also apply Zscore or SEAL (Simple Estimate of Activity by Leander)
% 
% All input is optional. If run without input files, the function will ask
% how you want to select files, giving two options for file selection
% Without the other input, they will default as explained next:
% AutoRetrieveDeconFun input:
%  filepaths (string or [n x 1] cell array with strings)
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
% Leander de Kraker
% 2026-9-4
% NOT FINISHED TO LIKING, look into AutoRetrieveDeconFun for details

filepaths = [];
filenames = [];
toRun = {'retrieveSignals', 'DeconvolveSignals'};
doParamEstimation = false;
sigToLoad = 'sigCorrected';
decToSave = 'deconCorrected';

AutoRetrieveDeconFun(filepaths, filenames, toRun, doParamEstimation, sigToLoad, decToSave)
 