function ZscoreSignals(varargin)
% Calculate z score of sigCorrected
% 
% Leander de Kraker
% 2024-7-4
% 
%%

if exist('varargin', 'var') && nargin == 1
    strFile = varargin{1};
else
    [strFileName, strFilePath] = uigetfile('*SPSIG.mat', 'Select SPSIG file to z score');
    strFile = [strFilePath, strFileName];
end

load(strFile, 'sigCorrected');
sigCorrected_Z = zscore(sigCorrected);
save(strFile, 'sigCorrected_Z', '-append')
fprintf('saved sigCorrected_Z\n')