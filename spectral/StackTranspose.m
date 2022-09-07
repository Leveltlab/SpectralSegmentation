function StackTranspose(varargin)
%% Stack transpose, spectral correlation and roi segmentation
% javaaddpath( [gitpath '\2Pimage\java\Jimutil\Jimutil.jar'])
%
% Optional input:
% 1: (string) pathname and filename of the sbx file (e.g. 'D:\2pdata\mouse3.sbx')
% 2: (string) pathname of where to save the transposed file
% 
% In case only one input is given the trans file will be saved where the
% sbx file is located.
% When input is not given a pop-up will ask for the sbx filename and the
% location of where to save
% Chris van der Togt, 2017, 
% Netherlands Institute for Neuroscience 

global info
%% Transpose Stack from SBX file
if exist('varargin', 'var') && nargin == 2
    strfp = varargin{1};
    [pn, filename] = fileparts(strfp);
    strSavepath = varargin{2};

elseif exist('varargin', 'var') && nargin == 1
    strfp = varargin{1};
    [pn, filename] = fileparts(strfp);
    strSavepath = pn;
    
else
    [filename , pn] = uigetfile('*.sbx');
    if ~isempty(filename) 
        fnsplit = strsplit(filename, '.');
        filename = fnsplit{1};
    end
    strSavepath = uigetdir(pn,'Where to save the output file?');
end

strfp = fullfile(pn, filename);
%load([strfp '.mat'])
sbxread(strfp, 1, 1);

%save frequency with data file
if isfield(info, 'Freq')
    Inf.Freq = info.Freq;
elseif isfield(info, 'Slices')
    Inf.Freq = info.resfreq/info.recordsPerBuffer/3;
else
    Inf.Freq = info.resfreq/info.recordsPerBuffer;
end

d = dir([strfp, '.sbx']);
max_idx = d.bytes/info.nsamples; 

%parameters for transposing
Inf.StrPath = strSavepath; %Path for temporary files
Inf.Save = fullfile(strSavepath, [filename '_Trans.dat']);
Inf.Source = [strfp '.sbx'];        %source input .sbx file
Inf.Dimensions = [info.Shape(1) info.Shape(2) info.Shape(3) max_idx];

Zorder(Inf); %mex function to transpose the sbx file


