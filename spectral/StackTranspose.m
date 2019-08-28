%% Stack transpose, spectral correlation and roi segmentation
% javaaddpath( [gitpath '\2Pimage\java\Jimutil\Jimutil.jar'])

%% Transpose Stack from SBX file
global info
[fn , pn] = uigetfile('*.mat');
if ~isempty(fn)
    fnsplit = strsplit(fn, '.');
    filename = fnsplit{1};
    strfp = [pn filename];
    load(strfp)
    strfn = [strfp '.sbx'];
end
strSavepath = uigetdir(pn,'Where to save the output file?');
%save frequency with data file
%freq = 30.1; %image sampling frequency
if isfield(info, 'Freq')
    Inf.Freq = info.Freq;
    
elseif isfield(info, 'Slices')
    Inf.Freq = info.resfreq/info.recordsPerBuffer/3;
    
else
    Inf.Freq = info.resfreq/info.recordsPerBuffer;
end


d = dir(strfn);
max_idx = d.bytes/info.nsamples;

%parameters for transposing
Inf.StrPath = strSavepath; %Path for temporary files
Inf.Save = [strSavepath '\' filename '_Trans.dat'];
Inf.Source = strfn;        %source input .sbx file
Inf.Dimensions = [info.Shape(1) info.Shape(2) info.Shape(3) max_idx];

Zorder(Inf); %mex function to transpose the sbx file


clear all
