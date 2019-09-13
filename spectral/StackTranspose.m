function StackTranspose(varargin)
%% Stack transpose, spectral correlation and roi segmentation
% javaaddpath( [gitpath '\2Pimage\java\Jimutil\Jimutil.jar'])
%
% Optional input: string with folder name and trans file name, (e.g. blafolder\mouse_Trans.dat)
% When input is not given a pop-up will ask for that filename
%

global info
%% Transpose Stack from SBX file
if exist('varargin', 'var') && nargin == 2
    strfn = varargin{1};
    [fp,fn] = fileparts(strfn);
    load(fullfile(fp, fn))
     
    strSavepath = varargin{2};

else
    [fn , pn] = uigetfile('*.sbx');
    if ~isempty(fn)
        strfn = fullfile(pn, fn);
        fnsplit = strsplit(fn, '.');
        fn = fnsplit{1};
        strfp = [pn fn];
        load(strfp)
    end

    strSavepath = uigetdir(pn,'Where to save the output file?');
end
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
Inf.Save = [strSavepath '\' fn '_Trans.dat'];
Inf.Source = strfn;        %source input .sbx file
Inf.Dimensions = [info.Shape(1) info.Shape(2) info.Shape(3) max_idx];

Zorder(Inf); %mex function to transpose the sbx file


clear all
