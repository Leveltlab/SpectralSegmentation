function [hz, scaleUm, FOVum, pixelAspectRatio, squareFOV, answers] = RequestRecInfo(varargin)
% Get some important (hz) and less important (scale and aspect ratio)
% information of a recording
% 
% Input (optional): 
%   - defInput ([6x1] cell array): default input to the 5 questions. and
%                                  the filename
%   - dimensions of image ([2x1] double): To propose aspect ratio
% 
% Output:
%   - hz (double): Imaging frequency
%   - scaleUm (double or empty): the width of one pixel in micrometer (um)
%   - FOVum   (double or empty): the width of the entire FOV in micrometer
%   - pixelAspectRatio (double or empty): How much higher a pixel is than
%           it is wide.
%           example: a square image recorded with dimensions 16x9 pixels 
%               (widthxheight) results in pixel aspect ratio 16/9 = 1.7778
%   - squareFOV (boolean): Was the imaged field square? True | false
%           if true, then the pixel aspect ratio is the width/height of the
%           image. And can be used if the pixel aspect ratio is not known
%   - answers ([5x1 cell array]) with the answers of the dialog.
% 
% See also: RequestRecInfoProcess
% 
% Leander de Kraker
% 2024-10-30
% 

if exist('varargin', 'var') && nargin>=1 && ~isempty(varargin{1})
    defInput = varargin{1};
    if ~isa(defInput{5}, 'char') & defInput{5} % if squareFOV == true
        defInput{5} = 'yes';
    elseif ~isa(defInput{5}, 'char')
        defInput{5} = 'no';
    end
else
    defInput = cell(5, 1);
    defInput{5} = 'yes';
end

for i = 1:4
    if ~ischar(defInput{i})
        defInput{i} = num2str(defInput{i});
    end
end

if exist('varargin', 'var') && isempty(defInput{4}) && nargin==2 && strcmp(defInput{5}, 'yes')
    defInput{4} = num2str(varargin{2}(2) / varargin{2}(1));
end

questions = {'Framerate for this dataset? (Hz)';...
             '[optional]   Horizontal scale of one pixel? (um). Or fill in the next one to calculate this';...
             '[optional]   Width of imaged field of view  (um)';...
             '[optional]   Aspect ratio of one pixel (how much higher than wide) (can be calculated if field is square = yes)';...
             '[yes / no]   Imaged field is square in reality?'};
dlgTitle = 'Give info of recording ';
if length(defInput)==6 && ischar(defInput{6})
    recTitle = strsplit(defInput{6}, '.mat');
    dlgTitle = [dlgTitle, recTitle{1}];
end
inDim = [ones(5,1), 120*ones(5,1)];

answers = inputdlg(questions, dlgTitle, inDim, defInput);


hz = str2double(answers{1});
scaleUm = str2double(answers{2});
if isnan(scaleUm)
    scaleUm = [];
end
FOVum = str2double(answers{3});
if isnan(FOVum)
    FOVum = [];
end
pixelAspectRatio = str2double(answers{4});
if isnan(pixelAspectRatio)
    pixelAspectRatio = [];
end
squareFOV = answers{5};
if strcmpi(squareFOV, 'yes')
    squareFOV = true;
else
    squareFOV = false;
end
