function [rgb, colors, colorbarVals] = SpectralColorImg(loadType, inData, varargin)
% [rgb, colors, colorbarVals] = SpectralColorImg('file', filename, [0 0.5], true)
% [rgb, colors, colorbarVals] = SpectralColorImg('data', {SPic, Sax, spar}, [], false)
% [rgb, colors, colorbarVals] = SpectralColorImg('data', {SPic, Sax, spar})
% 
% Create a color coded image of the spectral images
% 
% Input:
% 1 - inputType (string). load data or get data in input? 'file', or 'data'
% 2 - inData. In case of 'file': 
%     filename (string). SPSIG filename, which should contain SPic and Sax.
% 2  -inData. In case of 'data': 
%     ([3 x 1] cell array): 1st SPic, 2nd Sax, 3rd (optional) spar.
%   (optional inputs 3 and 4)
% 3  - freqToUse ([1 x 2] double). Frequencies in Hz to use.
% 4  - plotter (scalar logical). Plot the output? (in an existing axes)
% 
% Output:
%   - rgb ([h x w x 3] double): The colored spectral image.
%   - colors ([n x 3] double): The colors that correspond to each selected
%                              frequency.
%   - colorbarVals ([n x 1] double): The frequencies that are used.
%
%
% Leander de Kraker
% 2022-8-25
% 


norma = true; % normalize each frequency?
wbalance = false; % white balance red green and blue color channels?

% Load data or use inputted data
if strcmp(loadType, 'file')
    load(inData, 'Sax', 'SPic', 'spar')
    if ~exist('SPic', 'var')
        rgb = 'no SPic found in file';
        colors = 'no SPic found in file';
        colorbarVals = 'no SPic found in file';
        warning('SPic & Sax variables are essential for SpectralColorImg')
        return
    end
else
    SPic = inData{1};
    Sax = inData{2};
    if length(inData)==3
        spar = inData{3};
    end
end

% Get frequencies to select
if exist('varargin', 'var') && nargin >= 3 && ~isempty(varargin{1})
    freqToUse = varargin{1};
    if length(freqToUse) == 1 % In case of bad input
        freqToUse = [0, freqToUse];
    end
elseif exist('spar','var')
    if isfield(spar, 'cutOffHzMin')
        freqToUse = [spar.cutOffHzMin, spar.cutOffHzMax]; % Frequency needs to be within this range
        fprintf('Using frequencies selected by ROI selection: [%.2f %.2f Hz]\n',...
                                        spar.cutOffHzMin, spar.cutOffHzMax)
    else
        fprintf('spar seems outdated/ incomplete. Using entire specAxis\n')
        freqToUse = [-1 1000];
    end
else
    freqToUse = [-1 1000];
    fprintf('Using entire specAxis\n')
end

if exist('varargin', 'var') && nargin == 4
    plotter = varargin{2};
else
    plotter = false;
end

idxToUse = Sax>freqToUse(1) & Sax<freqToUse(2);
% If somehow no frequencies were selected, select some anyway. We want an image!
if sum(freqToUse)==0
    freqToUse(1:end) = true;
end

imgStack = SPic(:,:,idxToUse);
colorbarVals = Sax(idxToUse)';

vals = sort(unique(imgStack(:))); % remove -inf values
if vals(1)==-inf
    imgStack(imgStack==-inf) = vals(2);
end
imgStack = log1p(permute(imgStack,[2 1 3])); % transpose the SPic variable so it's same as BImg

% A spectral image for all the different frequencies
colors = flip(jet(size(imgStack,3)));


rgb = CreateRGB2_mat(imgStack, colors, norma, wbalance);


if plotter
    imagesc(rgb)
    % Set the correct colorbar for this image
    colormap(colors)

    h = colorbar;
    hTicks = linspace(colorbarVals(1),colorbarVals(end),length(h.Ticks));
    h.Ticks = h.Ticks; % prevent ticks from changing
    h.TickLabels = num2str(hTicks', '%2.2f');
    ylabel(h, 'spectral frequency', 'FontSize',12)
end

