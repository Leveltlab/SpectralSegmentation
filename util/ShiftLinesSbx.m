function ShiftLinesSbx(varargin)
% ShiftLinesSbx({strFilePathIn, strFileIn}, doTrans, s, strPathOut)
% ShiftLinesSbx([], false) % select file to correct for horizontal lines
% 
% Shift even or uneven lines to the left, to correct for bi-directional
% scanning mis-alignment
% Also removes very high pixel values 65500 with a more sensible value 
% given the data!
% 
% Input (all optional):
%   - {strFilePathIn, strFileIn}
%   - doTrans (boolean), transpose data? false(default): correct for
%   horizontal lines. True: Correct for vertical lines
%   - shift to apply (scalar value). if empty, tries to determine itself
%   - strPathOut (string): Where to save  output file (if empty, 
%                                               saves in current folder)
% 
% 
% !ONLY APPLY ON NON-MOTION-CORRECTED FILES!! Otherwise the lines are
% shifting around so the correction will be nonsense.
% 
% Leander de Kraker
% 2022-12-1 
% updated 2024-5-3, to implement possible vertical line correction
% 

%% Select file and load
if exist('varargin', 'var') && nargin>1 && ~isempty(varargin{1})
    if length(varargin{1})==2
        filepath = varargin{1}{1};
        filename = varargin{1}{2};
    else
        filepath = '';
        filename = varargin{1};
    end
else
    [filename, filepath] = uigetfile('.sbx', 'get your sbx file');
    filename = filename(1:end-4);
end
pnfn = [filepath filename];

clearvars -global info
clearvars info

load(pnfn) % load info

if ~isfield(info, 'max_idx')
    d = dir([pnfn '.sbx']); 
    switch info.channels
        case 1
            info.nchan = 2;      % both PMT0 & 1
        case 2
            info.nchan = 1;      % PMT 0
        case 3
            info.nchan = 1;      % PMT 1
    end
    nframes =  round(d.bytes/info.sz(2)/info.sz(1)/info.nchan/2);
else
    nframes = info.max_idx;
end

im = sbxread(pnfn, 10, 100);
val = prctile(im(:), 8);


%% Correct for horizontal lines or vertical lines (transpose data)
if exist('varargin', 'var') && nargin >= 2 && ~isempty(varargin{2})
    trans = varargin{2};
else
    trans = false; % Just do horizontal line correction
end


%% Set shift
if exist('varargin', 'var') && nargin >= 3 && ~isempty(varargin{3})
    s = varargin{3};
else % Determine best horizontal line shift
    s = ShiftLinesCheck(im, trans);
end

clearvars idx imS imSH imH1 nshifts shifts buf maxCor stri i d1  x y correl

%% output folder if needed
if exist('varargin', 'var') && nargin == 4
    filepathOut = varargin{4};
else
    filepathOut = '';
end
if ~isempty(filepathOut) & ~strcmp(filepathOut(end), '\')
    filepathOut = [filepathOut '\'];
end

%% Apply shift and save new file

filenameOut = [filepathOut, filename, '_shiftedLines.sbx'];

if s<0
    y = 1;
elseif s>0
    y = 2;
else
    warning('Was not going to apply shift. aborting execution')
    return
end
x = abs(s);


n = 25;
totake = 0:n:nframes;
nbatch = length(totake);
remaining = nframes-totake(end);

toreport = round(nbatch/10);

fileID = fopen(filenameOut, 'w'); % open sbx file

% Write the new file with shited lines
tic
for i = 1:nbatch-1
    im = sbxread(pnfn, totake(i), n);
    if trans
        im = permute(im, [2 1 3 4]);
    end
    im(im>65500) = val; % remove erroneous high values
    im(y:2:end, 1:end-x, :, :) = im(y:2:end, 1+x:end, :, :);
    im =  intmax('uint16')-permute(im, [2 1 3 4]);
    if trans % transpose back
        im = permute(im, [2 1 3 4]);
    end
    fwrite(fileID, im, 'uint16');

    if mod(i, toreport) == 0
        eta =  toc * nbatch / i - toc;
        fprintf('done with frame batch %d/%d in %.1f minutes ETA: %.1f minutes\n',...
                    i, nbatch, toc/60, eta/60)
    end
end
if remaining>0 % Write reamining frames
    im = sbxread(pnfn, totake(end), remaining);
    if trans
        im = permute(im, [2 1 3 4]);
    end
    im(im>65500) = val; % remove erroneous high values
    im(y:2:end, 1:end-x, :, :) = im(y:2:end, 1+x:end, :, :);
    im = intmax('uint16')-permute(im, [2 1 3 4]);
    if trans
        im = permute(im, [2 1 3 4]);
    end
    fwrite(fileID, im, 'uint16');
end
fprintf('Done\n')

fclose(fileID);


% save accompanying mat file
info.shiftedLines = s;
info.shiftedLinesTrans = trans;
if isfield(info, 'simon')
    % I think using sbxread made the data consistent regardless of simon flag
    % So simon can be ignored in the new file
    info = rmfield(info, 'simon');
end
save(filenameOut(1:end-4), 'info')


