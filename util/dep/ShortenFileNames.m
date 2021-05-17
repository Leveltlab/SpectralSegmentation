function filenamesShort = ShortenFileNames(filenames, parts, varargin)
% Just shorten the filename. The filenames are split on '_' or ' '
% 
% Input: 
%   filenames: [n x 1] cell array with strings
%   parts: 1D double, saying which parts of the filename to take
% Optional: [n x 1] datetime array. Last part of filename will be the date
% in yyyy-mm-dd format
% 
% Output: filenamesShort: [n x 1] cell array with strings that have just
% the first part of the filename and the date, or just the requested parts
% of the filename
% 
% Leander de Kraker
% 2021-5-14
%

nfiles = length(filenames);
filenamesShort = cell(nfiles, 1);
nparts = length(parts);

for i = 1:nfiles
    delims = regexp(filenames{i}, '_');
    delims2 = regexp(filenames{i}, ' ');
    delims = sort(unique([1, delims, delims2, (length(filenames{i})+1)]));
    for j = 1:nparts
        if length(delims)>= max(parts)
            filenamesShort{i} = [filenamesShort{i}, filenames{i}(delims(parts(j)):delims(parts(j)+1)-1)];
        else
            fprintf('filename part %d not found for name %d', parts(j), i)
        end
    end
    if nargin == 3 % Add date if given
        filenamesShort{i} = [filenamesShort{i}, ' ', datestr(varargin{1}(i), 'yyyy-mm-dd')];
    end
    filenamesShort{i} = strrep(filenamesShort{i}, '_', ' ');
end