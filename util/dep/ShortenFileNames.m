function filenamesShort = ShortenFileNames(filenames, parts, varargin)
% Shorten the filename. The filenames are split on '_' or ' ' and on '\' to
% split the filepath
% 
% filenameShort = ShortenFileNames{'bla\la\mousename_recNr_lol.mat',...
%                                                   [-1 1], datetimeValue)
% 
% Input: 
%   filenames: [n x 1] cell array with strings
%   parts:  [n x 1] double, saying which parts of the filename to take.
%                       Positive numbers: part of the filename
%                       Negative numbers: parts of the pathname (counting
%                                       from the end of the path name)
%                       zero: the optional filedate 
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
    delims0 = regexp(filenames{i}, '\');
    if ~isempty(delims0)
        pn = filenames{i}(1:(delims0(end)-1));
        fn = filenames{i}((delims0(end)+1):end);
    else
        fn = filenames{i};
        pn = [];
    end
    delims1 = regexp(fn, '_');
    delims2 = regexp(fn, ' ');
    delims1 = sort(unique([1, delims1, delims2, (length(fn)+1)]));
    for j = 1:nparts
        if length(delims1)>=(parts(j)+1) && parts(j)>0
            filenamesShort{i} = [filenamesShort{i}, fn(delims1(parts(j)):delims1(parts(j)+1)-1)];
        elseif parts(j)<0 && ~isempty(delims0)
            filenamesShort{i} = [filenamesShort{i}, pn(delims0(end+parts(j)):delims0(end+parts(j)+1)-1)];
        elseif parts(j)==0
            filenamesShort{i} = [filenamesShort{i}, ' ', datestr(varargin{1}(i), 'yyyy-mm-dd')];
        else
            fprintf('filename part %d not found for name %d', parts(j), i)
        end
    end
    filenamesShort{i} = strrep(strrep(filenamesShort{i}, '_', ' '), '\', ' ');
end

