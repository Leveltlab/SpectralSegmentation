function [sbxt, dim, freq] = transmemap(varargin) 
% [sbxt, dim, freq] = transmemap(filename)  returns transposed stack as 
% a memory mapped file
%
% input: string that contains the path- and filename of the transposed data
% When no input is given, the script will ask for the transposed data file
% 
% output: 
%   sbxt: memorymapped file
%   freq (double): the sampling rate of the data in Hz
%   dim: The dimensions of the data, [time, X, Y]
% 
% Transposed data can be decimated or non-decimated, whether the data is
% decimated should be in the filename! The filename should contain
% the text: DecTrans
% Decimated data is in data format double, and will not be subsampled by
% other scripts.
% 
% Chris v.d. Togt
% edited by Leander to include decimate recognition 2020-8-4
% 

% What file to ananlyse
if nargin > 0
    filenm = varargin{1};
else
   [fn, pn] = uigetfile('*Trans.dat');
    filenm = [pn fn];
end

% Find out whether the trans file contains decimated data
if nargin == 2
    decimated = varargin{2};
else
    % Find out whether the data is decimated
    decimated = false;
    strPos = regexp(filenm, 'Trans');
    if (strPos - 3) >= 1
        if strcmp(filenm(strPos-3:strPos-1), 'Dec')
            decimated = true;
        end
    end
end

%load transposed data
if decimated
    [dim, freq] = Zgetinfo(filenm); 
    sbxt = memmapfile(filenm, 'Format', {'double' [dim(1) dim(2) * dim(3)], 'y'}, 'Offset', 500);
    
else
   [dim, freq] = Zgetinfo(filenm);    
   sbxt = memmapfile(filenm, 'Format', {'uint16' [dim(1) dim(2) * dim(3)], 'y'}, 'Offset', 500);
   
end



