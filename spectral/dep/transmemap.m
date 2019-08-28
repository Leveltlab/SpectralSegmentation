function [sbxt, dim, freq] = transmemap(varargin) 
%[sbxt, dim, freq] = transmemap(varargin) 
%returns transposed stack as memory mapped file


if nargin > 0   
    filenm = varargin{1};
else
   [fn, pn] = uigetfile('*Trans*.dat');
   filenm = [pn fn];
end
    %load transposed data
   [dim, freq] = Zgetinfo(filenm);    
   sbxt = memmapfile(filenm, 'Format', {'uint16' [dim(1) dim(2) * dim(3)], 'y'}, 'Offset', 500);