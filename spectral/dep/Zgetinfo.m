function [Dim, Freq, Type] = Zgetinfo(strfn)
% Chris van der Togt, 2017, 
% Netherlands Institute for Neuroscience 

Fl = memmapfile(strfn);
Header = char(Fl.Data(1:64))';
strh = strsplit(Header,{' ', '\0'});
datestr = [strh{5} ' ' strh{6}];

Dim = [str2double(strh{1}) str2double(strh{2}) str2double(strh{3})];
Freq = str2double(strh{4});
Type = strh{7};


disp(datestr)
disp(['Dimensions : ' num2str(Dim)])
disp(['Frequency :  ' num2str(Freq)])
