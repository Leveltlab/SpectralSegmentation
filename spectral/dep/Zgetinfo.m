function [Dim, Freq, Type] = Zgetinfo(strfn)

% fid = fopen(strfn);
% Header = fread(fid, [1, 124], 'char');
% fclose(fid);
Fl = memmapfile(strfn);
Header = char(Fl.Data(1:64))';
strh = strsplit(Header,{' ', '\0'});

%idx = 1:min(length(Header), 100); % 
%strh = strsplit(Header(idx),{' ', '\0'});

datestr = [strh{5} ' ' strh{6}];

Dim = [str2double(strh{1}) str2double(strh{2}) str2double(strh{3})];
Freq = str2double(strh{4});
Type = strh{7};


disp(datestr)
disp(['Dimensions : ' num2str(Dim)])
disp(['Frequency :  ' num2str(Freq)])
