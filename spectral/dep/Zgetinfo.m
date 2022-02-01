function [Dim, Freq, Type] = Zgetinfo(strfn)

fid = fopen(strfn);
Header = fgetl(fid);
idx = 1:min(length(Header), 100); % 
strh = strsplit(Header(idx),{' ', '\0'});

datestr = [strh{5} ' ' strh{6}];

Dim = [str2double(strh{1}) str2double(strh{2}) str2double(strh{3})];
Freq = str2double(strh{4});
Type = strh{7};

fclose(fid);

disp(datestr)
disp(['Dimensions : ' num2str(Dim)])
disp(['Frequency :  ' num2str(Freq)])
