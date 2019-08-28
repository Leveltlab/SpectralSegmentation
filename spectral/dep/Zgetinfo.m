function [Dim, Freq] = Zgetinfo(strfn)

fid = fopen(strfn);
Header = fread(fid, 200, '*char');
%Header = Header(2:2:end);
strh = strsplit(Header(1:100)' , {' ','\0'});
datestr = [strh{5} ' ' strh{6}];
Dim = [str2double(strh{1}) str2double(strh{2}) str2double(strh{3})];
Freq = str2double(strh{4});
fclose(fid);

disp(datestr)
disp(['Dimensions : ' num2str(Dim)])
disp(['Frequency :  ' num2str(Freq)])