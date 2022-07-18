function strOut = FindNNumInStr(str, n)
% Find a number that is exactly n digits long in a string
% 
% 
% Leander de Kraker
% 2022-7-15
% 

% Get info about the digits in the string
pos = regexp(str, '\d');
a=diff(pos);
b=find([a inf]>1);
c=diff([0 b]); % length of the sequences
d=cumsum(c);   % endpoints of the sequences
good = pos(d(c==n));

strOut = str(good-n+1:good);