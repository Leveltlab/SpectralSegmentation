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
good = pos(b(c==n));

if isscalar(good)
    strOut = str(good-n+1:good);
else
    strOut = [];
end



