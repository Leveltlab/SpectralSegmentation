function t =  ETA(i, iend, telapsed)
% Calculated the Estimated Time to Arrival of your executing code
% 
% input: 
%   i (scalar double):        current iteration
%   iend (scalar double):     number of total iterations to do
%   telapsed (scalar double): how long it took to get to i
%
% output:
%   t (scalar double): time to completion probably
% 
% Leander de Kraker
% 2023-2-9
% 

t = (iend-i) / (i/telapsed);