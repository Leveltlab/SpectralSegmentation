function pos = GetCenteredFigPos(desiredSize)
% pos = GetCenteredFigPos(desiredSize)
% Get pixel positions for a figure that is centered on your screen with
% desired size
% 
% input: 
%   desiredSize ([2 x 1] double): width and height in pixels of desired 
%                                 figure size limited by screen size
% 
% output:
%   pos ([1 x 4] double): position of figure
%  
% 
% Leander de Kraker
% 2022-7-12
% 

scrnSz = get(0,'ScreenSize');
desiredSz = [min(desiredSize(1), scrnSz(3)),...
             min(desiredSize(2), scrnSz(4))];
p1 = round(scrnSz(3)/2)-round(desiredSz(1)/2); 
p2 = round(scrnSz(4)/2)-round(desiredSz(2)/2);
pos = [p1 p2 desiredSz(1) desiredSz(2)-25];