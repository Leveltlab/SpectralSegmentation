function [scaleBarx, scaleBary] = BurnScaleBarIdx(dims, len)
% At which pixel idx to burn the image to show scalebar
% 
% Leander de Kraker
% 2024-12-20
% 

if len>0
    scaleBary = round(dims(1)*0.9); % Set vertical position of scalebar
    scaleBary = scaleBary:scaleBary+2; % Set thickness of scalebar
    scaleBarx = round(dims(2)*0.9); % Setting the x position of scalebar
    scaleBarx = scaleBarx-round(len)+1:scaleBarx;
else
    scaleBary = [];
    scaleBarx = [];
end