function extend = cmapLimsAround0(img, n)
% extend = cmapLimsAround0(img, n)
% Get how many colors should be used for which color to center the colors 
% around the value 0 of the img
% 
% input:
%   - img (xD double): data to base colormap around
%   - n (number of colors in colormap)
% 
% Output:
%   - extend ([n-1 x 1] double): amount of colors to take for each
% 
% Example:
% img = randn(25,25)+2;
% colors1 = [0 1 0; 0 0 0; 1 0 0];
% n1 = 3;
% 
% n = cmapLimsAround0(img, n1); % Get how many colors to take for each color
% colorsMap = cmapL(colors1, n); % Use that to get the colormap
% imagesc(img)
% colormap(colorsMap);
% colorbar
% % Another set of colors
% colors1 = 'blue roast'; % italian roast has 7 colors
% n1 = 7;
% n = cmapLimsAround0(img, n1); % Get how many colors to take for each color
% colorsMap = cmapL(colors1, n);
% colormap(colorsMap)
% 
% Leander de Kraker
% 2024-4-23
% 

extend = [max(img(:)), min(img(:))];
extend = abs(round(extend./abs(diff(extend)).*256));
if n > 3
    n = floor(n/2);
    extend = floor([repmat(extend(1), [1 n]), repmat(extend(2),[1 n])]./n);
end