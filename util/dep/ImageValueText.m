function h = ImageValueText(img, color, rounding, ha)
% ImageValueText(img, color, rounding, ha)
% Input:
%   - img ([m x n] double)
%   - color (color for the text)
%   - rounding ([1 x 1] double): to how many decimals to round data
%   - ha ([1 x 1] axes handle): Where to plot the text
% 
% Leander de Kraker
% 2023-5-30
% 

dims = size(img);
x = repmat(1:dims(2), [dims(1), 1]);
y = repmat(1:dims(1), [dims(2), 1])';
x = x(:);
y = y(:);

valsText = round(img(:), rounding);
valsText = num2str(valsText);
valsText(sum(ismember(valsText, 'Na'),2)==3,:) = ' '; % replace nans with empty

h = text(ha, x, y, valsText, 'Color', color, 'HorizontalAlignment', 'center', 'FontName', 'calibri');


