function h = ImageValueText(img, color, rounding)
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

h = text(x, y, valsText, 'Color', color, 'HorizontalAlignment', 'center', 'FontName', 'calibri');


