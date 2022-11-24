function img = MidtoneBalance(img, m, varargin)
% img = MidtoneBalance(img, m)
% Balancing midtone/ contrast/ brightnesses of an image using a formula
% derived from the Bulirsch-Stoer algorithm. For evaluation of diagonal
% rational function..
% 
% Input: 
%   img (double): data to midtone balance
%   m (scalar, double): value between 0 and 1. 0.5 Means same image gets
%                       outputted. lower values result in brighter image, 
%                       higher values darker image.
%   doScale (scalar, boolean): Scale data between 0 and 1, for accurate
%                           midtone balancing. But perhaps your data values 
%                           is already nearly between 0 and 1, in which 
%                           case you could set soScale to false, for
%                           slightly faster
%
% Leander de Kraker
% 2022-7-20
% 

doScale = true;
if exist('varargin', 'var') && nargin == 3
    doScale = varargin{1};
end
if doScale
    mini = min(img(:));
    maxi = max(img(:));
    img = (img-mini) ./ (maxi-mini); % normalize to range [0 - 1]
end

% Do the midtone balance
img = ((m-1)*img) ./ ((2*m-1)*img - m);

if doScale
    img = (img .* (maxi-mini)) + mini; % Scale values to original limits
end
