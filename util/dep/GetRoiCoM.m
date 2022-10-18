function com = GetRoiCoM(Mask)
% com = GetRoiCoM(Mask)
% Get the Center of Masses (CoM) from ROIs in a Mask
% 
% Input:
%   Mask (2D double): field of view with numbers indicating which pixels
%                     are occupied by which ROI
% 
% Output:
%   com ([n x 2] double): Center of Masses of all n ROIs inside the Mask.
%                       first column = x values, second = y values.
% 
% Leander de Kraker
% 2022-4-10
% 

cnt = max(Mask(:));
com = zeros(2, cnt);
for i = 1:cnt
    [y, x] = find(Mask == i);
    com(1, i) = mean(x);
    com(2, i) = mean(y);
end