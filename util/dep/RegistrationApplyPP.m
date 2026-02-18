function PP = RegistrationApplyPP(PP, transformation, dims)
% Apply the transformation as created by Register2Imgs.m to either an image
% Or to contours data as described in the variable PP from getspectrois
% 
% input: 
%   tranformation: struct with fields x, y, and rotAng
%               x = vertical shift
%               y = horizontal shift
%               rotAng = rotation in degrees (positive = clockwise)
%   inData: struct with fields:
%             P, where P(1,i) = x coordinate of important point in ith ROI,
%                      P(2,i) = y coordinate of important point in ith ROI
%             Con, a struct array with x and y values describing contours
%             Cnt, number of ROIs
%   dims: [1x2 double] the size of the image that was being registered, to
%         which these contours belong.
%
% See also: RegistrationApply
% 
% Leander de Kraker
% 2020-6-11
% 

x = transformation.x;
y = transformation.y;

% Translate
for i = 1:PP.Cnt
    % CONTOUR TRANSLATION
    PP.Con(i).x = PP.Con(i).x + y;
    PP.Con(i).y = PP.Con(i).y + x;
end
% Also adjust the contour 'centres'
PP.P(1,:) = PP.P(1,:) + y;
PP.P(2,:) = PP.P(2,:) + x;

% Rotate
PP = RotatePPCoordinates(PP, transformation.rotAng, dims);

