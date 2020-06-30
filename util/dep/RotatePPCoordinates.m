function PP = RotatePPCoordinates(PP, rotation, dims)
% PP = RotatePPCoordinates(PP, rotation, dims)
% Rotate ROI contour coordinates from PP struct
% 
% Input: 
%   PP: struct with fields:
%       PP.Con: which has is a struct array with fields
%           PP.Con.x, PP.Con.y. Which hold the contour coordinates 
%           from all ROIs in a dataset.
%       PP.P: a matrix where the first row has the x coordinates for
%           the seed points from all the ROIs, and a second row for 
%           the y coordinates for the seedPoints from all the ROIs.
%       PP.Cnt: number of ROIs.
%   rotation: double. 
%       Rotation that has to be applied to the coordinates in degrees.
%   dims: dimension of the image that belongs to the contour coordinates
%       2 element double [height x width]
%
% Output:
%   PP: struct with updated coordinates
% 
% 
% Leander de Kraker
% 2020-5-22
% 

shift = dims./2;
rotRad = deg2rad(rotation);

for j = 1:PP.Cnt
    % shift the coordinate points so the centre coordinate is 0
    x = PP.Con(j).x - shift(1); % Adjust the contour coordinates
    y = PP.Con(j).y - shift(2);
    % Rotate the shifted coordinates (around the new centre [0, 0])
    x =  x  * cos(rotRad) - y  * sin(rotRad);
    y =  x  * sin(rotRad) + y  * cos(rotRad);
    % Shift them back where they should be
    PP.Con(j).x = x + shift(1);
    PP.Con(j).y = y + shift(2);
end

% Also adjust the contour 'centres' / seed points
xp = PP.P(1,:) - shift(1);
yp = PP.P(2,:) - shift(2);
xp = xp * cos(rotRad) - yp * sin(rotRad);
yp = xp * sin(rotRad) + yp * cos(rotRad);
PP.P(1,:) = xp + shift(1);
PP.P(2,:) = yp + shift(2);

end

