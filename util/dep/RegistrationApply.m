function outData = RegistrationApply(inData, transformation, varargin)
% Apply the transformation as created by Register2Imgs.m to either an image
% Or to contours data as described in the variable PP from getspectrois
% 
% input: 
%   inData: 2D or 3D matrix that needs to be registered
%   tranformation: struct with fields x, y, and rotAng
%               x = vertical shift (1st dimension)
%               y = horizontal shift (2nd dimension)
%               rotAng = rotation in degrees (positive = clockwise)
%   interpMethod (optional): interpolation method for rotation: 'bilinear', 
%                                           'linear', 'nearest' (default)
% 
% See also: RegistrationApplyPP, RotatePPCoordinates
% 
% Leander de Kraker
% 2019-7-31
%

if exist('varargin','var') && nargin==3
    interpMethod = varargin{1};
else
    interpMethod = 'nearest';
end

x = transformation.x;
y = transformation.y;
r = transformation.rotAng;

dims = size(inData);

outData = zeros(dims);
% Apply the transformation to the image
if x > 0
    if y > 0
        outData(1+x:end, 1+y:end, :) = inData(1:end-x, 1:end-y, :);
    else
        outData(1+x:end, 1:end+y, :) = inData(1:end-x, 1-y:end, :);
    end
else % if xoffset  < 0
    if y > 0                    
        outData(1:end+x, 1+y:end, :) = inData(1-x:end, 1:end-y, :);
    else
        outData(1:end+x, 1:end + y, :) = inData(1-x:end, 1-y:end, :);
    end
end
if r ~= 0
    outData = imrotate(outData, r, interpMethod, 'crop');
end

