function [image2, transformation, correl] = Register2ImgsNew(reference,image2,varargin)
% Register the second image to the first image, cut and/or pad image 2
% to the same size as image 1
% Uses the 2D cross-correlation then the 'check for best rotation'
% technique.
% 
% 
% input: reference: this 2D matrix is the reference
%        image2: this 2D matrix will be moved
%        (optional)interpMethod: interpolation for rotation(default=nearest)
%
% output: image2: 2D double with the size of image1. Registration
%           result
%         transformation: which x shift, y shift (in pixels) was applied.
%                         which rotation was applied(in degrees)
% 
%  
% See also, Register2Imgs. This new version can deal with different image
% sizes but I still have to check if that goes right in every case.
% 
% Leander de Kraker
% 2018-10-3


if exist('varargin','var') && nargin==3
    interpMethod = varargin{1};
else
    interpMethod = 'nearest';
end

dimsRef = size(reference);
dims2 = size(image2);

% Scale the values of the images
reference = (reference - mean(reference(:))) ./ range(reference(:));
image2 = (image2 - mean(image2(:))) ./ range(image2(:));

% Crop the image that has to be moved to make the proces faster
if min(dims2)>150
    buffer = 30; % Hardcoded buffer to remove from the image
else
    buffer = 1;
end

% Register the chronic dataset to this dataset
correl = xcorr2_fft(reference, image2(buffer:dims2(1)-buffer, buffer:dims2(2)-buffer));
[~,snd] = max(correl(:));
[ij, ji] = ind2sub(size(correl),snd);
x = dims2(1) - ij - buffer;
y = dims2(2) - ji - buffer;

% Cutting or padding the chornic image if necessary
if x>=1 % Cut top
    image2 = image2(x:end, :);
else % Padd top
    image2(-x+1:-x+size(image2,1),:) = image2;
end
if y>=1 % Cut left beginning
    image2 = image2(:, y:end);
else % Padd left beginning
    image2(:,-y+1:-y+size(image2,2)) = image2;
end
if size(image2,1) > dimsRef(1) % Cut bottom
    image2(dimsRef(1)+1:end,:) = [];
else % padd bottom
    image2(dimsRef(1),1) = 0;
end
if size(image2,2) > dimsRef(2) % Cut right side
    image2(:,dimsRef(2)+1:end) = [];
else % padd right side
    image2(1,dimsRef(2)) = 0;
end

% Rotation by checking rotations
rotation = -1:0.1:1; % rotations to apply: counterclockwise to clockwise degrees
correl  = zeros(1, length(rotation));
for i = 1:length(rotation)
    BImgRot = imrotate(image2, rotation(i), interpMethod,'crop');
    correl(i)= corr2(BImgRot, reference);
end

% Best rotations
[~, rotIdx] = max(correl);
rotAng = rotation(rotIdx);
if rotAng ~= 0
    image2 = imrotate(image2, rotAng, interpMethod, 'crop');
end

transformation = struct('rotAng',rotAng,'x',x,'y',y);

