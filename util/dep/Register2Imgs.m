function [image2, transformation] = Register2Imgs(reference,image2,varargin)
% Register the second image to the reference image, cut and/or pad image 2
% to the same size as reference
% Uses the 2D cross-correlation then the 'check for best rotation'
% technique.
%
%
% input: image1: this 2D matrix is the reference
%        image2: this 2D matrix will be moved
%
% output: image2: 2D double with the size of image1. Registration
%           result
%         transformation: which x shift, y shift (in pixels) was applied.
%                         which rotation was applied(in degrees)
%
% Leander de Kraker
% 2018-10-3


if exist('varargin','var') && nargin==3
    interpMethod = varargin{1};
else
    interpMethod = 'nearest';
end


dims = size(reference);

% Remove -infs if there are any
vals = unique(reference);
if vals(1) == -inf
    reference(reference<vals(2))= vals(2);
end
% Scale the values of the images
reference = (reference - mean(reference(:))) / range(reference(:));
img2mean = mean(image2(:));
img2range = range(image2(:));

% Crop the image that has to be moved to make the proces faster
if min(dims)>150
    buffer = 30; % Hardcoded buffer to remove from the image
else
    buffer = 1;
end

% Register the chronic dataset to this dataset
correl = xcorr2_fft(reference,...
                    (image2(buffer:dims(1)-buffer, buffer:dims(2)-buffer)-img2mean)/img2range);
[~,snd] = max(correl(:));
[ij, ji] = ind2sub(size(correl),snd);
x = dims(1) - ij - buffer;
y = dims(2) - ji - buffer;


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
if size(image2,1) > dims(1) % Cut bottom
    image2(dims(1)+1:end,:) = [];
else % padd bottom
    image2(dims(1),1) = 0;
end
if size(image2,2) > dims(2) % Cut right side
    image2(:,dims(2)+1:end) = [];
else % padd right side
    image2(1,dims(2)) = 0;
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
