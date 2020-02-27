function [imageOut, transformation] = Register2Imgs(reference,image2)
    % Register the second image to the first image, cut and/or pad image 2
    % to the same size as image 1
    % Uses the 2D cross-correlation then the 'check for best rotation'
    % technique.
    %
    %
    % input: image1: this 2D matrix is the reference
    %        image2: this 2D matrix will be moved
    %
    % output: imageOut: 2D double with the size of image1. Registration
    %           result
    %         transformation: which x shift, y shift (in pixels) was applied. 
    %                         which rotation was applied(in degrees)
    %
    % Leander de Kraker
    % 2018-10-3
    
    
    dims = size(reference);
    
    % Remove -infs if there are any
    vals = unique(reference);
    if vals(1) == -inf
        reference(reference<vals(2))= vals(2);
    end
    reference = (reference - mean(reference(:))) / range(reference(:));
    image2 = (image2 - mean(image2(:))) / range(image2(:));
    
    % Register the chronic dataset to this dataset
    if min(dims)>150
        buffer = 30; % Hardcoded buffer to remove from the image
    else
        buffer = 1;
    end
    im = image2(buffer:dims(1)-buffer, buffer:dims(2)-buffer);
    im = im - mean(im(:));
    correl = xcorr2(reference,im);
    [~,snd] = max(correl(:));
    [ij, ji] = ind2sub(size(correl),snd);
    x = size(reference,1) - ij - buffer;
    y = size(reference,2) - ji - buffer;
    % Cutting or padding the chornic image if necessary
    if x>=1 % Cut top
        image2 = image2(x:end, :);
    else % Padd top % Maybe add a +1; code untested because condition doesn't happen often
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
    simil  = zeros(1, length(rotation));
    % correl = zeros(nfiles, length(rotation));
    for j = 1:length(rotation)
        BImgRot = imrotate(image2, rotation(j),'bicubic','crop');
        simil(j)= corr2(BImgRot, reference);
    end

    % Best rotations
    [~, rotIdx] = max(simil);
    rotAng = rotation(rotIdx);

    imageOut = imrotate(image2, rotAng,'bicubic', 'crop');
    transformation = struct('rotAng',rotAng,'x',x,'y',y);
end

