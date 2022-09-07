function I = ZgetXY(zx, nm, sbxt, dim)
%return nm images given z coordinate
% Chris van der Togt, 2017, 
% Netherlands Institute for Neuroscience 

if ~isobject(sbxt)
    return
end
Stacksz = dim(1);
Colsz = dim(2);
Rowsz = dim(3);

Zx = (1:Stacksz:Stacksz*Colsz*Rowsz)+zx-1;

if nm > 1
    nx = (0:nm-1)';
    Imgx = repmat(Zx, nm, 1); %indices for each image
    Imgx = Imgx + repmat(nx, 1, length(Zx)); %indices for consecutive images
else
    Imgx = Zx;
end

I = sbxt.Data(Imgx(:));
I = reshape(I, nm, Colsz, Rowsz);
I = permute(I, [2 3 1]); %z order images