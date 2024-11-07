function info = RequestRecInfoProcess(info, hz, scaleUm, FOVum, pixelAspectRatio, squareFOV)
% Put the answers from RequestRecInfo into info
% 
% Input: 
%   - info (struct): with field info.Shap, which contains width and height
%                                          of recording in pixels
%   - hz (double): Imaging frequency
%   - scaleUm (double or empty): the width of one pixel in micrometer (um)
%   - FOVum   (double or empty): the width of the entire FOV in micrometer
%   - pixelAspectRatio (double or empty): How much higher a pixel is than
%           it is wide, so a square image with more pixels for the
%           horizontal direction results in >1.
%           example: a square image recorded with pixels so aspect ratio is
%           16:9 (width:height) results in pixel aspect ratio 16/9 = 1.7778
%   - squareFOV (boolean): Was the imaged field square? True | false
%           if true, then the pixel aspect ratio is the width/height of the
%           image.
%   
% Output: 
%   - info (struct): with new fields: Freq, 
%                                     scaleUm (optional),
%                                     pixelAspectRatio (optional)
% 
% See also: RequestRecInfo
% 
% Leander de Kraker
% 2024-10-30

if ~isempty(hz)
    if ~isfield(info, 'Freq')
        info.Freq = hz;
    elseif info.Freq ~= hz
        warning('info.Freq was already present. %.2fHz. Overwritten with %.2fHz',...
            info.Freq, hz)
    end
end

if ~isempty(scaleUm)
    info.scaleUm = scaleUm;
elseif ~isempty(FOVum) && isfield(info, 'Shape')
    info.scaleUm = FOVum / info.Shape(1);
end

if ~isempty(pixelAspectRatio)
    info.pixelAspectRatio = pixelAspectRatio;
elseif squareFOV && isfield(info, 'Shape')
    info.pixelAspectRatio = info.Shape(1) / info.Shape(2);
end

