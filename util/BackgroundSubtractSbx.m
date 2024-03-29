function BackgroundSubtractSbx(varargin)
% Circular or gaussian filtered background subtraction to improve
% vignetting and blurriness (for 1 photon data especially necessary)
% BackgroundSubtractSbx(fileIn, fileOut, fradius, method, shift, smoothDim, smoothSe, plotter)
% 
% 
% The code can be executed either as a function or as a script. If no 
% input is given prompts will ask for file and standard values will be used.
% 
% Optionally also does gaussian smoothing over 1 dimension to slightly
% decrease vertical or horizontal banding noise. Via variables smoothDim 
% & smoothSe
% 
% 
% Input (optional):
%   1 fileIn: string. filepath and filename of the Sbx file to load
%   2 fileOut:string. filepath and filename of the Sbx file to save results
%   3 fRadius:double [1x1]. the radius of the filter. Standard 55. 
%            final size of filter is fRadius * 2 + 1
%   4 method: string. which method? 'Circular Average' (default)
%                                   'Gaussian Average'
%   5 shift: double [1x1]. How much to shift the background, to prevent 
%           too much subtraction leading to underexposure
%   6 SmoothDim: double [1x1]. over which dimension to apply image 
%             smoothing. 0 (default) No smoothing
%                        1 horizontal smoothing to remove vertical banding
%                        2 vertical smoothing to remove horizontal banding
%   7 smoothSe: double [1x1]. Size of the smoothing gaussian. default 1.85
%   8 plotter: boolean. false (default). plot end results of final frame? 
% 
%
% Output: Saves a background subtracted .sbx file with fileOut name
%       in the accompanying .mat file there will be info about the
%       filtering parameters used in info.filtered
%           
% 
% Leander de Kraker
% 2020-11-9
% 
%% Use the input or ask files and use standard values

clearvars info
clearvars -global info


fRadius = [2, 25];
shift = 10000; % Shift background prevents underexposure
method = 'Donut Average'; % possiblities 'Gaussian Average' | 'Circular Average' | 'Donut Average'
subsample = 4;  % rescaling the data images for faster processing
smoothDim = 0;
smoothSe = 1.75;
plotter = false;

% How deep into the img edge to take data to average for padding values
padTaker = 15; 


if exist('varargin', 'var') && nargin >= 2 % In case both sbx and spsig file name are given
    fileIn = varargin{1};
    fileOut = varargin{2};
else % No files gives as input/ executing as script: ask for files
    [fileNameIn, filePathIn] = uigetfile('.sbx');
    fileIn = [filePathIn, fileNameIn(1:end-4)];
    [fileNameOut, filePathOut] = uiputfile([fileIn, '-Sub.sbx']);
    fileOut = [filePathOut, fileNameOut(1:end-4)];
    if all(fileNameOut == 0) || all(fileNameIn == 0) % if nothing was selected quit
        fprintf('No file selected!\n')
        return
    end
end

% Check filter Radius
if exist('varargin', 'var') && nargin >= 3
    fRadius = varargin{3};
end

% Check filter method
if exist('varargin', 'var') &&  nargin >= 4 && ~isempty(varargin{4})
    method = varargin{4};
end

% Shift background values so not too many values become 0
if exist('varargin', 'var') &&  nargin >= 5 && ~isempty(varargin{5})
    shift = varargin{5};
end

% Over which dimension to do the 1D smoothing to remove banding noise
if exist('varargin', 'var') && nargin >= 6 && ~isempty(varargin{6})
    smoothDim = varargin{6};
    if ~isempty(smoothDim) && smoothDim > 0
        smoothDim = varargin{6};
    end
end

% How big should that 1D filter be (standard deviation)
if exist('varargin', 'var') && nargin >= 6 && ~isempty(varargin{7})
    smoothSe = varargin{7};
end

% Plot final frame?
if exist('varargin', 'var') &&  nargin == 8
    plotter = varargin{8};
end


% Create filter using specified method
switch method
    case 'Circular Average'
        filt = fspecial('disk', fRadius);
        
    case 'Gaussian Average'
        sigma = fRadius/3;
        filt = fspecial('gaussian', fRadius*2+1, sigma);
        
    case 'Donut Average'
        filtInner = fspecial('disk', fRadius(1));
        filt = fspecial('disk', fRadius(2));
        filtInner = filtInner ./ ...
            (filtInner(fRadius(1)+1, fRadius(1)+1) ./ filt(fRadius(2), fRadius(2)));
        filtSize = size(filt, 1);
        fBuf = fRadius(2)-fRadius(1);
        idx =fBuf+1:filtSize-fBuf;
        filt(idx, idx) = filt(idx,idx)-filtInner;
        filt = filt ./ sum(filt(:));
        fRadiusInner = fRadius(1);
        fRadius = fRadius(2);
        
    otherwise
        warning('no valid method')
        return
end

%% Do the background subtraction
fileID = fopen([fileOut '.sbx'], 'w'); % open sbx file

% Thoroughly delete global info variable to prepare for the new file
clearvars -global info
clearvars info

imgOri = sbxread(fileIn, 0, 1);
global info
imgOriDim = size(imgOri);
img = imresize(imgOri, 1/subsample);
imgDim = size(img);
imgPadded = zeros(imgDim+fRadius*2);

imgOriginalMax = zeros(imgOriDim);
backgroundMax = zeros(imgOriDim);

tic
for i = 0:info.max_idx-1
    imgOri = sbxread(fileIn, i, 1);
    img = imresize(imgOri, 1/subsample);
%     imgPadded = padarray(img, [fRadius fRadius], 'replicate');
    
    % Creating nice padding that has averaged values of the outer edges of img, not only the most outer value
    imgPadded(1:fRadius, fRadius+1:end-fRadius)         = repmat(mean(img(1:padTaker, :)), [fRadius, 1]); % top
    imgPadded(end-fRadius+1:end, fRadius+1:end-fRadius) = repmat(mean(img(end-padTaker:end, :)), [fRadius, 1]); % bottom
    imgPadded(fRadius+1:end-fRadius, 1:fRadius)         = repmat(mean(img(:,1:padTaker), 2), [1, fRadius]); % left
    imgPadded(fRadius+1:end-fRadius, end-fRadius+1:end) = repmat(mean(img(:,end-padTaker:end), 2), [1, fRadius]); % right
    imgPadded(1:fRadius, 1:fRadius)                 = mean(imgPadded(fRadius+1:fRadius+padTaker,1)); % top left corner
    imgPadded(1:fRadius, end-fRadius+1:end)         = mean(imgPadded(fRadius+1:fRadius+padTaker,end)); % top right corner
    imgPadded(end-fRadius+1:end, 1:fRadius)         = mean(imgPadded(end-fRadius-padTaker:end-fRadius, 1)); % bottom left corner
    imgPadded(end-fRadius+1:end, end-fRadius+1:end) = mean(imgPadded(end-fRadius-padTaker:end-fRadius, end)); % bottom right corner
    imgPadded(fRadius+1:end-fRadius, fRadius+1:end-fRadius) = img; % Put in the actual image
    
    background = conv2(imgPadded, filt, 'same');
    background = background(fRadius+1:end-fRadius, fRadius+1:end-fRadius);
    background = imresize(background, imgOriDim);
    imgCorrected = uint16(double(imgOri) - (background - shift));
    
    backgroundMax = max(backgroundMax, background);
    imgOriginalMax = max(imgOriginalMax, double(imgOri));
    
    % Report status and check for over/under exposure
    if mod(i+1,1000)==0
        ETA = ((info.max_idx-i) / (i/toc)) / 60;
        fprintf('done with frame %d/%d: elapsed %.2f minutes, ETA %.2f minutes\n',...
                i+1, info.max_idx, toc/60, ETA)
        % warning if more than 1% of pixels are 0 / underexposed
        if sum(imgCorrected(:) == 0) > (numel(imgCorrected)/100) 
            warning('%.1f%% of pixels underexposured, increase the variable shift',...
                imgCorrected == 0 ./ numel(imgCorrected) * 100)
        end
        % Warning if more than 0.5% of pixels are 65536, 2^16, overexposed
        if sum(imgCorrected(:) >= 2^16 - 1) > (numel(imgCorrected) / 200)
            warning('%.1f%% of pixels overexposed, decrease the variable shift or/and check original data for overexposure',...
                imgCorrected >= 2^16 ./ numel(imgCorrected) * 100)
        end
    end
    
    % Banding diminishing & correlation increasing smooth
    if smoothDim == 1 % horizontal smoothing
        % Simple smoothing to reduce CMOS sensor vertical banding noise
        imgCorrected = smoothG(imgCorrected', smoothSe)';
    elseif smoothDim == 2 % vertical smoothing
        % Simple smoothing to reduce CMOS sensor horizontal banding noise
        imgCorrected = smoothG(imgCorrected, smoothSe);
    elseif smoothDim == 3
        % Simple smoothing in both directions to reduce noise
        imgCorrected = smoothG(smoothG(imgCorrected, smoothSe)', smoothSe)';
    end
    
    fwrite(fileID, imgCorrected', 'uint16');
end
fprintf('Completed background subtraction\n')
fclose(fileID);

% Save .mat file for .sbx file, with info about background subtraction
info.filtered.method = method;
info.filtered.size = fRadius;
if strcmp(method, 'Donut Average')
    info.filtered.sizeGap = fRadiusInner;
end
info.filtered.resized = 1/subsample;
info.filtered.ValueShift = shift;
info.filtered.filt = filt;
info.filtered.padTaker = padTaker;
info.filtered.smoothSe = smoothSe;
info.filtered.smoothDim = smoothDim;
info.filtered.backgroundMax = backgroundMax;
info.filtered.dataMax = imgOriginalMax;
save([fileOut '.mat'], 'info');

% Plot last frame
if plotter
    % Activate much tighter subplots
    % [subplot margin top&side],[figure bottomspace,topspace],[leftspace,rightspace]
    subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.04], [0.07 0.1], [0.07 0.01]);
    figure('Units', 'normalized', 'Position', [0.3 0.4 0.6 0.4])
    subplot(1,3,1)
    imagesc(img); title('last frame of data'); colorbar
    subplot(1,3,2)
    imagesc(background); title('background for last frame'); colorbar
    subplot(1,3,3)
    imagesc(imgCorrected); title('background subtracted frame'); colorbar
    colormap(gray)
    figtitle(fileOut)
end


