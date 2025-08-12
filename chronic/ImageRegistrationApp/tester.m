%{  
Fullpath()
IDname()
transform
transformedImage
transformedMask
%}

names = fullfile(registration_data.Path, registration_data.Name);
registration_data.tForm
disp(class(registration_data.Mask_raw{1}))
t = bwlabel(registration_data.Mask_raw{1});

% Apply area opening to remove small connections between large blobs
cleanedImage = bwareaopen(t, 5);  % Remove small regions (adjust threshold)

% Re-label after cleaning
finalLabels = bwlabel(cleanedImage);

coloredImage = label2rgb(finalLabels, 'jet', 'k', 'shuffle');

% Display the colored image
imshow(coloredImage);

bw = registration_data.Mask_raw{1};

% Step 1: Convert binary image to coordinates of non-zero pixels
[rows, cols] = find(bw);

% Combine the row and column indices into a data matrix
data = [rows, cols];

% Step 2: Apply DBSCAN to group the points into clusters (ROIs)
epsilon = 1;  % Maximum distance between points to be considered as neighbors
minPts = 4;   % Minimum number of points to form a cluster
[idx, corepts] = dbscan(data, epsilon, minPts);

% Step 3: Create a labeled image based on DBSCAN results
L = zeros(size(bw));  % Initialize an empty labeled image
for i = 1:max(idx)
    L(sub2ind(size(bw), rows(idx == i), cols(idx == i))) = i;
end

% Visualize the result with different colors for each cluster
coloredImage = label2rgb(L, 'jet', 'k', 'shuffle');
imshow(coloredImage);
title(['Detected ', num2str(max(idx)), ' ROIs using DBSCAN']);

%% 

allMasks = cat(3, registration_data.warpedImg{:});

% Sum over the 3rd dimension to get the combined mask
unifiedMask = sum(allMasks, 3);

coloredImage = label2rgb(unifiedMask, 'jet', 'k', 'shuffle');

imshow(coloredImage);
%% 
function warpedImg = warpImages(transforms, rawImages)
    n_imgs = max(size(rawImages));
    outView = @(x) affineOutputView(size(rawImages{x}),transforms(x),"BoundsStyle","FollowOutput");
    outviews = arrayfun(outView, 1:n_imgs)';
    
    finalView = getLargestImgRef(outviews);
    warpedImg = cell(size(rawImages));
    for i = 1:n_imgs
        warpedImg{i} = imwarp(rawImages{i},transforms(i), "OutputView", finalView, "interp", "nearest");
    end
end

function masks = fetchMasks(filepaths)
    arguments (Input)
            filepaths (:, 1) string
        end
        arguments (Output)
            masks (:, 1) cell
        end
        fetchStr = 'Mask';
           
        nFiles = length(filepaths);    
        masks = cell(size(filepaths));
        for i = 1:nFiles
            loadStruct = load(filepaths(i), fetchStr);
            masks{i} = loadStruct.(fetchStr);
        end
end

function showImages(images)
    colors = flipud(cmapL([0 0 1; 0 1 1; 0 1 0; 1 0.7 0; 1 0 0; 0.7 0 1], length(images)));
    finalImg  = CreateRGB2(images, colors);
    imagesc(finalImg);
end

masks = fetchMasks(tester.filepath);

t1 = warpImages(registration_data.tForm, masks);

showImages(t1)

%%

[linkmat, UnifiedMasks] = createUnifiedMask(tester.filepath, tester.transform);
%%
[n_Rois, nVars] = size(linkmat);

is_valid = true; % flag to track validity
for j = 1:nVars  % loop over each mask (column)
    col_vals = linkmat(:, j);
    col_vals(col_vals == 0) = []; % ignore zeros
    if numel(unique(col_vals)) ~= numel(col_vals)
        fprintf('Duplicate ROI ID found in column %d\n', j);
        is_valid = false;
    end
end

%%
clean_unified_masks = warpImages(tester.transform, UnifiedMasks);