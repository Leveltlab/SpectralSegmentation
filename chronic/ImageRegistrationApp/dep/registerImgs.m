function imgTForms = registerImgs(imgs, pivot, metric, optimizer, options)
    % 
    % 
    % Augustijn Vrolijk
    % 2025-7
    
    arguments (Input)
        imgs (:, 1) cell
        pivot double = []
        metric = []
        optimizer = []
        options.levels double = 4
        options.transform string = "affine"
        options.crossCTransform string = "translation"
        options.progressBar = 0;
    end
    arguments (Output)
        imgTForms (:, 1) affinetform2d
    end
    
    progExist = true;
    if isequal(options.progressBar, 0)
        disp("we are here and we shouldn't")
        progExist = false;
        options.progressBar = waitbar(0, "Initialising ...");
    end
    
    [defaultOptimizer, defaultMetric] = imregconfig("monomodal");
    if isempty(optimizer)
        optimizer = defaultOptimizer;
    end
    if isempty(metric)
        metric = defaultMetric;
    end
    
    nImages = length(imgs);
    % Get which recording to use as reference and which to move for every image
    if isempty(pivot) % Serially registering
        refs   = 1:(nImages-1);
        toMove = 2:nImages;
    else % Registering against fixed reference image
        refs   = pivot*ones(nImages-1, 1);
        toMove = 1:nImages;
        toMove(pivot) = [];
    end
    nToMove = length(toMove);
    
    % Remove NaNs because the imregcorr cannot deal with NaN :(
    for i = 1:nImages
        imgs{i}(isnan(imgs{i})) = min(imgs{i}(:));
    end

    imgTForms(nImages, 1) = affinetform2d(); %one extra to include the fixed image, it will have a blank transform
    
    for i=1:nToMove
        % Using the transform of the reference to build upon it
        fixedView = affineOutputView(size(imgs{refs(i)}), imgTForms(refs(i)),"BoundsStyle","FollowOutput");
        movingView = affineOutputView(size(imgs{toMove(i)}), imgTForms(refs(i)),"BoundsStyle","FollowOutput");
        outView = getLargestImgRef([fixedView, movingView]);
        fixedImg = imwarp(imgs{refs(i)}, imgTForms(refs(i)), "OutputView", outView);
        movingImg = imwarp(imgs{toMove(i)}, imgTForms(refs(i)), "OutputView", outView);
        
        % figure('OuterPosition',[248 42 1298 954])
        % x = outView.XWorldLimits;
        % y = outView.YWorldLimits;
        % imagesc(x, y, CreateRGB({fixedImg, movingImg}, 'gb r'))
        % xlim([-100 800])
        % ylim([-100 600])
        % title(sprintf('%d pre regi: rec %d - %d', i, refs(i), toMove(i)))
        
        progress = struct('iter',i,'total',nToMove,'progressBar',options.progressBar);
        
        newTForm = PyramidRegisterPair(fixedImg, movingImg, metric, optimizer, ...
                                        progress,...
                                        levels=options.levels,...
                                        transform=options.transform, ...
                                        crossCTransform=options.crossCTransform);
        
        imgTForms(toMove(i)) = affinetform2d(newTForm.A * imgTForms(refs(i)).A);

        % fixedView = affineOutputView(size(imgs{refs(i)}), imgTForms(refs(i)),"BoundsStyle","FollowOutput");
        % movingView = affineOutputView(size(imgs{toMove(i)}), imgTForms(toMove(i)),"BoundsStyle","FollowOutput");
        % outView = getLargestImgRef([fixedView, movingView]);
        % fixedImg = imwarp(imgs{refs(i)}, imgTForms(refs(i)), "OutputView", outView);
        % movingImg = imwarp(imgs{toMove(i)}, imgTForms(toMove(i)), "OutputView", outView);
        % figure('OuterPosition',[248 42 1298 954])
        % x = outView.XWorldLimits;
        % y = outView.YWorldLimits;
        % imagesc(x, y, CreateRGB({fixedImg, movingImg}, 'gb r'))
        % xlim([-100 800])
        % ylim([-100 600])
        % title(sprintf('%d post regi: rec %d - %d: correl%.3f', i, refs(i), toMove(i), corr(fixedImg(:), movingImg(:))))

    end
    
    if ~progExist
        close(options.progressBar)
    end
end