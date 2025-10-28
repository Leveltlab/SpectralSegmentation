function imgTForms = registerImgs(imgs, pivot, metric, optimizer, options)
    % 
    % 
    % Augustijn Vrolijk
    % 2025-7
    
    arguments (Input)
        imgs (:, 1) cell
        pivot (1, 1) double = []
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
    
    if isempty(pivot)
        refs   = 1:(nImages-1);
        toMove = 2:nImages;
    else
        refs   = pivot*ones(nImages-1, 1);
        toMove = 1:nImages;
        toMove(pivot) = [];
    end
    nToMove = length(toMove);
    
    imgTForms(nImages, 1) = affinetform2d(); %one extra to include the fixed image, it will have a blank transform
    
    for i=1:nToMove
        [fixedImg, movingImg] = warpPair(imgs{refs(i)}, imgs{toMove(i)}, imgTForms(refs(i)));
        
        progress = struct('iter',i,'total',nImages,'progressBar',options.progressBar);
        
        newTForm = PyramidRegisterPair(fixedImg, movingImg, metric, optimizer, ...
                                        progress,...
                                        levels=options.levels,...
                                        transform=options.transform, ...
                                        crossCTransform=options.crossCTransform);
        
        imgTForms(toMove(i)) = affinetform2d(imgTForms(toMove(i)).A * newTForm.A);
    end
    
    if ~progExist
        close(options.progressBar)
    end
    
    function [fixedImg, movingImg] = warpPair(fixed, moving, tForm, isPivot)
        noChange = affinetform2d();
        fixedView = affineOutputView(size(fixed),noChange,"BoundsStyle","FollowOutput");
        movingView = affineOutputView(size(moving),tForm,"BoundsStyle","FollowOutput");
        finalOutputView = getLargestImgRef([fixedView, movingView]);
        if isPivot
            fixedImg = imwarp(fixed, noChange, "OutputView", finalOutputView);
        else
            fixedImg = imwarp(fixed, tForm, "OutputView", finalOutputView);
        end
        movingImg = imwarp(moving, tForm, "OutputView", finalOutputView);     
    end
end