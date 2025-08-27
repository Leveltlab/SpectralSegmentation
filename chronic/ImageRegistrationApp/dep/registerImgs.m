function imgTForms = registerImgs(inputImages, inputFixed, metric, optimizer, options)
    % Augustijn Vrolijk
    % 2025-7
    
    arguments (Input)
        inputImages (:, 1) cell
        inputFixed (1, 1) cell = {0}
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

    [fixedCell, movingCell, isPivot] = resolveImagePairs(inputFixed, inputImages);
    nImages = length(movingCell);

    imgTForms((nImages + 1), 1) = affinetform2d(); %one extra to include the fixed image, it will have a blank transform
    
    for i=1:nImages
      
        [fixedWarp, movingWarp] = warpPair(fixedCell{i}, movingCell{i}, imgTForms(i), isPivot);

        progress = struct('iter',i,'total',nImages,'progressBar',options.progressBar);

        newTForm = PyramidRegisterPair(fixedWarp, movingWarp, metric, optimizer, ...
            progress,levels=options.levels,transform=options.transform, ...
            crossCTransform=options.crossCTransform);
        
        nextTForm = affinetform2d(imgTForms(i).A*newTForm.A);

        imgTForms(i+1) = nextTForm;
    end
    
    if ~progExist
        close(options.progressBar)
    end

    function [fixedClean, movingClean, isPivot] = resolveImagePairs(fixed, images)
        %{
        I want two arrays:
                    A:                                  B:
            FIXED : img1; img1; img1     -   FIXED : img1; img2; img3;
            MOVING: img2; img3; img4     -   MOVING: img2; img3; img4;
        
        Three scenarios:
            Fixed is not given - then we do option B
            Fixed is given - not duplicate in Images - option A
            Fixed is given - is duplicate in Images - option A
        %}
        %FLAG IS GIVEN AS THE MAIN LOOP NEEDS TO ACT DIFFERENTLY IF THE
        %PIVOT  IS ALWAYS THE SAME OR NOT, otherwise it registers to a
        %basearray of actually unregistered images.

        if length(fixed) >= 2
            error("only one image can be selected as the 'pivot'");
        elseif isequal(fixed{1}, 0)
            nPairs = length(images) - 1; 
            fixedClean = cell(nPairs, 1);
            movingClean = cell(nPairs, 1);
            for j = 1:nPairs
                fixedClean{j} = images{j};
                movingClean{j} = images{j + 1};
            end
            isPivot = false;
        else
            %THIS IS THE BUG; THE DUPLICATE IS NOT NECESSARILY Images{1} SO
            %NEED TO LOOP THROUGH AND CHECK
            toRemove = [];
            for j=1:length(images)
                if isequal(images{j}, fixed{:})
                    toRemove = [toRemove; j];
                end
            end
            for j=1:length(toRemove)
                warning("Removed image: %d; This image was the same as the pivot",toRemove(j));
                images(toRemove(j)) = [];
            end
            nPairs = length(images);
            fixedClean = cell(nPairs, 1);
            movingClean = cell(nPairs, 1);
            for j = 1:nPairs
                fixedClean{j} = fixed{:};
                movingClean{j} = images{j};
            end
            isPivot = true;
        end
    end

    function [fixedWarp, movingWarp] = warpPair(fixed, moving, tForm, isPivot)
        noChange = affinetform2d();
        fixedView = affineOutputView(size(fixed),noChange,"BoundsStyle","FollowOutput");
        movingView = affineOutputView(size(moving),tForm,"BoundsStyle","FollowOutput");
        finalOutputView = getLargestImgRef([fixedView, movingView]);
        if isPivot
            fixedWarp = imwarp(fixed, noChange, "OutputView", finalOutputView);
        else
            fixedWarp = imwarp(fixed, tForm, "OutputView", finalOutputView);
        end
        movingWarp = imwarp(moving, tForm, "OutputView", finalOutputView);     
    end
end