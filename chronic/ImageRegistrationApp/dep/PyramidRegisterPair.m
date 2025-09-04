function finalTForm = PyramidRegisterPair(fixed, moving, metric, optimizer, progress, options)
    % 
    % Augustijn Vrolijk
    % 2025-4
    % 
    arguments (Input)
        fixed (:, :) double 
        moving (:, :) double
        metric = []
        optimizer = []
        progress = struct();
        options.levels double = 4
        options.transform string = "affine"
        options.crossCTransform string = "translation"
    end
    arguments (Output)
        finalTForm affinetform2d
    end
    
    defaultProg = struct('iter',1,'total',1,'progressBar',0);
    progress = mergeStructs(defaultProg, progress);
    %LOOK into a default check function, for all the optimise funcs, which
    %does the optimizer/metric/levels/transform checking
    
    [defaultOptimizer, defaultMetric] = imregconfig("monomodal");
    if isempty(optimizer)
        optimizer = defaultOptimizer;
    end
    if isempty(metric)
        metric = defaultMetric;
    end
    if options.levels < 1
        error("Must select a level strength greater than or equal to 1")
    end    
    
    pyramidImgTable = table('Size', [options.levels, 2], 'VariableTypes',{'cell', 'cell'}, 'VariableNames',{'fixed', 'moving'});
    [pyramidImgTable.fixed{1}, pyramidImgTable.moving{1}] = trimToSize(fixed, moving);
    
    for i=2:options.levels
        pyramidImgTable.fixed{i} = impyramid(pyramidImgTable.fixed{i-1}, "reduce");
        pyramidImgTable.moving{i} = impyramid(pyramidImgTable.moving{i-1}, "reduce");
    end
    %can change the transform type for cross corr, but has a risk of
    %unpredictable behaviour
    if ~isequal(progress.progressBar, 0)
        level = 1;
        updateWaitbar(level);
    end  
    curTForm = imregcorr(pyramidImgTable.moving{options.levels},pyramidImgTable.fixed{options.levels},options.crossCTransform);
    
    for i = options.levels-1:-1:1
        if ~isequal(progress.progressBar, 0)
            level = options.levels+1 - i;
            updateWaitbar(level);
        end    
        nextTForm = pyramidUpTform(curTForm); %scale up transform for the next pyramid layer
        curTForm = imregtform(pyramidImgTable.moving{i}, pyramidImgTable.fixed{i}, options.transform, optimizer, metric, PyramidLevels=1, InitialTransformation=nextTForm);
    end
    
    finalTForm = curTForm;
    
    
    function tformUp = pyramidUpTform(inputTForm)
        %{
        I BELIEVE EVERYTHING SCALES OTHER THAN THE TRANSLATION, if its a 3x3
        matrix, the one in position, 1,3 and 2,3 are the translation vectors,
        therefore just need to double the .A (affine transform) in positions
        1,3 and 2,3 to get the upscaled transform
        %}
        transformMatrix = inputTForm.A;
        transformMatrix([1,2], 3) = 2*transformMatrix([1,2], 3);
        tformUp = affinetform2d(transformMatrix);
    end
    
    function [fixedTrimmed, movingTrimmed] = trimToSize(fixed, moving)
        fixedSize = size(fixed);
        movingSize = size(moving);
        minrow = min(fixedSize(1), movingSize(1)); 
        mincol = min(fixedSize(2), movingSize(2)); 
    
        fixedTrimmed = fixed(1:minrow, 1:mincol);
        movingTrimmed = moving(1:minrow, 1:mincol);
    end
    
    function updateWaitbar(tLevel)
        progressPercent = ((progress.iter - 1) * options.levels + tLevel) / (progress.total * options.levels);
        message = sprintf('Registering image pair: %d/%d on pyramid level: %d/%d', ...
            progress.iter, progress.total, tLevel, options.levels);
        waitbar(progressPercent, progress.progressBar, message);
    end
end