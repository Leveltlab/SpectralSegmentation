function outputTable = fetchData(fileNames, filePaths, vars, nFiles)
    arguments (Input)
        fileNames (:,1) cell
        filePaths (:,1) cell
        vars (1, :) cell = {'Mask', 'PP', 'BImg'}
        nFiles double = length(filenames)
    end
    arguments (Output)
        outputTable (:, :) table
    end

    for i = 1:nFiles
        fileName = fullfile(filePaths{i}, fileNames{i});
        temp(i) = load(fileName, vars{:});
       
    end
    outputTable = struct2table(temp);
end

function fileNames = mergeNames(filenames, filepaths, filedates)
 arguments (Input)
        filenames (:,1) cell
        filepaths (:,1) cell
        filedates (:,1) datetime
    end
    arguments (Output)
        fileNames (:, 4) table
    end
    filenamesShort = ShortenFileNames(filenames, 1, filedates);
    fileDates = cellstr(filedates);
    fileNames = table(filenames, filepaths, fileDates, filenamesShort);
end

function outputTable = cleanMultImg(Images, nFiles)
    arguments (Input)
        Images (:, 1) cell 
        nFiles double = length(Images)
    end
    arguments (Output)
        outputTable (:, 3) table
    end
    
    varTypes = {'double', 'double', 'cell'};
    colNames = {'minVal', 'maxVal', 'cleanImg'};
    outputTable = table('Size', [nFiles, length(varTypes)], 'VariableTypes',varTypes, 'VariableNames',colNames);

    for i = 1:nFiles
        [outputTable.cleanImg{i}, outputTable.minVal(i), outputTable.maxVal(i)] = cleanImage(Images{i});
    end
end

function outputTable = binarifyMultMask(Masks, nFiles)
    arguments (Input)
        Masks (:, 1) cell 
        nFiles double = length(Masks)
    end
    arguments (Output)
        outputTable (:, 1) table
    end
    varTypes = {'cell'};
    colNames = {'binaryMask'};
    outputTable = table('Size', [nFiles, length(varTypes)],'VariableTypes',varTypes, 'VariableNames',colNames);

    for i = 1:nFiles
        outputTable.binaryMask{i} = binaryMaskify(Masks{i});
    end
end

function binMask = binaryMaskify(mask)
    arguments (Input)
        mask (:, :) double
    end
    arguments (Output)
        binMask (:, :) double
    end

    binMask = mask;
    binMask(mask > 0) = 1;
end

function [cleanImg, minVal, maxVal] = cleanImage(Img)
    arguments (Input)
        Img (:,:) double
    end
    arguments (Output)
        cleanImg (:,:) double
        minVal double
        maxVal double
    end

    isClean = false;
    cleanImg = Img;
    
    while not (isClean) %trim outer image border

        isClean = true;
        if isequal(cleanImg(1,:),cleanImg(end,:))
            cleanImg = cleanImg(2:end-1,:);
            isClean = false;
        end
        if isequal(cleanImg(:,1),cleanImg(:,end))
            cleanImg = cleanImg(:,2:end-1);
            isClean = false;
        end
    end
    infMask = isinf(cleanImg); % create mask, 1 if value at index is inf, 0 if not
    [minVal, maxVal] = bounds(cleanImg(~infMask)); % check the bounds of the values not at inf
    cleanImg(infMask) = minVal; % clean inf values with true min vals
    cleanImg = (cleanImg - minVal) ./ range([minVal, maxVal]); % normalise image to have range [0, .., 1]
end

function [filenames, filepaths] = selectFolder()
    arguments (Output)
        filenames (:, :) cell
        filepaths (:, :) cell
    end
    fileInfo = dir(uigetdir());
    fileInfo([fileInfo.isdir]) = [];

    filenames = {fileInfo.name}';
    filepaths = {fileInfo.folder}';
end

function fileInfo = FetchFilesAndDates()
    arguments (Output)
        fileInfo table
    end
    
    [filenames, filepaths] = selectFolder();
    nFiles = length(filenames);
    filedates = NaT(nFiles, 1);

    for i=1:nFiles
        filedates(i) = FindDate(filepaths{i}, filenames{i});
    end
    
    filenamesShort = ShortenFileNames(filenames, 1, filedates);
    fileInfo = table(filenames, filepaths, filedates, filenamesShort);
    fileInfo = sortrows(fileInfo, 'filedates');   
end

function data = main()
    fileInfo = FetchFilesAndDates(); 
    nfiles = height(fileInfo);

    data = fetchData(fileInfo.filenames, fileInfo.filepaths, {'Mask', 'PP', 'BImg'}, nfiles);
    data = [data, cleanMultImg(data.BImg), binarifyMultMask(data.Mask)];
    
    optimiseFolder(data.binaryMask{1}, data.binaryMask(2:end), fileInfo.filenamesShort, originalImgs=data.cleanImg);
    optimiseFolder(data.cleanImg{1}, data.cleanImg(2:end), fileInfo.filenamesShort, originalImgs=data.binaryMask);

end

function finalImgRef = getLargestImgRef(ImgRefs)
    arguments (Input)
        ImgRefs (:, 1) imref2d
    end
    arguments (Output)
        finalImgRef imref2d
    end
    %first retrieve the world limit property into a 2*n matirx, then
    %reshape to an array of 1*2n. 
    XWorldLimitsArray = reshape(vertcat(ImgRefs.XWorldLimits), 1, []);
    YWorldLimitsArray = reshape(vertcat(ImgRefs.YWorldLimits), 1, []);
    %retrieve the range in coords to find the smallest and largest coord
    %needed to fit everything, round it to remove small fraction errors
    [xWorldmin, xWorldmax] = bounds(round(XWorldLimitsArray));
    [yWorldmin, yWorldmax] = bounds(round(YWorldLimitsArray));
    %calculate size
    xSize = xWorldmax - xWorldmin;
    ySize = yWorldmax - yWorldmin;
    finalImgRef = imref2d([ySize xSize], [xWorldmin, xWorldmax], [yWorldmin, yWorldmax]);
    if finalImgRef.PixelExtentInWorldX ~= 1 || finalImgRef.PixelExtentInWorldY ~= 1
        error('the calculated Pixel extent is wrong, this is a bug');
    end
end

function tformUp = pyramidUpTform(inputTForm)
    arguments (Input)
    inputTForm affinetform2d
    end
    arguments (Output)
    tformUp affinetform2d
    end
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
    arguments (Input)
        fixed (:, :) double 
        moving (:, :) double     
    end
    arguments (Output)
        fixedTrimmed (:, :) double 
        movingTrimmed (:, :) double   
    end
    fixedSize = size(fixed);
    movingSize = size(moving);
    minrow = min(fixedSize(1), movingSize(1)); 
    mincol = min(fixedSize(2), movingSize(2)); 

    fixedTrimmed = fixed(1:minrow, 1:mincol);
    movingTrimmed = moving(1:minrow, 1:mincol);

end

function [fixedWarp, movingWarp] = warpPair(fixed, moving, tForm)
    arguments (Input)
        fixed (:, :) double
        moving (:, :) double
        tForm affinetform2d
    end
    arguments (Output)
        fixedWarp (:, :) double
        movingWarp (: , :) double
    end
    noChange = affinetform2d();
    fixedView = affineOutputView(size(fixed),noChange,"BoundsStyle","FollowOutput");
    movingView = affineOutputView(size(moving),tForm,"BoundsStyle","FollowOutput");
    finalOutputView = getLargestImgRef([fixedView, movingView]);
    fixedWarp = imwarp(fixed, noChange, "OutputView", finalOutputView);
    movingWarp = imwarp(moving, tForm, "OutputView", finalOutputView);     
end

function finalTForm = compactOptimisePyramids(fixed, moving, metric, optimizer, options)
    arguments (Input)
        fixed (:, :) double 
        moving (:, :) double
        metric = []
        optimizer = []
        options.levels double = 4
        options.transform string = "affine"
    end
    arguments (Output)
        finalTForm affinetform2d
    end
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

    curTForm = imregcorr(pyramidImgTable.moving{options.levels},pyramidImgTable.fixed{options.levels}, "translation");

    for i = options.levels-1:-1:1
        nextTForm = pyramidUpTform(curTForm); %scale up transform for the next pyramid layer
        curTForm = imregtform(pyramidImgTable.moving{i}, pyramidImgTable.fixed{i}, options.transform, optimizer, metric, PyramidLevels=1, InitialTransformation=nextTForm);
    end

    finalTForm = curTForm;
end

function optimiseFolder(fixed, Images, imgNames, options)
    arguments (Input)
        fixed (:, :) double
        Images (:, :) cell
        imgNames (:, :) cell = arrayfun(@num2str, 1:length(imgs), 'UniformOutput', false);
        options.originalImgs = Images;
        options.nFiles double = length(Images)
    end
    [optimizer, metric] = imregconfig("monomodal");
    columnNames =  {'fixed', 'moving', 'fromPrevTForm', 'netTForm'};
    varTypes = {'cell', 'cell', 'affinetform2d','affinetform2d'};
    rowNames =  arrayfun(@(i) sprintf('transform from image: (%d -> fixed)', i), 1:options.nFiles, 'UniformOutput', false);

    outputTable = table('Size', [options.nFiles, length(columnNames)], 'VariableTypes',varTypes, 'VariableNames',columnNames, 'RowNames',rowNames);
    outputTable.moving = Images;
    outputTable.fixed = repmat({fixed}, options.nFiles, 1);

    nextTForm = affinetform2d();    
    noChange = nextTForm;

    for i=1:options.nFiles
        moving = Images{i};
        
        %plotImagesVarNames(fixed, moving);
    
        [fixedWarp, movingWarp] = warpPair(fixed, moving, nextTForm);

        %plotImagesVarNames(fixedWarp, movingWarp);

        newTForm = compactOptimisePyramids(fixedWarp, movingWarp, metric, optimizer,levels=4,transform="affine");
        
        outputTable.fromPrevTForm(i) = newTForm;
        nextTForm = affinetform2d(nextTForm.A*newTForm.A);
        
        %plotPairWithTForm(fixed, moving, nextTForm);
        disp(i);
        outputTable.netTForm(i) = nextTForm;
    end
    plotImagestFormBulk([{fixed}, outputTable.moving(:)'], [noChange, outputTable.netTForm(:)'], imgNames);
    plotImagestFormBulk(options.originalImgs, [noChange, outputTable.netTForm(:)'], imgNames);
end

function plotDisplaceField(moving, fixed)
    subsampleRatio = 10;
    [D, registeredMoving] = imregdemons(moving, fixed);
    plotImagesGetNames(registeredMoving, fixed);

    Ux = -D(:,:,1);
    Uy = -D(:,:,2);

    Ux_sub = Ux(1:subsampleRatio:end, 1:subsampleRatio:end);
    Uy_sub = Uy(1:subsampleRatio:end, 1:subsampleRatio:end);
    
    % Generate the grid for the displacement field
    [xGrid, yGrid] = meshgrid(1:size(Ux, 2), 1:size(Ux, 1));
    xGrid_sub = xGrid(1:subsampleRatio:end, 1:subsampleRatio:end);
    yGrid_sub = yGrid(1:subsampleRatio:end, 1:subsampleRatio:end);

    quiver(xGrid_sub, yGrid_sub, Ux_sub, Uy_sub);
    
end

function plotImagesVarNames(name)
    arguments (Input, Repeating)
        name 
    end
    len = nargin;
    names = cell(1, len);
    cellArr = cell(1, len);

    for i = 1:len
        names{i} = inputname(i);
        cellArr{i} = name{i};
    end
    plotImages(cellArr, names);

end

function plotImagestFormBulk(imgs, tForms, imgNames, nFiles)
    arguments (Input)
        imgs (:, :) cell
        tForms (:, :) affinetform2d
        imgNames (:, :) cell = arrayfun(@num2str, 1:length(imgs), 'UniformOutput', false);
        nFiles double = length(imgs)
    end

    if length(tForms) ~= length(imgs)
        error("different length images and tform arrays for plotImagestFormBulk");
    end
    
    outputViews = cell(nFiles, 1);
    warpedImgs = cell(nFiles, 1);

    for i = 1:nFiles
        outputViews{i} = affineOutputView(size(imgs{i}),tForms(i),"BoundsStyle","FollowOutput");
    end
    finalOutputView = getLargestImgRef([outputViews{:}]);
    
    for i = 1:nFiles
        warpedImgs{i} = imwarp(imgs{i}, tForms(i), "OutputView", finalOutputView);
    end
    
    plotImages(warpedImgs,imgNames);
end

function plotPairWithTForm(fixed, moving, tform)
    arguments (Input)
        fixed (:, :) double
        moving (:, :) double
        tform affinetform2d
    end
    fixedView = affineOutputView(size(fixed),affinetform2d(),"BoundsStyle","FollowOutput");
    movingView = affineOutputView(size(moving),tform,"BoundsStyle","FollowOutput");
    finalOutputView = getLargestImgRef([fixedView, movingView]);
    fixedWarp = imwarp(fixed, affinetform2d(), "OutputView", finalOutputView);
    movingWarp = imwarp(moving, tform, "OutputView", finalOutputView);
    plotImagesVarNames(fixedWarp, movingWarp);
end

function plotImgDistribution(img, accuracy)
    arguments (Input)
        img (:, :) double
        accuracy double = 2
    end
    roundedVals = round(img, accuracy);
    flattenedValues = roundedVals(:);
    [uniqueVals,~,idx] = unique(flattenedValues);
    counts = accumarray(idx, 1);
    
    figure;
    bar(uniqueVals, counts, 'FaceColor',[0.2 0.6 0.8]);
    xlabel('Value (Rounded to 2 Decimals)');
    ylabel('Frequency');
    title('Histogram of Rounded Values');
    grid on;

end

function plotImages(images, fileNames)
    arguments (Input)
        images (:, 1) cell
        fileNames (:, :) cell = cell(1, length(images));
    end
    nFiles = length(images);
    toplot = 1:nFiles; % sessions to plot
    
    colors = flipud(cmapL([0 0 1; 0 1 1; 0 1 0; 1 0.7 0; 1 0 0; 0.7 0 1], nFiles));
    if size(colors,1)==3
        colors = [1 0 0; 0 1 0; 0.3 0.3 1];
    elseif size(colors,1)==2
        colors = [1 0 0; 0 1 0.5];
    end

    % Plot image of unregistered overlayed recordings
    figure('Units','Normalized', 'Position', [0.1766 0.1407 0.4203 0.7074])
    subplot('Position', [0 0 1 1])
    RGB = CreateRGB2(images(toplot), colors(toplot,:));
    imagesc(RGB)

    for i = 1:length(toplot)
        text(20, 20+i*17,fileNames(toplot(i)),'color',colors(toplot(i),:),...
            'fontweight','bold','fontsize',12)
    end
    hold on
    %{
    ------------ 
    NEED TO ADD FUNCTIONS TO PLOT PPs etc... maybe do
    seperately? pass plot object and contain PPs plotting there
    ------------
    for r = 1:nfiles
        for i = 1:PPs(r).Cnt
            plot(PPs(r).Con(i).x, PPs(r).Con(i).y, 'color', [colors(r,:) 0.2])
        end
    plot(PPs(r).P(1, :), PPs(r).P(2, :), '.', 'color', colors(r,:))
    end
    %}
    figtitle('non-registered images', 'Color', 'w')
end

function finalTForm = compactOptimise(fixed, moving, options)
    arguments (Input)     
        fixed (:, :) double 
        moving (:, :) double
        options.optimizer = []
        options.metric = []
        options.initialTform affinetform2d = affinetform2d()
        options.transform string = 'affine'
        options.levels double = 1
    end
    arguments (Output)
        finalTForm affinetform2d
    end

    [optimizer, metric] = imregconfig("monomodal");
    if isempty(options.optimizer)
        options.optimizer = optimizer;
    end
    if isempty(options.metric)
        options.metric = metric;
    end

    validTForms = {'translation', 'rigid', 'similarity', 'affine'};
    isValid = ismember(options.transform, validTForms);
    if ~ isValid
        error("the specified transform for compactOptimise is not valid")
    end
    if options.levels < 1
        error("the specified levels amount for compactOptimise must be greater than or equal to 1")
    end 
    if options.levels > 6
        error("the specified levels amount for compactOptimise must be less than 6, (can modify this requirement but not recommended)")
    end

    finalTForm = imregtform(moving, fixed, options.transform, options.optimizer, options.metric, PyramidLevels=options.levels, InitialTransformation=options.initialTform);
    
end

function [outputTable, finalTForm] = optimise(optimizer, metric, transforms, fixed, moving, initialTform)
    arguments (Input)
        optimizer
        metric
        transforms (:,1) string
        fixed (:, :) double 
        moving (:, :) double
        initialTform affinetform2d = affinetform2d()
    end
    arguments (Output)
        outputTable table
        finalTForm affinetform2d
    end
    columnNames =  {'tForm', 'outputView', 'warpedImg'};
    rowNames =  [{'fixed'},{'moving'},transforms(:)'];
    nTransforms = length(transforms)+2; %first two rows are the fixed and moving with the new output view
    outputTable = table('Size', [nTransforms, 3], 'VariableTypes',{'affinetform2d', 'imref2d','cell'}, 'VariableNames',columnNames, 'RowNames',rowNames);
    

    dynamicTforms = cellfun(@(i) imregtform(moving, fixed, i, optimizer, metric, PyramidLevels=1, InitialTransformation=initialTform),transforms, 'UniformOutput', false);
    noChange = affinetform2d();
    allTforms = [{noChange}, {noChange}, dynamicTforms(:)'];
    outputTable.tForm = allTforms';

    outputTable.outputView(1) = affineOutputView(size(fixed),outputTable.tForm{1},"BoundsStyle","FollowOutput");
    for i=2:nTransforms
        outputTable.outputView(i) = affineOutputView(size(moving),outputTable.tForm{i},"BoundsStyle","FollowOutput");
    end

    finalOutputView = getLargestImgRef(outputTable.outputView(:));
    outputTable.warpedImg{1} =  imwarp(fixed, outputTable.tForm{1} , "OutputView", finalOutputView);

    for i=2:nTransforms
        outputTable.warpedImg{i} =  imwarp(moving, outputTable.tForm{i} , "OutputView", finalOutputView);
    end
    
    finalTForm = outputTable.tForm{transforms{end}};
end

function [outputTable, finalTForm] = optimisePyramids(metric, optimizer, fixed, moving, options)
    arguments (Input)
        metric
        optimizer
        fixed (:, :) double 
        moving (:, :) double
        options.levels double = 4
        options.transform string = "affine"
    end
    arguments (Output)
        outputTable table
        finalTForm affinetform2d
    end

    if options.levels < 1
        error("Must select a level strength greater than or equal to 1")
    end
    
    columnNames =  {'fixed', 'moving', 'tForm', 'ogtForm'};
    varTypes = {'cell', 'cell', 'affinetform2d', 'affinetform2d'};
    rowNames =  arrayfun(@(i) sprintf('level %d', i), 1:options.levels, 'UniformOutput', false);
    outputTable = table('Size', [options.levels, length(columnNames)], 'VariableTypes',varTypes, 'VariableNames',columnNames, 'RowNames',rowNames);
    [fixedTrimmed, movingTrimmed] = trimToSize(fixed, moving);
    outputTable.fixed{1} = fixedTrimmed;
    outputTable.moving{1} = movingTrimmed;
        
    for i=2:options.levels
        outputTable.fixed{i} = impyramid(outputTable.fixed{i-1}, "reduce");
        outputTable.moving{i} = impyramid(outputTable.moving{i-1}, "reduce");
    end
    %plotImagesVarNames(outputTable.fixed{options.levels}, outputTable.moving{options.levels})
    [outputTable.tForm(options.levels), peakCor] = imregcorr(outputTable.moving{options.levels},outputTable.fixed{options.levels}, "translation");
    %plotPairWithTForm(outputTable.fixed{options.levels}, outputTable.moving{options.levels}, outputTable.tForm(options.levels));

    for i = options.levels-1:-1:1
        nexttForm = pyramidUpTform(outputTable.tForm(i+1)); %scale up transform for the next pyramid layer
        outputTForm = compactOptimise(outputTable.fixed{i}, outputTable.moving{i}, optimizer=optimizer, metric=metric, initialTform=nexttForm, transform=options.transform, levels=1);
        outputTable.tForm(i) = outputTForm;
        %plotPairWithTForm(outputTable.fixed{i}, outputTable.moving{i}, outputTable.tForm(i))
    end

    %plotPairWithTForm(outputTable.fixed{1}, outputTable.moving{1}, outputTable.tForm(1))
    finalTForm = outputTable.tForm(1);
end

data = main();

%{
          --------------------- TODO: ----------------------
In optimisePyramids, the images are trimmed to be the same size
I then return the final transform which is applied, but this warping may
shift the image around, which then gets excacerbated when trimmed again,
needd to think about how to handle this.

%}