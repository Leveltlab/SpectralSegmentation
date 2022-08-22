function ViewVolumeRGB1(img, scalexy, z)
% ViewVolumeRGB1({mat1 mat2}, scalexc, z) visualizes the 3D matrixes
%
% First input can either be one 3D matrix, or a cell array filled with
% 3D matrices.
% Second input: the scale of the pixels in x y direction in um. (a double)
% Third input: the axis of the z direction in um. (vector)
% 
% The first figure will show the max projection through the entire z-stacks
% with the reference lines of where the zstack will be cut. 
% The data in these cuts is visualized in second figure by max projection.
% 
% Go through the zstacks by clicking on the reference image, or by using
% the keyboard: 'a'(left), 'd'(right), 
%               's'(down/to back/ y direction), 'w'(up/ to front), 
%               'r'(shallower/ z direction), 'f'(deeper),
%               'g'(decrease cut size), 't'(increase cut size)
%               
% different colors options by different amount of inputs.
% Possible inputs: One 3D matrix is will be visualized in black and white
%                  Two 3D matrices will be red and green
%                  Three 3D matrices will be red green and blue
%                  Four 3D matrices will be red green blue and cyan
%                  Five 3D matrices will be red green blue cyan and magenta
%         other amount of matrices in the cell array input will be refused
% 
%  
% 
% Leander de Kraker
% 2019-11-18
%

    z = double(z);
    cutSize = 5; % width for maxprojection slice in idx
    cutSizeOld = 0; % Track if the cutsize changes
    
    % 3D matrix should be in a cell
    if ~iscell(img)
        img = {img};
    end
    
    nimgs = length(img); % amount of z-stacks
    % which colors to use:
    if nimgs == 1
        colors = 'rgb';
    elseif nimgs == 2
        colors = 'r g';
    elseif nimgs == 3
        colors = 'r g b';
    elseif nimgs == 4
        colors = 'r g b gb';
    elseif nimgs == 5
        colors = 'r g b gb rg';
    else
        warning('unrecognized/unsupported number of data arrays was given')
        return
    end
    
    [heig, wid, ~]= size(img{1});
    x = (1:wid)*scalexy; % The position of pixels in um
    y = (1:heig)*scalexy;

    % Check z axis
    if length(z) ~= size(img{1},3)
        warning('z axis length is not equal to data z lenght!!!!!')
        if length(z) > size(img{1},3)
            warning('adjusting the z axis length to correspond to data')
            z = z(1:size(img{1},3));
        else
            warning('make z axis length the same as 3th dimension of the data!')
            return
        end
    end
        
    
    m = 0.05; % margin to the sides of the figure (relative)
    ys = double(max(y)); % dimensions of the stack in um
    xs = double(max(x));
    zs = double(max(z));
    h = (ys + zs);
    w = (xs + zs); % Total width of the image when they are plotted next to each other

    % position = [horstart, vertstart, width, height] (start at bottomleft)
    position1 = [m,    ys/h, xs/w-m, zs/h-m]; % The first vertical slice on top of the x direction
    position2 = [m,    m,    xs/w-m, ys/h-m]; % The image as we usually see it.
    position3 = [xs/w, m,    zs/w-m, ys/h-m]; % The vertical slice 90 degrees rot of 1st one
    axposition = [position1; position2; position3];
    
    % The figure with the example data and slice selector
    % With this figure everything that is plotted has a handle that gets
    % kept up to date by deleting the graphic object when it has to update
    hF1 = figure(1);
    set(hF1,'units','normalized', 'OuterPosition',[0 0.05 0.4 0.7], 'KeyPressFcn',@KeyPresser)
    h1xz = subplot('position', axposition(1,:),...
                   'ButtonDownFcn', {@axesXZfcn}, 'nextplot', 'add');
    h1xy = subplot('position', axposition(2,:),...
                   'ButtonDownFcn', {@axesXYfcn}, 'nextplot', 'add');
    h1zy = subplot('position', axposition(3,:),...
                   'ButtonDownFcn', {@axesYZfcn}, 'nextplot', 'add');
               
    strExplain = {'w s: move up,   down';...
                  'a d: move left, right';...
                  'r  f: move z up, z down';...
                  't: enlarge cut';... 
                  'g: smaller cut'};
    annotation('textbox', [xs/w, 1-m, zs/w, 0], 'string', strExplain)
               
    % The figure that has the sliced data
    % With this figure only slices are plotted, those are replaced every
    % update
    hF2 = figure(2);
    set(hF2,'units','normalized', 'OuterPosition',[0.5 0.05 0.4 0.7], 'KeyPressFcn',@KeyPresser)
    h2xz = subplot('position', axposition(1,:),...
                   'ButtonDownFcn', {@axesXZfcn}, 'nextplot', 'replaceChildren');
    h2xy = subplot('position', axposition(2,:),...
                   'ButtonDownFcn', {@axesXYfcn}, 'nextplot', 'replaceChildren');
    h2zy = subplot('position', axposition(3,:),...
                   'ButtonDownFcn', {@axesYZfcn}, 'nextplot', 'replaceChildren');
    
    % Default slice is the middle of the data
    % Cut is calculated based on mu
    xCut = x(end)/2;
    yCut = y(end)/2;
    zCut = max(z)/2;
    
    % Initialize the x,y,zCutIdx variable
    % This cut is used to select data, because it has actual indexes
    xCutIdx = 0;
    yCutIdx = 0;
    zCutIdx = 0;
    
    % Keep track of changes
    xCutIdxOld = xCutIdx; 
    yCutIdxOld = yCutIdx;
    zCutIdxOld = zCutIdx;
    FindCutIdx
    
    % draw the xz maxprojection (top plot)
    projections = {};
    for i = 1:nimgs
        projection = squeeze(max(img{i}, [], 1))';
        projections{end+1} = projection;
    end
    RGBXZ = CreateRGB(projections, colors);
    imagesc(x,z,RGBXZ, 'Parent', h1xz, 'hittest', 'off')
    % draw the xy maxprojection (center plot)
    projections = {};
    for i = 1:nimgs
        projection = squeeze(max(img{i}, [], 3));
        projections{end+1} = projection;
    end
    RGBXY = CreateRGB(projections, colors);
    imagesc(x,y,RGBXY, 'Parent', h1xy, 'hittest', 'off')
    % draw the zy maxprojection (right plot)
    projections = {};
    for i = 1:nimgs
        projection = squeeze(max(img{i}, [], 2));
        projections{end+1} = projection;
    end
    RGBZY = CreateRGB(projections, colors);
    imagesc(z,y,RGBZY, 'Parent', h1zy, 'hittest', 'off')
    
    % Set the correct axis limits for all subplots
    xlim([h1xz, h2xz, h1xy h2xy], x([1 end]))
    xlim([h1zy, h2zy],            z([1 end]))
    ylim([h1xy, h2xy, h1zy, h2zy], y([1 end]))
    ylim([h1xz, h2xz],             z([1 end]))
    set([h1xy h1xz h1zy], 'YDir', 'normal', 'XDir', 'normal')
    
    
    
    % Plot the cut position lines
    hxylx = plot(x([1 end]), [yCut yCut], '-w', 'Parent', h1xy, 'hittest', 'off');
    hxyly = plot([xCut xCut], y([1 end]), '-y', 'Parent', h1xy, 'hittest', 'off');
    hzylx = plot(z([1 end]), [yCut yCut], '-w', 'Parent', h1zy, 'hittest', 'off');
    hzyly = plot([zCut zCut], y([1 end]), '-y', 'Parent', h1zy, 'hittest', 'off');
    hxzlx = plot(x([1 end]), [zCut zCut], '-w', 'Parent', h1xz, 'hittest', 'off');
    hxzly = plot([xCut xCut], z([1 end]), '-y', 'Parent', h1xz, 'hittest', 'off');

    % Draw the initial slices in figure 2
    DrawData

    clearvars m xs zs ys h w

    % --- Executes on mouse press over axes background.
    function axesXYfcn(~, eventdata)
        % hObject    handle to axesXY (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % h          structure with handles and user data (see GUIDATA)
        
        xCut = round(eventdata.IntersectionPoint(1));
        yCut = round(eventdata.IntersectionPoint(2));
        
        % Make sure the click doesn't cause out of bound data selections
        if xCut + cutSize >= wid
            xCut = xCut - cutSize - 1;
        end
        if yCut + cutSize >= heig
            yCut = yCut - cutSize - 1;
        end
        
        slicer
    end

    % --- Executes on mouse press over axes background.
    function axesXZfcn(~, eventdata)
        % hObject    handle to axesXZ (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % h          structure with handles and user data (see GUIDATA)

        xCut = round(eventdata.IntersectionPoint(1));
        zCut = round(eventdata.IntersectionPoint(2));
        
        % Make sure the click doesn't cause out of bound data selections
        if xCut + cutSize >= wid
            xCut = xCut - cutSize -1;
        end
        if zCut + cutSize >= length(z)
            zCut = zCut - cutSize -1;
        end
        
        slicer
    end

    % --- Executes on mouse press over axes background.
    function axesYZfcn(~, eventdata)
        % hObject    handle to axesYZ (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % h          structure with handles and user data (see GUIDATA)
        
        zCut = round(eventdata.IntersectionPoint(1));
        yCut = round(eventdata.IntersectionPoint(2));
        
        % Make sure the click doesn't cause out of bound data selections
        if zCut + cutSize >= length(z)
            zCut = zCut - cutSize -1;
        end
        if yCut + cutSize >= heig
            yCut = yCut - cutSize -1;
        end
        
        slicer
    end

    function KeyPresser(~,event)
        % This callback is used when a keyboard press hapened while one of
        % the UIs figures is activated.
        switch event.Key
            case 'a' % to the left
                if xCutIdx > 1
                    xCutIdx = xCutIdx - 1;
                end
            case 'd' % to the right
                if xCutIdx < wid-cutSize
                    xCutIdx = xCutIdx + 1;
                end
            case 's' % to up (y)
                if yCutIdx > 1
                    yCutIdx = yCutIdx - 1;
                end
            case 'w' % to down (y)
                if yCutIdx < heig-cutSize
                    yCutIdx = yCutIdx + 1;
                end
            case 'f' % to up (z)
                if zCutIdx > 1
                    zCutIdx = zCutIdx - 1;
                end
            case 'r' % to down (z)
                if zCutIdx < length(z)-cutSize
                    zCutIdx = zCutIdx + 1;
                end
            case 't' % increase cutsize 
                cutSize = cutSize + 1;
                if xCutIdx + cutSize > wid
                    xCutIdx = xCutIdx - 1;
                end
                if yCutIdx + cutSize > heig
                    yCutIdx = yCutIdx - 1;
                end
                if zCutIdx + cutSize > length(z)
                    zCutIdx = zCutIdx - 1;
                end
            case 'g' % decrease cutsize
                if cutSize >= 1
                    cutSize = cutSize - 1;
                end
        end
        if ismember(event.Key, 'asdwrfgt')
            % If a keys that does something is pressed update slices
            FindCut
            DrawSliceRefLines
            DrawData
        end
    end

    function slicer
        % Update the slices because a new position is requested
        FindCutIdx
        DrawSliceRefLines
        DrawData
    end

    function FindCutIdx
        % Given a certain x, y, z value on the axes, which index in the
        % data does that correspond to?
        [~, xCutIdx] = min(abs(x-xCut));
        [~, yCutIdx] = min(abs(y-yCut));
        [~, zCutIdx] = min(abs(z-zCut));
    end

    function FindCut
        % Given a certain x, y, z index, return the value for the axes
        xCut = x(xCutIdx);
        yCut = y(yCutIdx);
        zCut = z(zCutIdx);
    end

    function DrawSliceRefLines
        % Update the position of the slice lines in fig1
        xCut2 = x(xCutIdx+cutSize); % until where the cut cuts
        yCut2 = y(yCutIdx+cutSize);
        zCut2 = z(zCutIdx+cutSize);
        
        hxylx.XData = x([1 end end 1]);
        hxylx.YData = [yCut yCut yCut2 yCut2];
        hxyly.XData = [xCut xCut xCut2 xCut2];
        hxyly.YData = y([1 end end 1]);
        hzylx.XData = z([1 end end 1]);
        hzylx.YData = [yCut yCut yCut2 yCut2];
        hzyly.XData = [zCut zCut zCut2 zCut2];
        hzyly.YData = y([1 end end 1]);
        hxzlx.XData = x([1 end end 1]);
        hxzlx.YData = [zCut zCut zCut2 zCut2];
        hxzly.XData = [xCut xCut xCut2 xCut2];
        hxzly.YData = z([1 end end 1]);
    end

    function UpdateOldIdx
        % Update the 'old idx' to correspond to the current idx
        xCutIdxOld = xCutIdx;
        yCutIdxOld = yCutIdx;
        zCutIdxOld = zCutIdx;
        cutSizeOld = cutSize;
    end


    function DrawData
        % Calculate and draw the data selected by the slices
        % Only update a slice if the data selection for it has changed
        
        % Calculate for the xz slice (top plot)
        if (yCutIdx ~= yCutIdxOld) || (cutSize ~= cutSizeOld)
            projections = {};
            for j = 1:nimgs
    %             projection = squeeze(img{j}(yCutIdx,:,:))';
                projection = squeeze(max(img{j}(yCutIdx:yCutIdx+cutSize,:,:),[],1))';
                projections{end+1} = projection;
            end
            RGBXZ = CreateRGB(projections, colors);
            imagesc(x,z,RGBXZ, 'Parent', h2xz, 'hittest', 'off')
            set(h2xz, 'YDir', 'normal', 'XDir', 'normal')
        end
        
        % draw the xy maxprojection (center plot)
        if (zCutIdx ~= zCutIdxOld) || (cutSize ~= cutSizeOld)
            projections = {};
            for j = 1:nimgs
    %             projection = squeeze(img{j}(:,:,zCutIdx));
                projection = squeeze(max(img{j}(:,:,zCutIdx:zCutIdx+cutSize),[],3));
                projections{end+1} = projection;
            end
            RGBXY = CreateRGB(projections, colors);
            imagesc(x,y,RGBXY, 'Parent', h2xy, 'hittest', 'off')
            set(h2xy, 'YDir', 'normal', 'XDir', 'normal')
        end
        
        % draw the zy maxprojection (right plot)
        if (xCutIdx ~= xCutIdxOld)  || (cutSize ~= cutSizeOld)
            projections = {};
            for j = 1:nimgs
    %             projection = squeeze(img{j}(:,xCutIdx,:)); % simply select one
                projection = squeeze(max(img{j}(:,xCutIdx:xCutIdx+cutSize,:), [], 2));
                projections{end+1} = projection;
            end
            RGBZY = CreateRGB(projections, colors);
            imagesc(z,y,RGBZY, 'Parent', h2zy, 'hittest', 'off')
            set(h2zy, 'YDir', 'normal', 'XDir', 'normal')
        end
        
        UpdateOldIdx
    end
end

