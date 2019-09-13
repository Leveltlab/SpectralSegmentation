function ChronicViewer(BImg2, Masks2, names, nLinksMask, linkMat, contours, score, inRoi)
% UI for going through chronically analyzed recordings
%
% The masks view only shows masks that are chronically linked to another
% ROI in at least 1 other recording.
%
% Todo:
%  http://undocumentedmatlab.com/blog/customizing-listbox-editbox-scrollbars
%  
%  recalculate the scores after the cells are edited. 
%  import the variable inRoi
% 


nfiles = length(Masks2);
dims = size(Masks2{1});

sShow = true(1, nfiles);
sCont = 1:nfiles;
rShow = 1:length(linkMat);
selCells = [2,2]; % which cells of the table are selected
selFiles = [1]; % Which files does that selection correspond to
selRoi = {1}; % Which ROIs does that selection correspond to
contAlpha = 0.4;
autoZoom = false;
oldZoom = [1, dims(1), 1, dims(2)];
viewToggle = 1;
% colors = strsplit('r b g gb rb rg r g b r'); % The colors for the sessions
% colors = colors(1:nfiles);
colors = flipud(cmapL([0 0 1; 0 1 1; 0 1 0; 1 1 0; 1 0 0; 1 0 1], nfiles));

% Activate tight subplotting
% [subplot margin side&top],[figure bottomspace,topspace],[leftspace,rightspace]
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.01], [0.025 0], [0.025 0]);

hFig = figure('name','Chronic Viewer', 'units', 'normalized', 'position', [0.1 0.05 0.8 0.85]);

hImgAx = subplot(1, 3, [1 2]); hold on
hImgAx.ButtonDownFcn = @hImgDown;
hImg = imagesc(zeros(dims([1 2])),'hittest','off');
axis([1, dims(2), 1, dims(1)])

% Draw contours and legend
hCons = {};
hLegendText = gobjects(nfiles,1);
for i = 1:nfiles
    % Draw contours and save handles
    hCons{i} = gobjects(1,contours(i).Cnt);
    for j = 1:contours(i).Cnt
        hCons{i}(j) = plot(contours(i).Con(j).y, contours(i).Con(j).x,...
            'color',[colors(i,:), contAlpha],'hittest','off');
    end
    
    % Write colored legend text
    names{i} = names{i}(1:16); % Shorten the filename
    names{i} = strrep(names{i},'_',' '); % replace underscores with spaces
    hLegendText(i) = text(3,dims(1)-12-i*12,names{i},'color',colors(i,:));
end

% Positions of the other UI elements
hRoisTAx= subplot(8, 3, [3 6 9 12 15]); % table axis
hViewAx = subplot(16, 6, 77); 
hAlphAx = subplot(16, 6, 83);
hZoomAx = subplot(16, 6, 89);
hSaveAx = subplot(16, 6, 95);
hClickTAx= subplot(8, 3, 18); % clicked on image roi information table axis
hSessTAx = subplot(4, 6, 24); % which sessions are selected table axis
hRoisTAx.Units = 'pixels';
hClickTAx.Units = 'pixels';
hViewAx.Units  = 'pixels';
hAlphAx.Units  = 'pixels';
hZoomAx.Units  = 'pixels';
hSaveAx.Units  = 'pixels';
hSessTAx.Units = 'pixels';
hViewAx.Visible  = 'off';
hRoisTAx.Visible = 'off';
hClickTAx.Visible = 'off';
hAlphAx.Visible  = 'off';
hZoomAx.Visible  = 'off';
hSaveAx.Visible  = 'off';
hSessTAx.Visible = 'off';

% Select sessions table
hSessionTable = uitable('Position',hSessTAx.Position,'Units','normalized',...
    'CellSelectionCallback',@SessionSelect);
hSessionTable.ColumnName = [{'img'},{'cont'}];
hSessionTable.ColumnWidth = [{40}, {40}];
hSessionTable.Data = repmat({true},[nfiles,2]);

% View dropdown menu
hViewPopup = uicontrol('Style', 'popup',...
    'String', {'Backgrd','Masks','nLinks','mask overlap'},...
    'Position', hViewAx.Position,'units','normalized','Callback',@BackgroundView);  

% Alpha value slider for contours
hAlph = uicontrol('Style','slider', 'Min',0, 'Max',1, 'Value',contAlpha,...
    'sliderstep',[0.1 0.1], 'Position',hAlphAx.Position,'Callback', @AlphCont);
hAlph.Position(4) = 15;
hAlphTitle = uicontrol('Style','text',...
    'position',hAlphAx.Position,'horizontalAlignment','left',...
    'string','Contours alpha');
hAlphTitle.Position([2, 4]) = [hAlphTitle.Position(2)+15, 15];
% Normalized units allows resizing of the UI
hAlph.Units   = 'normalized';
hAlphTitle.Units = 'normalized';

% automatic Zoom in functionality
hZoomButton = uicontrol('Style','togglebutton','String','auto zoom',...
    'Position',hZoomAx.Position,'Units','normalized','Callback',@AutoZoomButton);

% Delete button functionality
hSaveButton = uicontrol('Style','pushbutton','String','save table changes',...
    'Position',hSaveAx.Position,'Units','normalized','Callback',@SaveRois);
hSaveButton.FontWeight = 'bold';

% Rois selectable info table
hRoisTable = uitable('Position',hRoisTAx.Position,'Units','normalized',...
    'CellSelectionCallback',@RoiSelect, 'CellEditCallback',@RoiEdit);
titles = cell(1,nfiles);
for i = 1:nfiles
    titles{i} = sprintf('rec %d', i);
end
nLinks = sum(linkMat~=0, 2);
% length(score)
% length(nLinks)
% [nLinks, score]
[~, sorted]=  sortrows([nLinks score], 'descend');
linkMat = [nLinks(sorted), linkMat(sorted,:), score(sorted,:)];
hRoisTable.ColumnWidth = [{35}, num2cell(repmat(35, [1, nfiles]))];
hRoisTable.ColumnName = [{'n links'}, titles, {'score'}];
hRoisTable.ColumnEditable = [false, true(1,nfiles), false];
hRoisTable.Data = num2cell(linkMat);

hClickedTable = uitable('Position',hClickTAx.Position,'Units','normalized');
hClickedTable.Data = num2cell(zeros(4,nfiles));
hClickedTable.ColumnName = [titles];
hClickedTable.ColumnWidth = [num2cell(repmat(35, [1, nfiles]))];
hClickedTable.RowName = [{'ROI'},{'row'},{'n-links'},{'score'}];

% Start the UI with nice background image
UpdateView

% Callback functions and such

function BackgroundView(source,~)
    % Callback function to change the background view
    if source.Value ~= viewToggle
        viewToggle = source.Value;
        UpdateView
    end
end


function SessionSelect(source, event)
    % Callback to select and deselect sessions
    if ~isempty(event.Indices)
        idx = sub2ind([nfiles 2], event.Indices(1), event.Indices(1,2));
        source.Data{idx} = ~source.Data{idx};
    end
    sShow = cell2mat(source.Data(:,1));
    sCont = find(cell2mat(source.Data(:,2)));
    UpdateView
    UpdateCont
end


function RoiSelect(~, event)
    % Callback to select and deselect ROI
    
    % Set the previous selection linewidth and alpha back to normal
    for i = 1:size(selCells,1)
        if length(selRoi{i})==1
            if selRoi{i}>0
                % Set properties of the contour of session, roi
                hCons{selFiles(i)}(selRoi{i}).Color(4) = contAlpha;
                hCons{selFiles(i)}(selRoi{i}).LineWidth = 1;
            end
        else % the nLinks or score column was clicked -> color all ROIs
            for j = 1:nfiles
                if selRoi{i}(j)>0
                    hCons{j}(selRoi{i}(j)).Color(4) = contAlpha;
                    hCons{j}(selRoi{i}(j)).LineWidth = 1;
                end
            end
        end
    end
    
    % Which ROIs were clicked? 
    selCells = event.Indices; % Selected cells
    selFiles = selCells(:,2)-1; % -1 because first column is nLinks
    selRoi = cell(1,length(selFiles)); % Which ROIs need to get highlighted  
    for i = 1:length(selFiles)
        if selFiles(i)>0 && selFiles(i)<=nfiles % A ROI in a session was clicked
             selRoi{i} = hRoisTable.Data{selCells(i,1),selCells(i,2)};
        else % multiple rois should be selected because cell nLinks or score was clicked
            % Get the ROI number out of the table
            for j = 1:nfiles
                selRoi{i}(j) = hRoisTable.Data{selCells(i,1),j+1};
            end
        end
    end
    
    % Set the clicked rois to fat contours
    for i = 1:size(selCells,1)
        if length(selRoi{i})==1 
            if selRoi{i}>0
                hCons{selFiles(i)}(selRoi{i}).Color(4) = 1;
                hCons{selFiles(i)}(selRoi{i}).LineWidth = 2;
            end
        else % the nLinks or score column was clicked -> color all ROIs
            for j = 1:nfiles
                if selRoi{i}(j)>0
                    hCons{j}(selRoi{i}(j)).Color(4) = 1;
                    hCons{j}(selRoi{i}(j)).LineWidth = 2;
                end
            end
        end
    end
    
    rShow = unique(event.Indices(:,1));
%     rois = linkMat(rShow, 2:end)
%     [~, idx] = find(rois)
%     idx = idx(1:size(rois,1))

    % With certain views the image changes as well
    if viewToggle == 2 || viewToggle == 4
        UpdateView
    end
    if autoZoom
        AutoZoomFunc
    end
end


function RoiEdit(source,event)
    % Edits the ROI after editing a cell and then also recalculates:
    % The new number of links & the new score for the edited row
    % A ROI was already present in the table, but no ROI is allowed to be
    % present twice in the table, so the previous spot/row should be deleted,
    % and also needs recalculated number of links and new score
    
    row = event.Indices(1);
    column = event.Indices(2);
    oldData = event.PreviousData;
    newData = event.NewData;
    
    if newData > 0 
        % Where this ROI was previously should be deleted
        oldRow = find(linkMat(:,column)==event.NewData);
        linkMat(oldRow,column) = 0;
        
        % Update the number of links for the row where the new ROI was
        % previously
        n = sum(linkMat(oldRow,2:end-1)>0);
        if n == 0
            linkMat(oldRow,:) = []; % delete old row if it is empty now
        else
            linkMat(oldRow,1) = n; % new number of links
            linkMat(oldRow,end) = OverlapScore(inRoi, linkMat(:,2:end-1), oldRow);
        end
        
    elseif isnan(newData)
        newData = 0;
    end
    
    if oldData > 0
        % Previous ROI gets removed, so it has be added as a single ROI
        linkMat(end+1,column) = oldData;
        linkMat(end,1) = 1;
    end
    
    % Update the number of links for the edited row
    linkMat(row,column) = newData;
    n = sum(linkMat(row,2:end-1)>0);
    if n == 0
        linkMat(row,:) = []; % delete old row if it is empty now
    else
        linkMat(row,1) = n; % new number of links
    end
    
    linkMat(row,end) = OverlapScore(inRoi, linkMat(:,2:end-1), row); % Recalculate score
    
    source.Data = num2cell(linkMat);
    jScrollPane = findjobj(source);
    jScrollPane.getVerticalScrollBar.setValue(row)
end



function UpdateView
    % Updates the image after a change is requested
    
    sShow2 = find(sShow);
    if ~isempty(sShow2) && ~isempty(rShow)
        colors2 = colors;
        nShow = length(sShow2);
        if nShow < 9
            colors2 = colors(sShow2,:);
        end
        
        % ViewToggle selection
        if viewToggle == 1 || viewToggle == 2 
            
            if viewToggle == 1 % Background images
                hImg.CData = CreateRGB2(BImg2, sShow2, colors2);
            else % Mask view!
                MasksShow = Masks2;
                for x = sShow2(:)' % ROI selection
                    idx = linkMat(rShow, x+1);
                    MasksShow{x}(~ismember(MasksShow{x},idx)) = 0;
                end
                hImg.CData = CreateRGB(MasksShow, sShow2, colors2, 'binary');
            end
            
            % legend Text
            for x = 1:nfiles % All strings need to check for update
                if ismember(x,sShow2)
                    hLegendText(x).String = names{x};
                else
                    hLegendText(x).String = {''};
                end
            end
            
        elseif viewToggle == 3 
            % Number of links image. color based on ROI number of links
            hImg.CData = CreateRGB(nLinksMask, sShow2, colors2);
            
            % legend Text
            for x = 1:nfiles
                if ismember(x,sShow2)
                    hLegendText(x).String = sprintf('ROI found across %d recordings', x);
                else
                    hLegendText(x).String = {''};
                end
            end
           
        elseif viewToggle == 4 
            % Mask pixel overlap determines intensity
            linkHeat = Masks2;
            for x = sShow2(:)' % ROI selection possible
                idx = linkMat(rShow, x+1);
                linkHeat{x}(~ismember(linkHeat{x},idx)) = 0;
            end
            linkHeat = cellfun(@logical, linkHeat(sShow2), 'UniformOutput', false);
            linkHeat = reshape(cell2mat(linkHeat), [dims, nShow]);
            hImg.CData = sum(linkHeat,3);
            colormap hot
            
            hotcolor = hot(nShow);
            for x = 1:nfiles
                if x <= nShow
                    hLegendText(x).Color = hotcolor(x,:);
                    hLegendText(x).String = sprintf('overlap of %d recordings',x);
                else
                    hLegendText(x).String = {''};
                end
            end
        end
        
        % Put the normal colors if the view is not linkview 2
        if viewToggle ~= 4
            for x = 1:nfiles
                hLegendText(x).Color = colors(x,:);
            end
        end
        
        if autoZoom
            AutoZoomFunc
        end
    else
        hImg.CData = repmat(zeros(dims),[1, 1, 3]);
    end
end


function hImgDown(~,event)
    % Plots text with fprintf haha
%     source
%     event
    point = round(event.IntersectionPoint);
    hClickedTable.Data = num2cell(zeros(4,nfiles));
    for i = 1:nfiles
        roi = Masks2{i}(point(2), point(1));
        if roi>0
            table = cell2mat(hRoisTable.Data(:,i+1));
            row = find(table==roi);
            hClickedTable.Data{1,i} = roi;
            if ~isempty(row)
    %             fprintf('recording %d, roi %d, row %d\n', i, roi, row)
                hClickedTable.Data{2,i} = row;
                hClickedTable.Data{3,i} = hRoisTable.Data{row,1}; % number of links
                hClickedTable.Data{4,i} = round(hRoisTable.Data{row,nfiles+2},1); % score
            end
    	end
    end
end


function AlphCont(source,~)
    % Adjust alpha value for contours
    contAlpha = round(source.Value,1);
    source.Value = contAlpha;
    UpdateCont
end


function UpdateCont
    % Update contours
    
    % Make the contours of some sessions visible
    for x = sCont(:)'
        for y = 1:contours(x).Cnt
            hCons{x}(y).Color(4) = contAlpha;
        end
    end
    
    % Make the other contours invisible
    others = 1:nfiles;
    others = others(~ismember(others, sCont));
    for x = others(:)'
        for y = 1:contours(x).Cnt
            hCons{x}(y).Color(4) = 0;
        end
    end
end


function AutoZoomFunc
    % zoom to the final selected contour
    [~, idx] = find(linkMat(rShow(end),2:end)>0);
    roi = linkMat(rShow(end),idx(1)+1);
    hImgAx.XLim = [min(contours(idx(1)).Con(roi).y)-30,...
                   max(contours(idx(1)).Con(roi).y)+30];
    hImgAx.YLim = [min(contours(idx(1)).Con(roi).x)-30,...
                   max(contours(idx(1)).Con(roi).x)+30];
end


function AutoZoomButton(source, ~)
    % Toggle autozooming
    autoZoom = source.Value;
    if autoZoom % if turned on, remember current axes zoom
    	oldZoom = [hImgAx.XLim, hImgAx.YLim];
        hZoomButton.FontWeight = 'bold';
        hZoomButton.ForegroundColor = [0 0 1];
    else
        hImgAx.XLim = oldZoom(1:2);
        hImgAx.YLim = oldZoom(3:4);
        hZoomButton.FontWeight = 'normal';
        hZoomButton.ForegroundColor = [0 0 0];
    end
end


function SaveRois(source, ~)
    % Saving the changes of the edited linkMat table
    
    hSaveButton.String = 'saving the table';
    assignin('base', 'linkMat2', linkMat(:,2:end-1));
    assignin('base', 'score', linkMat(:,1));
    
    try
        pause(1)
        hSaveButton.String = 'saved the edited table';
    catch
        fprintf('the savebutton disapeared before the text could be edited\n')
    end
end


end