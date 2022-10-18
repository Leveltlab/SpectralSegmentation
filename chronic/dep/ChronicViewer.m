function ChronicViewer(BImgs, Masks, names, nLinksMask, linkMat, contours, score, inRoi)
% UI for going through chronically linked recordings
% 
% It is possible to edit the match matrix. Editing the match matrix will
% recalculate relevant variables: number of links in the match & overlap
% score of the match.
% When editing the match matrix, duplicate ROI numbers will be removed 
% automatically.
% 
% ChronicViewer2(BImgs, Maskss, names, nLinksMask, linkMat, contours, score, inRoi)
% BImgs: [1 x n] cell array with 2D double. Registered spectral images.
% Masks: [1 x n] cell array with 2D double. Registered ROI masks.
% names: [1 x n] string names of the different recordings.
% nLinksMask: [1 x n] cell array with 2D doubles: ROI masks, with nr of
%                      links instead of ROI number per ROI
% linkMat: [nmatches x n] double. Which ROI in each match
% contours: [1 x n] struct with registered contour information (PP)
% score: [nmatches x 1] double. The overlap score of each match
% inRoi: [1 x n] cell array. all overlap info for all ROIs in all recordings
% 
% n = number of linked recordings
% 
% 
% The masks view only shows masks that are chronically linked to another
% ROI in at least 1 other recording.
% 
% 
% Leander de Kraker
% 2021-5-12
% 2022-9-24: Added ability to deselect ROIs by pressing the empty 1st row
% 2022-9-28: to-edit-ROI-row takes into account changing row from editing
% 
%

% The height of the table individual ROIs according to the java properties
tableRowHeight =  17.8349; 

nfiles = length(Masks);

dims = size(Masks{1});


sShow = true(1, nfiles);
sCont = 1:nfiles;
rShow = 1:length(linkMat);
selCells = [2,2; 3,2]; % which cells of the table are selected
selFiles = [1, 1]; % Which files does that selection correspond to
selRoi = [{1} {2}]; % Which ROIs does that selection correspond to
contAlpha = 0.4;
autoZoom = false;
oldZoom = [1, dims(1), 1, dims(2)];
viewToggle = 1;
colors = flipud(cmapL([0 0 1; 0 1 1; 0 1 0; 1 1 0; 1 0 0; 1 0 1], nfiles));

% Activate tight subplotting
% [subplot margin side&top],[figure bottomspace,topspace],[leftspace,rightspace]
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.01], [0.025 0], [0.025 0]);

hFig = figure('name','Chronic Viewer', 'units', 'normalized', 'position', [0.1 0.05 0.825 0.85]);

hImgAx = subplot(1, 3, [1 2]); hold on
hImgAx.ButtonDownFcn = @hImgDown;
hImgAx.YDir = 'reverse';
hImg = imagesc(zeros(dims([1 2])),'hittest','off');
axis([1, dims(2), 1, dims(1)])

% Draw contours and legend
hCons = {};
hLegendText = gobjects(nfiles,1);
for i = 1:nfiles
    % Draw contours and save handles
    hCons{i} = gobjects(1,contours(i).Cnt);
    for j = 1:contours(i).Cnt
        hCons{i}(j) = plot(contours(i).Con(j).x, contours(i).Con(j).y,...
            'color',[colors(i,:), contAlpha],'hittest','off');
    end
    
    % Write colored legend text
    names{i} = strrep(names{i},'_',' '); % replace underscores with spaces
    hLegendText(i) = text(3,12+i*12,names{i},'color',colors(i,:));
end

% Positions of the other UI elements
hRoisTAx= subplot(16, 3, [3 24]); % chronic table table axis
hInfoAx = subplot(16, 3, [27 30]);  % Info axis
hViewAx = subplot(16, 6, 89); % select background image list axis
hAlphAx = subplot(16, 6, 90); % Alpha slider axis
hZoomAx = subplot(16, 6, 95); % Auto zoom button axis
hSaveAx = subplot(16, 6, 96); % Save button axis
hClickTAx= subplot(16, 3, [33 36]); % clicked on image roi information table axis
hSessTAx = subplot(16, 3, [39 42]); % which sessions are selected table axis
hRoisTAx.Units = 'pixels';
hInfoAx.Units  = 'pixels';
hViewAx.Units  = 'pixels';
hAlphAx.Units  = 'pixels';
hZoomAx.Units  = 'pixels';
hSaveAx.Units  = 'pixels';
hClickTAx.Units= 'pixels';
hSessTAx.Units = 'pixels';
hRoisTAx.Visible = 'off';
hInfoAx.Visible  = 'off';
hViewAx.Visible  = 'off';
hAlphAx.Visible  = 'off';
hZoomAx.Visible  = 'off';
hSaveAx.Visible  = 'off';
hClickTAx.Visible= 'off';
hSessTAx.Visible = 'off';

% Select sessions table
hSessionTable = uitable('Position',hSessTAx.Position,'Units','normalized',...
    'CellSelectionCallback',@SessionSelect);
hSessionTable.RowName = [{'img'},{'cont'}];
hSessionTable.ColumnWidth = [num2cell(repmat(35, [1, nfiles]))];
hSessionTable.Data = repmat({true},[2,nfiles]);

% View dropdown menu
hViewPopup = uicontrol('Style', 'popup',...
    'String', {'Backgrd','Masks selected','Masks all','nLinks','mask overlap'},...
    'Position', hViewAx.Position,'units','normalized','Callback',@BackgroundView);  

% Alpha value slider for contours
hAlph = uicontrol('Style','slider', 'Min',0, 'Max',1, 'Value',contAlpha,...
    'sliderstep',[0.1 0.1], 'Position',hAlphAx.Position,'Callback', @AlphCont);
hAlph.Position(4) = 15;
hAlphTitle = uicontrol('Style','text',...
    'position',hAlphAx.Position,'horizontalAlignment','left',...
    'string','Contours alpha (opacity)');
hAlphTitle.Position([2, 4]) = [hAlphTitle.Position(2)+tableRowHeight, 15];
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

% TABLE Chronic matching ROI info table
hRoisTable = uitable('Position',hRoisTAx.Position,'Units','normalized',...
            'CellSelectionCallback',@RoiSelect, 'CellEditCallback',@RoiEdit);
titles = cell(1,nfiles);
for i = 1:nfiles
    titles{i} = sprintf('rec %d', i);
end
nLinks = sum(linkMat~=0, 2);
[~, sorted]=  sortrows([nLinks score], 'descend');
linkMat = [nLinks(sorted), linkMat(sorted,:), score(sorted,:)];
linkMat = [zeros(1,size(linkMat,2)); linkMat];
hRoisTable.ColumnWidth = [{35}, num2cell(repmat(35, [1, nfiles]))];
hRoisTable.ColumnName = [{'n links'}, titles, {'score'}];
hRoisTable.ColumnEditable = [false, true(1,nfiles), false];
hRoisTable.Data = num2cell(linkMat);


% TABLE that shows info on ROIs which you clicked on in main image
hClickedTable = uitable('Position',hClickTAx.Position,'Units','normalized',...
                'CellSelectionCallBack', @clickTableCallback);
hClickedTable.Data = num2cell(zeros(4,nfiles));
hClickedTable.ColumnName = [titles];
hClickedTable.ColumnWidth = [num2cell(repmat(35, [1, nfiles]))];
hClickedTable.RowName = [{'ROI'},{'row'},{'n-links'},{'score'}];


% Text that says something about the UI status
hInfoText = uicontrol('Style', 'text', 'String', 'hello',...
            'Position',hInfoAx.Position, 'Units','normalized');
hInfoText.HorizontalAlignment = 'Left';

% Start the UI with nice background image
UpdateView
RefreshInfoText


% Callback functions and such % % % % % % % % % % % % % % % % % % % % % %

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
        idx = sub2ind([2, nfiles], event.Indices(1), event.Indices(1,2));
        source.Data{idx} = ~source.Data{idx};
    end
    sShow = cell2mat(source.Data(1,:));
    sCont = find(cell2mat(source.Data(2,:)));
    UpdateView
    UpdateCont
end


function RoiSelect(~, event)
    % Callback to select and deselect ROI
    if ~isempty(event.Indices)
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
        
        % Which ROIs were clicked? Update selected ROIs
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
        if length(rShow)==1 && rShow==0
            rShow = 1:length(linkMat);
        end
        
        RefreshInfoText
        
        % With certain views the image changes as well
        if viewToggle == 2 || viewToggle == 5
            UpdateView
        end
        if autoZoom
            AutoZoomFunc
        end
    end
end


function RoiEdit(source,event)
    % Edits the ROI after editing a cell and then also recalculates:
    % The new number of links & the new score for the edited row
    % A ROI was already present in the table, but no ROI is allowed to be
    % present twice in the table, so the previous spot/row should be deleted,
    % and also needs recalculated number of links and new score
    
    % Remember scrollposition
    jscrollpane = javaObjectEDT(findjobj(source));
    viewport    = javaObjectEDT(jscrollpane.getViewport);
    P = viewport.getViewPosition();
    
    
    row = event.Indices(1);
    column = event.Indices(2);
    oldData = event.PreviousData;
    newData = event.NewData;
    
    
    if ((P.y-1) /tableRowHeight) > row 
        % scroll position is bigger than row to go to, which has to mean it
        % is not in view: change scroll position
        P.y = (rowToGoTo-1) *tableRowHeight;
    end
    
    if newData > 0 
        % Where this ROI was previously should be deleted
        oldRow = find(linkMat(:,column)==event.NewData);
        linkMat(oldRow,column) = 0;
        
        % Update the number of links for the row where the new ROI was
        % previously
        n = sum(linkMat(oldRow,2:end-1)>0);
        if n == 0
            linkMat(oldRow,:) = []; % delete old row if it is empty now
            if row > oldRow % If that deleted old row is deleted, shift the selected row up.
                row = row - 1;
                selCells(1) = selCells(1) - 1;
                RefreshInfoText
            end
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
    
    % If the added ROI was a new match enable new match by adding empty entry
    if row == 1
        linkMat = [zeros(1, size(linkMat,2)); linkMat];
        selCells(1) = selCells(1)+1; % the selected cell is thus shifted
        RefreshInfoText
    end
    
    source.Data = num2cell(linkMat);    
    
    % Reset Set scroll position
    drawnow() % This is necessary to ensure the view position is set after matlab hijacks it
    viewport.setViewPosition(P);
    
    eventToGive.Indices = selCells;
    RoiSelect(nan, eventToGive)
end


function clickTableCallback(source, event)
    % Clicked in the clickedTable, 
    % Option 1: If clicked on ROI number add ROI in RoiTable
    %   but only if clicked on ROI number in clickedTable and if only
    %   one match is selected in the RoiTable.
    % Option 2: If clicked on row scroll to match in RoiTable
    
    clickedCells = event.Indices; % Selected cells
    
    if size(clickedCells,1) == 1 % Only do things if clicked on one cell

        if length(unique(selCells(:,1)))==1 && clickedCells(1)==1
            % Add ROI into roiTable
            RoiToAdd = source.Data{clickedCells(1), clickedCells(2)};
            if RoiToAdd ~= 0
                row = selCells(1);
                rec = clickedCells(2)+1; % +1 because first row chronic table is nlinks.
                eventToGive.Indices = row;
                eventToGive.Indices(2) = rec;
                eventToGive.PreviousData = hRoisTable.Data{row,rec};
                eventToGive.NewData = source.Data{clickedCells(1), clickedCells(2)};
                
                if eventToGive.PreviousData ~= eventToGive.NewData
                    RoiEdit(hRoisTable, eventToGive)
                end
            end
        elseif clickedCells(1)==2
            % Scroll to clicked row
            rowToGoTo = source.Data{clickedCells(1), clickedCells(2)};
            if rowToGoTo ~= 0
                jscrollpane = javaObjectEDT(findjobj(hRoisTable));
                viewport    = javaObjectEDT(jscrollpane.getViewport);
                P = viewport.getViewPosition();
                P.y = (rowToGoTo-1) *tableRowHeight;
                viewport.setViewPosition(P);
            end
        end
    end
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
        if ismember(viewToggle, [1 2 3])
            
            if viewToggle == 1 % Background images
                hImg.CData = CreateRGB2(BImgs(sShow2), colors2);
            else % Mask view
                MasksShow = Masks;
                if rShow>1 & viewToggle == 2
                    for x = sShow2(:)' % ROI selection
                        idx = linkMat(rShow, x+1);
                        MasksShow{x}(~ismember(MasksShow{x},idx)) = 0;
                        MasksShow{x}(MasksShow{x}>0) = 1;
                    end
                end
                hImg.CData = CreateRGB2(MasksShow(sShow2), colors2);
                
            end
            
            % legend Text
            for x = 1:nfiles % All strings need to check for update
                if ismember(x,sShow2)
                    hLegendText(x).String = names{x};
                else
                    hLegendText(x).String = {''};
                end
            end
            
        elseif viewToggle == 4
            % Number of links image. color based on ROI number of links
            hImg.CData = CreateRGB2(nLinksMask(sShow2), colors2);
            
            % legend Text
            for x = 1:nfiles
                if ismember(x,sShow2)
                    hLegendText(x).String = sprintf('ROI found across %d recordings', x);
                else
                    hLegendText(x).String = {''};
                end
            end
           
        elseif viewToggle == 5
            % Mask pixel overlap determines intensity
            linkHeat = Masks;
            for x = sShow2(:)' % ROI selection possible
                idx = linkMat(rShow, x+1);
                linkHeat{x}(~ismember(linkHeat{x},idx)) = 0;
            end
            linkHeat = cellfun(@logical, linkHeat(sShow2), 'UniformOutput', false);
            linkHeat = cat(3, linkHeat{:});
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
        if viewToggle ~= 5
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
    % Highlights ROI in the table and main image after clicking on main image
    point = round(event.IntersectionPoint);
    hClickedTable.Data = num2cell(zeros(4,nfiles));
    for i = 1:nfiles
        roi = Masks{i}(point(2), point(1));
        if roi>0
            table = cell2mat(hRoisTable.Data(:,i+1));
            row = find(table==roi);
            hClickedTable.Data{1,i} = roi;
            if ~isempty(row)
    %             fprintf('recording %d, roi %d, row %d\n', i, roi, row)
                hClickedTable.Data{2,i} = row;
                hClickedTable.Data{3,i} = hRoisTable.Data{row,1}; % number of links
                hClickedTable.Data{4,i} = round(hRoisTable.Data{row,nfiles+2},3); % score
            end
    	end
    end
end


function RefreshInfoText
    % Refresh the info text
    if length(unique(selCells(:,1)))==1 
        str = sprintf('match   %4d   selected    ', selCells(1));
    else
        str = 'If you select one match ';
    end
    str = {[str...
        'in the "match matrix ↑" you can add a ROI into that match '...
        'by clicking on the ROI number in the "clicked table ↓"'...
        'Get ROIs into the "clicked table ↓" by clicking on ROIs in the "main image".']};
    str{2} = ['You can move to the match matrix row by clicking on the row nr in the "clicked table ↓".'];
    str{3} = ['You can remove ROIs in the match matrix by deleting or typing 0 '...
              'in the desired ROI entry in the "match matrix ↑".'];
    str{4} = ['You can also edit the match matrix by typing a ROI nr in the desired ROI entry'...
              'in the "match matrix" itself.'];
    hInfoText.String = str;
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
    [~, idx] = find(linkMat(rShow(end),2:end-1)>0);
    if ~isempty(idx)
        roi = linkMat(rShow(end),idx+1);
        x = [];
        y = [];
        for i = 1:length(idx)
            x = cat(2, x, contours(idx(i)).Con(roi(i)).x);
            y = cat(2, y, contours(idx(i)).Con(roi(i)).y);
        end
        hImgAx.XLim = [min(x)-35, max(x)+35];
        hImgAx.YLim = [min(y)-35, max(y)+35];
    end
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
    assignin('base', 'linkMat', linkMat(2:end, 2:end-1));
    assignin('base', 'score', linkMat(:,1));
    
    try
        pause(1)
        hSaveButton.String = 'saved the edited table';
    catch
        fprintf('the savebutton disapeared before the text could be edited\n')
    end
end


end