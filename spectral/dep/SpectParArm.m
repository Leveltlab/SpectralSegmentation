function varargout = SpectParArm(varargin)
%
% SpectParArm(imgstack, sax) creates a popup where you can select the 
% "spectral parameters" (spar) that are used in the automatic ROI selection
% procedure.
% 
% 
% input:
%   imgstack: 3D double, with different images of the 2P data
%   sax: vector, spectral axis, should have same size as 3rd dimension of
%                imgstack
% 
% output: output is generated when clicking on the accept & exit button
%       spar: struct with the parameters to inform the auto ROI selection
% 
% Now the matlab standard info for guide GUIs.
% SPECTPARARM MATLAB code for SpectParArm.fig
%      SPECTPARARM, by itself, creates a new SPECTPARARM or raises the existing
%      singleton*.
% 
%      H = SPECTPARARM returns the handle to a new SPECTPARARM or the handle to
%      the existing singleton*.
%
%      SPECTPARARM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPECTPARARM.M with the given input arguments.
% 
%      SPECTPARARM('Property','Value',...) creates a new SPECTPARARM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SpectParArm_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SpectParArm_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
% 
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SpectParArm

% Last Modified by GUIDE v2.5 04-May-2020 22:05:33

% Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @SpectParArm_OpeningFcn, ...
                       'gui_OutputFcn',  @SpectParArm_OutputFcn, ...
                       'gui_LayoutFcn',  [] , ...
                       'gui_Callback',   []);
    if nargin && ischar(varargin{1})
        gui_State.gui_Callback = str2func(varargin{1});
    end

    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end
    % End initialization code - DO NOT EDIT
end


% --- Executes just before SpectParArm is made visible.
function SpectParArm_OpeningFcn(hObject, eventdata, h, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to SpectParArm (see VARARGIN)
    
    
    if nargin <= 3
        warning('seriously give more input plz')
        return
    elseif nargin > 3
        specImg = varargin{1};
        dims = size(specImg);
    end
    if nargin == 5
        sax = varargin{2};
    else
        sax = 1:dims(3);
        fprintf('Created placeholder spectral axis')
    end 
    
    global DISPLAY % to activate plotting in roisfromlocalmax.m
    DISPLAY = false; % global variable, which lowercase name is also a function..

    
    % Check if a spar is present in current folder.
    if exist('spar.mat', 'file')
        % Use the parameter values if spar is already present
        load('spar.mat')
        fprintf('using parameters from spar.m file in current folder\n')
        if isfield(spar, 'border')
            h.sparBorderBox.String = num2str(spar.border);
        end
        if isfield(spar, 'areasz')
            h.sparAreaMinBox.String = num2str(spar.areasz(1));
            h.sparAreaMaxBox.String = num2str(spar.areasz(2));
        end
        if isfield(spar, 'roundedness')
            h.sparRoundednessBox.String = num2str(spar.roundedness);
        end
        if isfield(spar, 'voxel')
            h.sparVoxelBox.String = num2str(spar.voxel);
        end
        if isfield(spar, 'cutOffHz')
            h.sparCutOffHzMaxBox.String = num2str(spar.cutOffHz);
        end
        if isfield(spar, 'cutoffcorr')
            h.sparCutOffCorrBox.String = num2str(spar.cutoffcorr);
        end
    end
    
%     % Allow the usual toolbars
    set(h.hGUI,'toolbar','figure');
    set(h.hGUI,'menubar','figure');
    
    % Change colorbar for 1st axes
    h.colorbarAxes1.YAxisLocation = 'right';
    h.colorbarAxes1.XTickLabel = '';
    h.colorbarAxes1.NextPlot = 'add';
    h.colorbarAxes1Img = imagesc(1, sax, flip(rot90(1:dims(3))), 'Parent', h.colorbarAxes1);
    h.colorbarAxes1.YLim = sax([1 end]);
    h.colorbarAxes1.YLabel.String = 'frequency (Hz)';
    
    % Colorbar axes for the 2nd axes
    h.colorbarAxes2.YAxisLocation = 'right';
    h.colorbarAxes2.XTickLabel = '';
    h.colorbarAxes2.NextPlot = 'add';
    h.colorbarAxes2Img = gobjects(1);
    h.colorbarAxes1.YLim = sax([1 end]);
    h.colorbarAxes2.YLabel.String = 'frequency (Hz)';
    
    
    colors = flipud(jet(length(sax)));
    
    % Plot the spectral in seperate color
    norma = true;
    specImCell = squeeze(num2cell(specImg, [1 2]));
    specImRGB = CreateRGB2(specImCell, colors, norma);
    h.axes1.NextPlot = 'add';
    h.axes1Img = imagesc(specImRGB, 'Parent', h.axes1);
    
    % Also add the imagesc to axes2 to add the img graphics object handle
    h.axes2.NextPlot = 'add';
    h.axes2Img = imagesc(specImRGB, 'Parent', h.axes2);
    
    colormap(colors);
    
    % Fill in the spar struct 
    % (might have different fields then loaded spar)
    spar = struct();
    spar.cutOffHzMin = 0;
    spar.cutOffHzMax = 0; % will be retrieved by the callback in a moment
    spar.border = str2double(h.sparBorderBox.String); 
    spar.areasz = [str2double(h.sparAreaMinBox.String), str2double(h.sparAreaMaxBox.String)];
    spar.roundedness = str2double(h.sparRoundednessBox.String);  
    spar.voxel = str2double(h.sparVoxelBox.String);
    spar.cutoffcorr = str2double(h.sparCutOffCorrBox.String);
    % Save the spar for future use in the GUI
    setappdata(h.hGUI, 'spar', spar)
    
    % Save the sax and specImg for future use in the GUI
    data.sax = sax;
    data.dims = dims;
    data.specImg = specImg;
    data.norma = norma; % normalization of colors?
    setappdata(h.hGUI, 'data', data)
    
    % Plot the spectral image with the cutoffHz selection
    sparCutOffHzMaxBox_Callback(h.sparCutOffHzMaxBox, 0, h)
    h = guidata(h.hGUI); % Update the handles, was edited by slider callback
    
    % Plot the border for the buffer area    
    h.borderLine = plot(nan, nan, 'Parent', h.axes1);
    guidata(h.hGUI, h)
    PlotBorderLine(spar.border, dims, h)
    h = guidata(h.hGUI); % Update the handles, was edited by PlotBorderLine
    
    % Plot the circle for the roudnedness example
    guidata(h.hGUI, h)
    h.roundednessText = text(dims(2)/20, (dims(1)/20)*3,...
                            'Minimum ROI roundedness:', 'Color', [1 1 1],...
                            'HorizontalAlignment', 'left');
    h.roundednessLine = plot(nan, nan, 'Parent', h.axes1);
    PlotRoundednessLine(spar, h)
    h = guidata(h.hGUI);
    
    % Plot the circles for the minimum and maximum size examples
    guidata(h.hGUI, h)
    h.areaSzPatch = plot(nan, nan, 'Parent', h.axes1);
%     h.areaSzMinLine = plot(nan, nan, 'Parent', h.axes1);
    h.areaSzText = text(dims(2)/20, dims(1)-dims(1)/20,...
                            'Allowed ROI size range:', 'Color', [1 1 1],...
                            'HorizontalAlignment', 'left');
    PlotAreaSzLine(spar, dims, h)
    h = guidata(h.hGUI);
    
    linkaxes([h.axes1 h.axes2],'xy')
    
    % Choose default command line output for SpectParArm
    h.output = hObject;
    
    % Update handles structure
    guidata(hObject, h);

    % UIWAIT makes SpectParArm wait for user response (see UIRESUME)
    % uiwait(handles.hGUI);
end


% --- Outputs from this function are returned to the command line.
function varargout = SpectParArm_OutputFcn(hObject, eventdata, h) 
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % h    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    varargout{1} = h.output;
end


%% graphical object Creation functions

% --- Executes during object creation, after setting all properties.
function box_CreateFcn(hObject, eventdata, h) %#ok<*DEFNU>
    % hObject    handle to sparBorderBox (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % h    empty - handles not created until after all CreateFcns called
    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end


%% Plotting functions
function PlotBorderLine(border, dims, h)
    % plot the border parameter. Outside the border no ROIs can be created
    % this border is a buffer area
    
    x = [border, dims(2)-border, dims(2)-border, border,         border];
    y = [border, border          dims(1)-border, dims(1)-border, border];
    delete(h.borderLine) % delete old line
    h.borderLine = plot(x, y, 'w', 'Parent', h.axes1);
    guidata(h.hGUI, h) % Update the handles
end


function PlotAreaSzLine(spar, dims, h)
    % Plots circles that have the minimum and maximum size for a ROI
    
    distx = dims(2)/10; % distance from border
    disty = dims(1) - dims(1)/10;
    
    % The circle with the minimum size
    r = sqrt(spar.areasz(1)/pi);
    t = 0:pi/25:2*pi;
    x = r .* cos(t) + distx;
    y = r .* sin(t) + disty;    
    % The circle with the maximum size
    r = sqrt(spar.areasz(2)/pi); % radius of circle
    t = 0:pi/25:2*pi;
    x2 = r .* cos(t) + distx;
    y2 = r .* sin(t) + disty;
    delete(h.areaSzPatch)
	h.areaSzPatch = patch([x x2],[y y2], 'w', 'faceAlpha', 0.6, 'Parent', h.axes1);
%     fprintf('area of circle = %.3f\n', polyarea(x2,y2))


%     delete(h.areaSzMinLine)
%     h.areaSzMinLine = plot(x,y, 'w', 'linewidth', 2, 'Parent', h.axes1);
    guidata(h.hGUI, h)
end


function PlotRoundednessLine(spar, h)
    % plot a circle that is roughly as round as the minimum roundedness
    % permits
    r = sqrt(spar.areasz(2)/pi);
    t = 0:pi/25:2*pi;
    rounded = spar.roundedness;
    % morph = approximate circle squeezing to get requested roundedness
    morph = 1./(rounded+0.05)+0.05 - sin(2*pi * rounded)* 1/pi;
    dims = [h.axes1.YLim(2), h.axes1.XLim(2)];
    x = r .* cos(t) .* morph + dims(2)/10;
    y = r .* sin(t) ./ morph + dims(1)/10;

%     fprintf('requested roundndess %.3f, got %.4f, morph=%.3f\n',...
%                                   rounded, perimarea(x,y),morph)
    delete(h.roundednessLine)
    h.roundednessLine = plot(x, y, 'w', 'Linewidth',2, 'Parent', h.axes1);
    guidata(h.hGUI, h) % Update the handles
end


function PlotCutOffHz(data, selectedHz, selectedSax, h)
    % Plots the color image of the spectral components and adjust the
    % colorbar
    % plot the updated spectral frequency
    nSelectedHz = sum(selectedHz);
    colors = flipud(jet(nSelectedHz));
    
    specImgCell = squeeze(num2cell(data.specImg(:,:,selectedHz), [1 2]));
    specImgRGB = CreateRGB2(specImgCell, colors, data.norma);
    delete(h.axes2Img) % Delete previous image
    h.axes2Img = imagesc(specImgRGB, 'Parent', h.axes2);

    % Adjust colorbar
    delete(h.colorbarAxes2Img)
    h.colorbarAxes2Img = imagesc(1, selectedSax, flip(rot90(selectedSax)), 'Parent', h.colorbarAxes2);
    h.colorbarAxes2.YDir = 'normal';
    h.colorbarAxes2.YLim = selectedSax([1 end]);
    h.colorbarAxes2.Colormap = colors;
    
    guidata(h.hGUI, h)
end

%% Callback functions

function sparBorderBox_Callback(hObject, eventdata, h)
    % hObject    handle to sparBorderBox (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % h    structure with handles and user data (see GUIDATA)
    spar = getappdata(h.hGUI, 'spar');
    data = getappdata(h.hGUI, 'data');
    
    border = str2double(hObject.String);
    if (spar.border ~= border) % border should have changed
        % border shouldn't be too big, or negative
        if border < min(data.dims([1 2]))/2 && border >= 0
            % Plot the new requested border
            PlotBorderLine(border, data.dims, h)
            h = guidata(h.hGUI); % Update the handles
            
            % Save the border
            spar.border = border;
            setappdata(h.hGUI, 'spar', spar)
        else
            % Reset the users input
            hObject.String = spar.border;
            fprintf('border input was not good (maybe too large or nan)\n')
        end
    end
end


function sparAreaMinBox_Callback(hObject, ~, h)
    % Set the minimum ROI area value
    spar = getappdata(h.hGUI, 'spar');
    data = getappdata(h.hGUI, 'data');
    spar.areasz(1) = str2double(hObject.String);
    setappdata(h.hGUI, 'spar', spar)
    
    PlotAreaSzLine(spar, data.dims, h)
	h = guidata(h.hGUI); % Update the handles
end


function sparAreaMaxBox_Callback(hObject, ~, h)
    % Set the maximum ROI area value
    spar = getappdata(h.hGUI, 'spar');
    data = getappdata(h.hGUI, 'data');
    spar.areasz(2) = str2double(hObject.String);
    setappdata(h.hGUI, 'spar', spar)
    
    PlotAreaSzLine(spar, data.dims, h)
	h = guidata(h.hGUI); % Update the handles
end


function sparRoundednessBox_Callback(hObject, ~, h)
    % hObject    handle to sparRoundednessBox (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % h    structure with handles and user data (see GUIDATA)
    spar = getappdata(h.hGUI, 'spar');
    
    roundedness = str2double(hObject.String);
    
    % If roundedness value is nonsense ignore input and reset to spar value
    if roundedness < 0 || roundedness > 1 || isnan(roundedness)
        hObject.String = num2str(spar.roundedness);
        return
    end
    
    % If roudnedness changed update the plot
    if spar.roundedness ~= roundedness 
        spar.roundedness = roundedness;
        PlotRoundednessLine(spar, h)
        h = guidata(h.hGUI); % Update the handles
        setappdata(h.hGUI, 'spar', spar)
    end
end


function sparVoxelBox_Callback(hObject, ~, h)
    % Update the spar when the voxel value gets edited
    spar = getappdata(h.hGUI, 'spar');
    spar.voxel = str2double(hObject.Value);
    setappdata(h.hGUI, 'spar', spar)
end



function sparCutOffHzMinBox_Callback(hObject, ~, h)    
    % Sets the minimum spectral frequency to use in ROI search
    spar = getappdata(h.hGUI, 'spar');
    data = getappdata(h.hGUI, 'data');
    sax = data.sax;
    cutOffHzMin = str2double(hObject.String);
    
    cutOffHzMax = spar.cutOffHzMax;
    
    % If the cutoffHz is changed calculate and save the update
    if spar.cutOffHzMin ~= cutOffHzMin
        % Find which frequencies are now selected
        selectedHz = (sax < cutOffHzMax) & (sax > cutOffHzMin);
        selectedSax = sax(selectedHz);
        nSelectedHz = sum(selectedHz);
        
        if nSelectedHz > 0
            PlotCutOffHz(data, selectedHz, selectedSax, h)
            h = guidata(h.hGUI); % Update the handles
            
            % Save the updated spar and handle
            spar.cutOffHzMin = cutOffHzMin;
            setappdata(h.hGUI, 'spar', spar)
        else
            fprintf('no spectral components selected with this input, aborting change\n')
            hObject.String = num2str(spar.cutOffHzMin);
        end
    end
end


function sparCutOffHzMaxBox_Callback(hObject, ~, h)
    % Sets the maximum spectral frequency to use in ROI search
    spar = getappdata(h.hGUI, 'spar');
    data = getappdata(h.hGUI, 'data');
    sax = data.sax;
    cutOffHzMin = spar.cutOffHzMin;
    
    cutOffHzMax = str2double(hObject.String);
    
    % If the cutoffHz is changed calculate and save the update
    if spar.cutOffHzMax ~= cutOffHzMax
        % Find which frequencies are now selected
        selectedHz = (sax < cutOffHzMax) & (sax > cutOffHzMin);
        selectedSax = sax(selectedHz);
        nSelectedHz = sum(selectedHz);
        
        if nSelectedHz > 0
            PlotCutOffHz(data, selectedHz, selectedSax, h)
            h = guidata(h.hGUI); % Update the handles
            
            % Save the updated spar and handle
            spar.cutOffHzMax = cutOffHzMax;
            setappdata(h.hGUI, 'spar', spar)
        else
            fprintf('no spectral components selected with this input, aborting change\n')
            hObject.String = num2str(spar.cutOffHzMax);
        end
    end
end


% --- Executes on button press in DisplayCheckBox.
function DisplayCheckBox_Callback(hObject, ~, h)
    % Sets the global DISPLAY boolean
    global DISPLAY
    DISPLAY = hObject.Value;
end


function sparCutOffCorrBox_Callback(hObject, ~, h)
    % edits spar when cutoffcorr gets edited
    spar = getappdata(h.hGUI, 'spar');
    spar.cutOffCorr = str2double(hObject.Value);
    setappdata(h.hGUI, 'spar', spar)
end


% --- Executes on button press in acceptButton.
function acceptButton_Callback(hObject, ~, h)
    % Saves the spar and exits the UI
    
    spar = getappdata(h.hGUI, 'spar');
    assignin('base', 'spar', spar);
    
    % Saving the spar into the current folder
    save('spar.mat', 'spar')
    fprintf('saved the spar into current folder and workspace\n')
    
    % exit
    closereq();
end
