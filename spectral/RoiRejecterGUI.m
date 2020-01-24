function varargout = RoiRejecterGUI(varargin)
%
% RoiRejecterGUI enables the visualisation of two-photon data, and it
% can show the signal of any place in the two-photon microscopy image over
% time.
%
% when no input is given the RoiRejecter will ask for the spectral file,
%   and the transposed data file
% In case that input is given, the input should be:
% RoiRejecterGUI(Mask, PP, SPic, Transfile, Sax)
%   input 1: Mask, the 2D matrix which says which pixels are occupied by
%               which ROI
%   input 2: PP, struct with contour information, created by getspectrois.m
%   input 3: SPic, SPectral images components, 3D matrix
%   input 4: The filename of the transposed datafile
%   input 5: Sax, the spectral axis, as created by spectral.m, that should
%            correspond to the SPic variable
%
% 
% The functionality of the RoiRejecter can be accessed with the 6 different
% tabs on the right.
% 
% ROIs can be refined by:
%	deletion: either by sliders that select based on ROI properties
%             or by manual clicking.       TAB 'Reject ROIs'
%	splitting (with pixel k-means clustering). Click on a ROI that
%             needs to be split.           TAB 'Split ROIs'
%	creation: click on a place in the image and the computer will try
%             to find a contour, the contour will be based on the
%             current background image in the RoiRejecter, in case that
%             that is a color image it will average the three colorchannels
%             to create a grayscale image. TAB 'CreateROI'
%	manual creation: manually draw a ROI by clicking on the main
%                    background image.     TAB 'Manual ROI'
%
% Data can be visualised in the following ways:
% FLUORESCENCE TIMETRACE (signal in bottom main axis)
%   Show the timetrace fluorescence data of an ROI: 
%   	in TAB 'Reject ROIs' toggle 'Plot ROI signal on ROI toggling', 
%       and then white or blacklist a ROI.
%   Show the data from a square selection of the data:
%       in TAB 'Show data' toggle the 'select signal' button
%   
% MAIN IMAGE OPTIONS
%   Load in any variable (that has the same size as the main background 
%       image), by selecting the 'variable from workspace' option in the
%       background view option in TAB 'Show data'
%       In case you load in a 3D matrix a colormap will be created, that is
%       defined in the function BackGrdView
%   Comparing the current recording with other recordings is possible
%       by selection the 'chronic' option in the background view option 
%       in the TAB 'show data'. The background images in the chronic
%       file will be averaged together, and then registered to the current
%       background image in the RoiRejecter.
%   the ROIs can be made colorful by selecting 'inner ROI corr' in the
%       background view option in the TAB 'show data'. This will do the
%       same analysis that 'Split ROIs' will do, for every ROI. Calculation
%       can take a while.
%
%
%
%      H = ROIREJECTERGUI returns the handle to a new ROIREJECTERGUI or the handle to
%      the existing singleton*.
%
%      ROIREJECTERGUI('Property','Value',...) creates a new ROIREJECTERGUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to RoiRejecterGUI_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      ROIREJECTERGUI('CALLBACK') and ROIREJECTERGUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in ROIREJECTERGUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
% Made by: Leander de Kraker
% 2018

% Last Modified by GUIDE v2.5 04-Oct-2019 19:28:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RoiRejecterGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @RoiRejecterGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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

%% Opening and creation functions %%
% --- Executes just before RoiRejecterGUI is made visible.
function RoiRejecterGUI_OpeningFcn(hObject, ~, h, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    
    % nargin == varargin + hObject, evendata and handles.
    if nargin == 6 || nargin == 7 || nargin == 8 || nargin == 9
        % Use variables if given 
        Mask = varargin{1};
        PP = varargin{2};
        SPic = varargin{3};
        if nargin > 6
            datafile = varargin{4};
        else
            [fn, pn] = uigetfile('*Trans*.dat', 'Select the file with 2P data');
            datafile = [pn fn];
        end
        if nargin == 8
            Sax = varargin{5};
        else
            Sax = 1:size(SPic,3);
        end
        if nargin == 9
            SpatialCorr = varargin{6};
        end
        
        if ~isstruct(PP) % Little checks for variable correctness
            warning('\n\nWrong variables are given probably! PP should be struct\n\n')
            return
        end
        if ~isnumeric(Mask)
            warning('\n\nWrong variables are given probably! Mask should be numeric\n\n')
            return
        end
        if ~isnumeric(SPic)
            warning('\n\nWrong variables are given probably! SPic should be numeric\n\n')
            return
        end
    else
        % load data otherwise
        [fn, pn] = uigetfile('*SPSIG.mat');
        fprintf('Loading...\n')
        SPSIGfile = [pn fn];
        load(SPSIGfile,'Mask','PP','Sax','SPic', 'SpatialCorr');
        fprintf('Done loading\n')
        
        % Try to find the 2P transposed datafile based on given SPSIG file
        datafile = [pn fn(1:end-9) 'Trans.dat'];
        if exist(datafile,'file') ~= 2 % if not found, request it
            [fn, pn] = uigetfile('*Trans*.dat', 'Select the file with 2P data');
            if any(fn ~= 0)
                datafile = [pn fn];
            else
                fprintf('cannot run RoiRejecterGUI without transposed data file\n')
                return
            end
        end
    end
    
    % Memorymap the 2P signal file
    [sbxt, dim, freq] = transmemap(datafile);
    
    % Some colors we will use
    colors = struct();
    colors.parula = parula(256);
    colors.parula(:,2) = (colors.parula(:,2)-0.15) * 1.15;
    
    % switches are things that determine what to do
    switches = struct();
    switches.currentMajor = 1; % start with ROI rejector panel (panel 1).
    switches.currentMinor = 'none'; % which minor function button is on
    switches.viewToggle = 0; % which background view to show on mainAxes
    switches.ThresToggle = h.thresToggle.Value; % 1 Which of the two threshold sliders to use
    switches.ThresCorToggle = h.thresCorToggle; % 2 threshold via mean correlation value
    switches.alph    = h.AlphaSlider.Value; % alpha value for ROI contours
    switches.plotListing = h.plotListing.Value; % Plot signal of clicked ROI? false default
    switches.roiCorrIm = false; % roi correlation image is not made yet.
    switches.minSize = round(min(PP.A)-1); % The current minimum ROI size for rejection
    switches.maxSize = round(max(PP.A),-1)+10; % set maximum ROI size value above biggest ROI
    switches.originalROIs = PP.A; % The current list of ROI sizes (to check changes before save)
    switches.selSize = h.selSize.Value; % size of data selection
    switches.selNumMax = h.numSelectSlider.Value; % Number of max data tracers
    switches.selNumCur = 1; % Current datatrace to draw
    switches.ylims = zeros(2); % the extreme values of selected signal
    switches.splitStarted = false; % Is an ROI being split at the moment?
    switches.splitRoi = 0; % Which ROI is being split at the moment
    switches.nCluster = h.nClusterSlider.Value; % in how many clusters to split the ROI
    switches.creationAllow = false; % allow creation of new ROI with current data
    switches.creationPos = [100, 100]; % position where the user wants a new ROI
    switches.creationThres = 80; % threshold for percentile pixelvalue new ROI inclusion
    switches.chronic = false;
    switches.backgrdCLim = [0 1]; % color axis limits of the background image in mainAx
    switches.activeSignalAx = 2;
    
    
    % Setting more useful limits for the  'minimum size rejection' slider
    h.minSizeSlider.Min = switches.minSize; % minimum value
    h.minSizeSlider.Max = switches.maxSize; % maximum value
    h.minSizeSlider.Value = switches.minSize; % current value
    h.minSizeTitle.String = {sprintf('Minimum size: %4.0f px', switches.minSize)};

    % The maximum ROI sizes can be very different per dataset, so the
    % slider value gets edited here
    h.maxSizeSlider.Max = switches.maxSize;
    h.maxSizeSlider.Value = switches.maxSize;
    h.maxSizeTitle.String = {sprintf('Maximum size: %4.0f px', switches.maxSize)};
    
    sliderNames = {'minSize', 'maxSize', 'threshold' 'innerCorr', 'roundedness'};
    sliderDefaults = nan([1, length(sliderNames)]);
    switches.sliderSettings = array2table(sliderDefaults, 'variableNames',sliderNames);
    
    % idx are ROI indexes that tells which ROIs to reject
    idx = struct();
    idx.White = false(1,PP.Cnt);
    idx.Black = false(1,PP.Cnt); 
    idx.Size  = false(1,PP.Cnt);
    idx.Thres = false(1,PP.Cnt);
    idx.ThresCor = false(1,PP.Cnt);
    idx.Round = false(1,PP.Cnt);
    
    Sax = Sax(2:end);
    imgStackT = permute(SPic(:,:,2:end),[2 1 3]);
    nSpect = size(imgStackT, 3);
    % If the spectral profiles are not yet there, calculate them
    if ~isfield(PP,'SpectralProfile')
        PP.SpecProfile = zeros(nSpect, PP.Cnt); % For each ROI save the spectral density values
        for i = 1:PP.Cnt
            % Retrieve the average spectral profile of the ROI
            ROIi = Mask == i; 
            npixels = sum(ROIi(:)); % number of pixels for this ROI
            ROIi3D = repmat(ROIi, [1,1,nSpect]);
            specProfile = reshape(imgStackT(ROIi3D), [npixels, nSpect]);
            PP.SpecProfile(:,i) = mean(log(specProfile));
        end
        % Calculate which frequency had the highest spectral density values for
        [~, peakFreq] = max(PP.SpecProfile);
        PP.peakFreq = Sax(peakFreq);
    end
    % If the creation method for the ROIs doesn't exist yet, initialize
    if ~isfield(PP,'creationMethod')
        % assume the ROIs that are present were created automatically
        PP.creationMethod = repmat({'auto'}, [PP.Cnt, 1]);
        warning("PP.creationMethod didn't exist, assuming existing ROIs were automatically created")
    end
    % If the SpatialCorr is not there yet, initialze empty placeholder and warn user
    if ~exist('SpatialCorr', 'var')
        SpatialCorr = zeros(size(Mask));
        warning('SpatialCorr not present, calculate it later with SpatialCorrCalcRun.m!')
    end
    if ~isfield(PP, 'Rvar') % if Rvar doesn't exist, warn user and initialze zeros
        PP.Rvar = zeros(PP.Cnt,1);
        warning('PP.Rvar not present, calculate it later with SpatialCorrCalcRun.m!...')
    end
    
    % Data to keep accessable in every function
    data = struct();
    if exist('SPSIGfile','var')
        data.SPSIGfile = SPSIGfile; % filename to save data to
    end
    data.xas = (1:dim(1))./freq;
    data.dim = dim;
    data.freq = freq;
    data.imgStackT = imgStackT; % transpose spectral
    data.Sax = Sax;
    % BImg = the main background image. get only lowest 5th of frequencies
    data.BImg = squeeze(max(log(SPic),[],3))'; % Get the max (results in neurons popping out most)
    vals = unique(data.BImg);
    data.BImg(data.BImg<vals(2)) = vals(2);
    data.Mask = Mask;
    data.PP = PP;
    data.SpatialCorr = SpatialCorr;
    data.pos = []; % pixel positions of a ROI that is being splitted
    data.corners = []; % Corner coordinates from which to start correlations
    data.clIdx = []; % to which cluster each pixel belongs for ROI splitting
    data.clCon = []; % Contours of clustered ROIs
    data.newMask = []; % new mask after a ROI is created
    data.newRoi = []; % contour and ROI data for newly created ROI
    data.manRoi = struct('Con',struct('x',[],'y',[]),'mask',[],'A',0); % data of manually created ROI
    
    % legend for Roi Creation axes
    x = [0.8 0.2]; y = [0 1];
    lines = gobjects(1,4);
    h.createRoiLegend.NextPlot = 'add';
    h.createRoiLegend.Visible = 'off';
    lines(1) = plot(x,y,'xk','Parent',h.createRoiLegend);
    lines(2) = plot(x,y,'g','Parent',h.createRoiLegend);
    lines(3) = plot(x,y,'r','Parent',h.createRoiLegend);
    lines(4) = plot(x,y,'xr','Parent',h.createRoiLegend);
    legend(h.createRoiLegend,...
        {'clicked point', 'contour above threshold','chosen ROI',...
        'maximum of new ROI'},'Location','North');
    h.createRoiAx2.YDir = 'reverse';
    
    % Start with the ROI rejection turned on.
    h.majorToggle1 = TurnOn(h.majorToggle1);
    % Place major uiPanels where they should be: position of 1st panel
    h.majorPanel2.Position = h.majorPanel1.Position;
    h.majorPanel3.Position = h.majorPanel1.Position;
    h.majorPanel4.Position = h.majorPanel1.Position;
    h.majorPanel5.Position = h.majorPanel1.Position;
    h.majorPanel6.Position = h.majorPanel1.Position;
    h.majorPanel2.Visible = 'off'; % And make them invisable
    h.majorPanel3.Visible = 'off';
    h.majorPanel4.Visible = 'off';
    h.majorPanel5.Visible = 'off';
    h.majorPanel6.Visible = 'off';
    
    % Set the split signal axes to the same sizes as the normal signalAx1
    h.signalAx2.Position = h.signalAx1.Position;
    h.signalAx1.Visible  = 'off';
    
    h.clusterPlot1.XTick = [];  h.clusterPlot1.YTick = [];
    h.clusterPlot2.XTick = [];  h.clusterPlot2.YTick = [];
    
    % Choose default command line output for RoiRejecterGUI
    h.output = hObject;
    
    % Save the input data for access in later functions
    setappdata(hObject, 'sbxt', sbxt)
    setappdata(hObject, 'data', data)
    setappdata(hObject, 'idx', idx)
    setappdata(hObject, 'switches', switches)
    setappdata(hObject, 'colors', colors)
    
    % Start up the background image    data = getappdata(h.hGUI, 'data');
    h.mainAx.NextPlot = 'add';
    h.mainAx.YDir = 'reverse';
    h.im = imagesc(data.Mask, 'Parent',h.mainAx);
    
    % Check if background size is the same as sbx data size.
    if ~(dim(2)==size(data.BImg,2) && dim(3)==size(data.BImg,1))
        t = cell(1,4);
        t{1} = 'Size of Background image and TRANS datafile are not the same!';
        t{2} = 'Wrong files loaded';
        t{3}= sprintf('Trans datafile:_______%d x %d pixels',dim(2),dim(3));
        t{4} = sprintf('Background image:___%d x %d pixels',size(data.BImg,2),size(data.BImg,1));
        t{5} = sprintf('TRANS data file: %s', datafile);
        h.rejectInfo.String = {strjoin(t,'\n')};
        h.rejectInfo.ForegroundColor = [1 0 0];
        guidata(hObject, h);
        return
    end
    
    backGrdView(1, h);
    h.mainAx.XLim = [1 size(h.im.CData,2)];
    h.mainAx.YLim = [1 size(h.im.CData,1)];
    h = guidata(h.hGUI); % Update the handles, was edited by Backgroundview
    DrawRois(h);
    h = guidata(h.hGUI); % Update the handles, was edited by DrawRois
    UpdateRois(h);
    
    RejectInfoFunc(h);
    
    % Plot a datatrace, and create handles for that stuff
    numSelectSlider_Callback(h.numSelectSlider, 0, h)
    h = guidata(h.hGUI); % Update the handles, was edited by slider callback
    
    % Plot the 2Ps signal of the first neuron as an example and to fix that
    % YDir of the plot
    signal = mean(sbxt.Data.y(:,Mask==1),2);
    plot(data.xas, signal, 'k','Parent', h.signalAx2)
    h.signalAx2.YDir = 'normal';
    h.signalAx2.XLim = [data.xas(1), data.xas(end)];

    % Select 2P signal from mouse click surrounding selsize x selsize pixels
    % The GUI starts with selecting the data on the positon of ROI 1
    [xmesh, ymesh, signal] = SelectData(data.PP.P([1,2],1), switches.selSize, h, 'mean');
    % Plot the select data
    DrawData(xmesh, ymesh, signal, h)
    h = guidata(h.hGUI); % Update the handles, was edited by DrawRois
    
    % Update handles structure
    guidata(hObject, h);
end
% UIWAIT makes RoiRejecterGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RoiRejecterGUI_OutputFcn(hObject, eventdata, handles)
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    varargout{1} = handles.output;
end


% --- Executes during object creation, after setting all properties.
function hGUI_CreateFcn(hObject, eventdata, handles) %#ok
    % hObject    handle to hGUI (see GCBO)
    % handles    empty - handles not created until after all CreateFcns called
end


% --- Executes during object creation, after setting all properties.
function slider_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
    % hObject    handle to the slider calling at the moment
    % handles    empty - handles not created until after all CreateFcns called
    % Every slider calls this function for asthetic reasons
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end
end


% --- Executes during object creation, after setting all properties.
function backGrdView_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
    % hObject    handle to backGrdView (see GCBO)
    % handles    empty - handles not created until after all CreateFcns called

    hObject.String = {'Spectral', 'Spectral color', 'SpatialCorr', 'ROI corr',...
                      'Mask', 'Chronic (red) & spectral (green)',...
                      'Chronic', 'peak frequency', 'workspace variable'};
    
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

%% Button activation functions %%
% --- Executes on button press of major toggle buttons
function majorToggles_Callback(hObject, ~, h) %#ok<DEFNU>
    % Turns a functionality of the UI on or off. switches panels.

    buttonN = str2double(hObject.Tag(end)); % Which button was pressed? Important
    switches = getappdata(h.hGUI, 'switches');
    
    % If an inactive major toggle was clicked, change panels etc.
    if switches.currentMajor ~= buttonN
        % Turn previous panel off
        switch switches.currentMajor
            case 1 % ROI rejector panel
                h.majorToggle1 = TurnOff(h.majorToggle1);
                uistack(h.majorPanel1, 'bottom');
                h.majorPanel1.Visible = 'off';
            case 2 % ROI splitting panel
                h.majorToggle2 = TurnOff(h.majorToggle2);
                uistack(h.majorPanel2, 'bottom');
                h.majorPanel2.Visible = 'off';
            case 3 % ROI creation panel
                h.majorToggle3 = TurnOff(h.majorToggle3);
                uistack(h.majorPanel3, 'bottom');
                h.majorPanel3.Visible = 'off';
            case 4 % Show Data panel
                h.majorToggle4 = TurnOff(h.majorToggle4);
                uistack(h.majorPanel4, 'bottom');
                h.majorPanel4.Visible = 'off';
            case 5 % ROI Manual creation Panel
                h.majorToggle5 = TurnOff(h.majorToggle5);
                uistack(h.majorPanel5, 'bottom');
                h.majorPanel5.Visible = 'off';
        end
        
        % Turn requested panel on
        switch buttonN
            case 1 % ROI rejector panel
                h.majorToggle1 = TurnOn(h.majorToggle1);
                h.majorPanel1.Visible = 'on';
                uistack(h.majorPanel1, 'top');
            case 2 % ROI splitting panel
                h.majorToggle2 = TurnOn(h.majorToggle2);
                h.majorPanel2.Visible = 'on';
                uistack(h.majorPanel2, 'top');
            case 3 % ROI creation panel
                h.majorToggle3 = TurnOn(h.majorToggle3);
                h.majorPanel3.Visible = 'on';
                uistack(h.majorPanel3, 'top');
            case 4 % Show Data panel
                h.majorToggle4 = TurnOn(h.majorToggle4);
                h.majorPanel4.Visible = 'on';
                uistack(h.majorPanel4, 'top');
            case 5 % ROI manual creation panel
                h.majorToggle5 = TurnOn(h.majorToggle5);
                h.majorPanel5.Visible = 'on';
                uistack(h.majorPanel5, 'top');
        end
        
        if buttonN ~= 4 && switches.activeSignalAx == 1
            % If the sigselect signal axis was on, but the major toggle was
            % not set to 4, replace the visible signal axis with signalAx2
            h.signalAx1.Visible = 'off'; 
            uistack(h.signalAx1, 'bottom');
            h.signalAx2.Visible = 'on';
            uistack(h.signalAx2, 'top');
        end
        
        % Save new current panel
        switches.currentMajor = buttonN;
        setappdata(h.hGUI, 'switches', switches)
    end
end
    

function minorToggles_Callback(hObject, ~, h) %#ok<DEFNU>
    % Turns a functionality of the UI on or off
    % different buttons are:
    % corrButton
    % sigSelectButton
    % whiteButton
    % blackButton
    % addPointbutton
    % delClusterButton
    
    buttonName = hObject.Tag; % Which button was pressed?
    toggled = hObject.Value;
    switches = getappdata(h.hGUI, 'switches');
    
    % Turn the previous one that was on -> off
    if ~strcmp(switches.currentMinor, 'none')
        % Check for every possible button because we don't want to use eval
        if strcmp(switches.currentMinor, 'plotROIsButton')
            h.plotROIsButton = TurnOff(h.plotROIsButton);
        elseif strcmp(switches.currentMinor, 'corrButton')
            h.corrButton = TurnOff(h.corrButton);
        elseif strcmp(switches.currentMinor, 'sigSelectButton')
            h.sigSelectButton = TurnOff(h.sigSelectButton);
        elseif strcmp(switches.currentMinor, 'whiteButton')
            h.whiteButton = TurnOff(h.whiteButton);
        elseif strcmp(switches.currentMinor, 'blackButton')
            h.blackButton = TurnOff(h.blackButton);
        elseif strcmp(switches.currentMinor, 'addPointButton')
            h.addPointButton = TurnOff(h.addPointButton);
        elseif strcmp(switches.currentMinor, 'delClusterButton')
            h.delClusterButton = TurnOff(h.delClusterButton);
        else
            fprintf('\n\n\n\nsomething went terribly wrong!\n')
            fprintf('Tried to find a button called "%s" to disable, ',switches.currentMinor)
            fprintf('but no button was found with that name\n\n')
        end
    end
    
    % If the button is toggled on, save it as the current minor toggle
    if toggled
        hObject = TurnOn(hObject);
        switches.currentMinor = buttonName;
        
        % In case of sigSelect button, activate SignalAx1
        if strcmp(buttonName, 'sigSelectButton')
            switches.activeSignalAx = 1;
            % also turn the correct signal axis on
            h.signalAx1.Visible = 'on';
            uistack(h.signalAx1, 'top');
            h.signalAx2.Visible = 'off';
            uistack(h.signalAx2, 'bottom');
            
        elseif switches.activeSignalAx == 1
            % If active button is not sigSelectButton but active signal
            % Axes is 1, activate signalAx2
            switches.activeSignalAx = 2;
            % also turn the correct signal axis on
            h.signalAx1.Visible = 'off'; 
            uistack(h.signalAx1, 'bottom');
            h.signalAx2.Visible = 'on';
            uistack(h.signalAx2, 'top');
            
        end
    else
        switches.currentMinor = 'none';
    end
    % Save new current panel
    setappdata(h.hGUI, 'switches', switches)
end


% --- Executes on button press in plotListing.
function plotListing_Callback(hObject, ~, h) %#ok<DEFNU>
    % Save the state of the toggle/radio button
    % If turned on the signal of clicked ROIs while white and blacklisting
    % will be plotted. This can be slow sometimes
    switches = getappdata(h.hGUI, 'switches');
    switches.plotListing = hObject.Value;
    setappdata(h.hGUI, 'switches', switches)
end


function button = TurnOff(button)
    % Changes to toggle button when turned off
    button.Value = false;
    button.FontWeight = 'normal';
    button.ForegroundColor = [0 0 0];
end


function button = TurnOn(button)
    % Changes to toggle button when turned on
    button.Value = true;
    button.FontWeight = 'bold';
    button.ForegroundColor = [0 0 1];
end


% --- Executes on button press in thresCorToggle.
function thresToggle_Callback(~, ~, h)
    % Which of the two threshold rejection criteria to use?
    % MATLAB already turns the value of the clicked button automatically.
    % So this function makes sure one of the two is activated
    if h.thresToggle.Value && h.thresCorToggle.Value
        h.thresCorToggle.Value = false;
    else
        h.thresCorToggle.Value = true;
    end
    RejectInfoFunc(h)
    UpdateRois(h)
end


% --- Executes on button press in thresCorToggle.
function thresCorToggle_Callback(~, ~, h)
    % Which of the two threshold rejection criteria to use?
    if h.thresToggle.Value && h.thresCorToggle.Value
        h.thresToggle.Value = false;
    else
        h.thresToggle.Value = true;
    end
    RejectInfoFunc(h)
    UpdateRois(h)
end


%% Select ROIs for deletion functions %% 

function SelectROI(pos, h, do)
    % This function is called when either the blacklist, whitelist or plot
    % ROI button is toggled and there has been a click on the main axes.
    
    data = getappdata(h.hGUI, 'data');
    idx = getappdata(h.hGUI, 'idx');
    
    id = data.Mask(round(pos(2)), round(pos(1))); % Which ROI was clicked on?
    if id > 0
        % If the plotlisting is turned on plot the signal of this ROI
        switches = getappdata(h.hGUI, 'switches');
        if switches.plotListing || strcmp(switches.currentMinor, 'plotROIsButton')
            % Plot the signal of the ROI you just clicked
            sbxt = getappdata(h.hGUI, 'sbxt');
            signal = mean(sbxt.Data.y(:,data.Mask'==id),2);
            h.signalAx2.NextPlot = 'replacechildren';
            plot(data.xas, signal, 'color', [0 0 0], 'Parent', h.signalAx2);
            ylims = [round(min(signal),-2)-100, round(max(signal),-2)+100];
            xlims = h.signalAx2.XLim;
            h.signalAx2.YLim = ylims;
            % Say which ROI is being plotted, at 90% height of the signalAx 
            str = sprintf('ROI %d', id);
            text(xlims(1)+diff(xlims)/100, ylims(2)-diff(ylims)/20, str,...
                'Parent',h.signalAx2, 'Clipping', 'on')
        end
        
        switch do
            case 'black' % Add and remove ROIs to the black list to ensure deletion
                idx.Black(id) = mod(idx.Black(id)+1,2); % Toggle ROI to blacklist
                if idx.White(id) == 1 % Get out of whitelist
                    idx.White(id) = 0;
                end
            case 'white' % Add and remove ROIs to the white list to ensure survival
                idx.White(id) = mod(idx.White(id)+1,2); % Toggle ROI to whitelist
                if idx.Black(id) == 1 % Get out of blacklist
                    idx.Black(id) = 0;
                end                
        end
    end
    
    setappdata(h.hGUI, 'idx', idx)
    RejectInfoFunc(h)
    UpdateRois(h)
end


% --- Executes on slider movement.
function minSizeSlider_Callback(hObject, ~, h) %#ok<DEFNU>
    % Callback function to control minimum size limit
    switches = getappdata(h.hGUI, 'switches');
    switches.minSize = hObject.Value;
    switches.sliderSettings.minSize = hObject.Value;
    setappdata(h.hGUI, 'switches', switches)
    RejectRoiSizes(h)
    RejectInfoFunc(h)
    UpdateRois(h)
    
    % Update slider title text
    h.minSizeTitle.String = {sprintf('Min size: %4.0f px', hObject.Value)};
end


% --- Executes on slider movement.
function maxSizeSlider_Callback(hObject, ~, h) %#ok<DEFNU>
    % Callback function to control minimum size limit
    switches = getappdata(h.hGUI, 'switches');
    switches.maxSize = hObject.Value;
    switches.sliderSettings.maxSize = hObject.Value;
    setappdata(h.hGUI, 'switches', switches)
    RejectRoiSizes(h)
    RejectInfoFunc(h)
    UpdateRois(h)
    
    % Update slider title text
    h.maxSizeTitle.String = {sprintf('Max size: %4.0f px', hObject.Value)};
end


function RejectRoiSizes(h)
    % Reject ROIs based on minimum and maximum sizes.
    
    switches = getappdata(h.hGUI, 'switches');
    data = getappdata(h.hGUI, 'data');
    idx = getappdata(h.hGUI, 'idx');
    % Calculates ROI indexes to delete, based on size limits
    idx.Size = (data.PP.A <switches.minSize | data.PP.A > switches.maxSize);
    setappdata(h.hGUI, 'idx', idx)
end


% --- Executes on slider movement.
function thresSlider_Callback(hObject, ~, h) %#ok<DEFNU>
    % changes how many ROIs are rejected based on creation threshold
    threshold = hObject.Value;
    h.thresTitle.String = sprintf('Threshold: %.2f',threshold);
    
    switches = getappdata(h.hGUI, 'switches');
    idx = getappdata(h.hGUI, 'idx');
    data = getappdata(h.hGUI, 'data');
    PP = data.PP;
    switches.sliderSettings.threshold = hObject.Value;
    
    if threshold > 0
        % Calculate ROI indexes to delete based on given threshold and data
        Df = PP.P(3,:) - PP.P(4,:); % max ? - threshold
        RngD = Df - min(Df);
        RngD = RngD./max(RngD);
        idx.Thres = (RngD < threshold);
    else
        idx.Thres = false(1,PP.Cnt);
    end
    
    setappdata(h.hGUI, 'switches', switches)
    setappdata(h.hGUI, 'idx', idx) % save new indexes
    RejectInfoFunc(h)
    UpdateRois(h) % apply coloring to ROIs
end


% --- Executes on slider movement.
function thresCorSlider_Callback(hObject, ~, h) %#ok<DEFNU>
    % Changes how ROI are rejected based on 'mean inner correlation
    % coefficients of raw signals' threshold
    threshold = hObject.Value;
    h.thresCorTitle.String = sprintf('Mean inner corr: %.2f',threshold);
    
    switches = getappdata(h.hGUI, 'switches');
    idx = getappdata(h.hGUI, 'idx');
    data = getappdata(h.hGUI, 'data');
    PP = data.PP;
    switches.sliderSettings.innerCorr = hObject.Value;
    
    if size(PP.P,1) == 5 % Check if the inner correlation values are in the data
        if threshold > 0
            idx.ThresCor = (PP.P(5,:) < threshold);
        else
            idx.ThresCor = false(1,PP.Cnt);
        end
    else
        h.thresCorTitle.String = 'required data missing!';
        h.thresCorTitle.ForegroundColor = [1 0 0];
    end
    
    setappdata(h.hGUI, 'switches', switches)
    setappdata(h.hGUI, 'idx', idx) % save new indexes
    RejectInfoFunc(h)
    UpdateRois(h) % apply coloring to ROIs
end



% --- Executes on slider movement.
function roundnessSlider_Callback(hObject, ~, h) %#ok<DEFNU>
	% Remove ROIs based on their roundness
    threshold = hObject.Value;
    h.roundnessSliderTitle.String = sprintf('Min roundness: %.2f',threshold);
    
    switches = getappdata(h.hGUI, 'switches');
    data = getappdata(h.hGUI, 'data');
    idx = getappdata(h.hGUI, 'idx');
    
    switches.sliderSettings.roundedness = hObject.Value;
    
    % Calculate the roundness of all ROIs
    roundness = zeros(1,data.PP.Cnt);
    for i = 1:data.PP.Cnt
        roundness(i) = perimarea(data.PP.Con(i).x, data.PP.Con(i).y);
    end
    
    idx.Round = roundness < threshold; % true for not round enough ROIs
    
    setappdata(h.hGUI, 'switches', switches)
    setappdata(h.hGUI, 'idx', idx) % save new indexes
    RejectInfoFunc(h)
    UpdateRois(h) % apply coloring to ROIs
end


function RejectInfoFunc(h)
    % update the text that tells how many ROIs will be rejected for what
    % reasons
    idx = getappdata(h.hGUI, 'idx');
    data = getappdata(h.hGUI, 'data');
    
    t = cell(7,1);
    if h.thresToggle.Value % The creation threshold should be used
        selectedPoints = sum((idx.Black|idx.Size|idx.Thres) & ~idx.White);
        t{2} = sprintf('%3d ROIs rejected because threshold',sum(idx.Thres));
    else % The correlation threshold should be used
        selectedPoints = sum((idx.Black|idx.Size|idx.ThresCor) & ~idx.White);
        t{2} = sprintf('%3d ROIs rejected because correlation threshold',sum(idx.ThresCor));
    end
    
    t{1} = sprintf('%3d ROIs rejected because size limits',sum(idx.Size));
    t{3} = sprintf('%3d ROIs rejected because of roundness threshold', sum(idx.Round));
    t{4} = sprintf('%3d ROIs manually rejected',sum(idx.Black));
    t{5} = sprintf('%3d ROIs white listed',sum(idx.White));
    t{6} = '_________________';
    t{7} = sprintf('%3d of %3d ROIs will be deleted:\n %3d remain',...
        selectedPoints, data.PP.Cnt, (data.PP.Cnt-selectedPoints));
    summary = strjoin(t,'\n');
    h.rejectInfo.String = {summary};
end


%% Edit ROIs based on selection
% --- Executes on button press in deleteButton.
function deleteButton_Callback(~, ~, h) %#ok<DEFNU>
    % Callback function for the apply delete button
    % Will remove ROIs from data
    idx = getappdata(h.hGUI, 'idx');
    data = getappdata(h.hGUI, 'data');
    PP = data.PP;
    Mask = data.Mask;
    
    if h.thresToggle.Value
        selectedPoints = (idx.Black|idx.Size|idx.Round|idx.Thres) & ~idx.White;
    else
        selectedPoints = (idx.Black|idx.Size|idx.Round|idx.ThresCor) & ~idx.White;
    end
    
    if ~isempty(find(selectedPoints, 1))
        idxDel = logical(selectedPoints);
        PP.Con(idxDel) = [];
        PP.A(idxDel) = [];
        PP.P(:,idxDel) = [];
        PP.SpecProfile(:,idxDel) = [];
        PP.peakFreq(idxDel) = [];
        PP.creationMethod(idxDel) = [];
        PP.Cnt = size(PP.P,2);
        
        % Update the selected points list
        idx.Black = false(1, PP.Cnt);
        idx.Size = false(1, PP.Cnt);
        idx.Thres = false(1, PP.Cnt);
        idx.ThresCor = false(1, PP.Cnt);
        idx.Round = false(1, PP.Cnt);
        
        % Clean up the mask
        for i = find(idxDel)
            Mask(Mask==i) = 0;
        end
        idx.White(idxDel) = [];
        
        v = unique(Mask(:));
        v = v(2:end);        
        for i = 1:length(v)
            if v(i) ~= i
                Mask(Mask == v(i)) = i;
            end
        end
        
        % Save edited variables
        data.PP = PP;
        data.Mask = Mask;
        setappdata(h.hGUI, 'data', data)
        setappdata(h.hGUI, 'idx',  idx) 
        
        UpdateMainImg(h)
    end
end
    

%% ROI creation functions
% --- Executes on slider movement.
function creationThresholdSlider_Callback(hObject, ~, h) %#ok<DEFNU>
    % hObject    handle to creationThresholdSlider (see GCBO)
    % handles    structure with handles and user data (see GUIDATA)
    switches = getappdata(h.hGUI, 'switches');
    switches.creationThres = round(hObject.Value);
    hObject.Value = switches.creationThres;
    h.creationThresholdTitle.String = sprintf('Threshold: %d%%', hObject.Value);
    setappdata(h.hGUI, 'switches', switches); % update the switches
    CreateRoi(h,'update')
end


function CreateRoi(h, do)
    % Creating new ROIs, using the spectral image and a threshold set by
    % the user
    data = getappdata(h.hGUI, 'data');
    switches = getappdata(h.hGUI, 'switches');
    MaskTemp = data.Mask; % temporary mask that gets edited
    pos = switches.creationPos;
    
    if MaskTemp(pos(2),pos(1))==0
        voxelSz = round(64./2); % SIZE OF THE SEARCH AREA
        % Make sure the searchfield isn't outside of the image
        piecex = pos(2)-voxelSz:pos(2)+voxelSz; % x coordinates of the patch of data
        piecey = pos(1)-voxelSz:pos(1)+voxelSz; % y coordinates of the patch of data

        piecex(piecex<1) = [];
        piecey(piecey<1) = [];
        piecex(piecex>data.dim(3)) = [];
        piecey(piecey>data.dim(2)) = [];
        pieceWidth = length(piecex);
        pieceHeight = length(piecey);
        
        % get the coordinate of the clicked point. Coordinate is not the
        % center of the searchfield when the searchfield is cut because it
        % is close to the edges
        centre = [pos(2)-piecex(1), pos(1)-piecey(1)];
        
        img = h.im.CData;
        dims = size(img);
        % if the current background image is a color image average colors
        if length(dims) == 3
            img = mean(img, 3);
        end
        vals = sort(unique(img(:))); 
        
        % A new place for a new ROI is being requested: create and plot search area
        if strcmp(do,'new')
            plot(pos(1),pos(2), 'x', 'color', 'w', 'hittest', 'off', 'Parent', h.mainAx);
            
            % expand ROI masks so no ROIs can be made close to existin ROIs
            [~,~,MaskTemp]= BufferMask(data.Mask, 2);
            MaskTemp = MaskTemp + data.Mask; 
            img(MaskTemp>0) = vals(2); % remove pieces that are already ROIs
            
            % Get a 'voxel' of the data and filter the search space to reduce noise
            piece = imgaussfilt(img(piecex, piecey),1.3);

            % Plot the 'voxel'/ searchfield with the spectral data
            h.createRoiAx2.NextPlot = 'replacechildren';
            imagesc(piece, 'Parent', h.createRoiAx2)
            h.createRoiAx2.NextPlot = 'add';
            h.createRoiAx2.XLim = [1, size(piece,2)];
            h.createRoiAx2.YLim = [1, size(piece,1)];
            plot(centre(2),centre(1), 'xk', 'Parent', h.createRoiAx2)
        elseif strcmp(do,'update')
            % A different contour threshold is being requested. Try to
            % extract the search area that was calculated previously
            piece = getimage(h.createRoiAx2);
            
            if isempty(piece) % If no searchfield had been created
                fprintf('changed threshold while not creating ROI\n')
                return
            else
                h.createRoiAx2.NextPlot = 'replacechildren';
                imagesc(piece, 'Parent', h.createRoiAx2)
                h.createRoiAx2.NextPlot = 'add';
                plot(centre(2),centre(1), 'xk', 'Parent', h.createRoiAx2)
            end
        end
        
        plot(centre(2),centre(1), 'xk', 'Parent', h.createRoiAx2)
        found = false;
        
        threshval = switches.creationThres;
        threshold = prctile(piece(:), threshval); % calculate percentile threshold
        con = contourc(piece, [threshold threshold]); % Get the contours
        con = getcontourlines(con);
        
        % Only use the contour if it has the centre point inside of itself
        for i = 1:length(con) % Check every contour
            MaskTemp = poly2mask(con(i).x, con(i).y, pieceWidth,  pieceHeight);
            plot(con(i).x, con(i).y, 'color', [0 1 0 0.5], 'Parent', h.createRoiAx2)
            if MaskTemp(centre(1),centre(2))
                 % the contour that includes the clicked point is found
                num = i;
                found = true;
                plot(con(i).x, con(i).y, 'color', 'r', 'Parent', h.createRoiAx2)
            end
        end
        
        if found
            con = con(num);
            switches.creationAllow = true; % This could be a good ROI
            sbxt = getappdata(h.hGUI, 'sbxt');
            
            % Create a mask that includes this new ROI
            MaskTemp = poly2mask(con.x, con.y, pieceWidth,  pieceHeight);
            MaskNew = data.Mask;
            MaskNew(MaskNew>0) = MaskNew(MaskNew>0) + 1; % make room for the new ROI
            MaskNew(piecex, piecey) = double(MaskTemp) + MaskNew(piecex, piecey);
            
            % Calculate the extra info
            piece(MaskTemp==0) = min(piece(:));
            [maxval, y] = max(piece); % maximum value for 1st dimension
            [maxval, x] = max(maxval); % maximum value of complete matrix
            y = y(x); % so which row was it again? - the xth row
            xcord = x+piecey(1)-1; % the coordinates in the complete mask
            ycord = y+piecex(1)-1;
            
            % Plot the local maximum
            plot(x,y, 'xr', 'Parent', h.createRoiAx2)
            
            % Create the contour data for this new ROI
            newRoi = struct;
            newRoi.Con.x = con.x + piecey(1) - 1;
            newRoi.Con.y = con.y + piecex(1) - 1;
            newRoi.A = sum(MaskTemp(:));
            newRoi.P = [x+piecey(1)-1; y+piecex(1)-1; threshold; maxval];
            
            % If the correlation values are saved (in PP.P(5,:)) calculate
            % for this new ROI as well
            if size(data.PP.P, 1) == 5
                signal = sbxt.Data.y(:, MaskNew'==1);
                shorten = floor(data.freq);
                signaldecim = zeros(ceil(size(signal,1)/shorten), size(signal,2));
                cord1D = sub2ind(size(data.Mask'),xcord, ycord);
                signalCenter = mean(sbxt.Data.y(:, cord1D-2:cord1D+2),2);
                % signals for all pixels in ROI decimated by freq
                for i = 1:size(signal,2) 
                    signaldecim(:,i) = decimate(double(signal(:,i)), shorten);
                end
                signalCenterDecim = decimate(double(signalCenter), shorten);
                corVal = mean(corr(signalCenterDecim,signaldecim,'rows','pairwise'));
                newRoi.P = [newRoi.P; corVal];
            end
            
            % Plot the signal of the found ROI
            h.signalAx2.NextPlot = 'replacechildren';
            if size(data.PP.P,1) == 5
                str = cell(5,1);
                % If the correlation with the max was calculated, plot both
                plot(data.xas, signalCenter, 'color', [0 0 0.5 0.2], 'Parent', h.signalAx2)
                h.signalAx2.NextPlot = 'add';
                plot(data.xas, mean(signal,2), 'color', [0.5 0 0], 'Parent', h.signalAx2)
                legend(h.signalAx2, {'signal of the brightest pixels','signal of the new ROI'})
                h.signalAx2.YLim = [round(min(signaldecim(:)),-2)-100, round(max(signalCenter),-2)+100];
                str{5} = sprintf('average signal correlation=%.2f',corVal);
            else
                str = cell(4,1);
                signal = mean(sbxt.Data.y(:, MaskNew'==1),2);
                plot(data.xas, signal, 'color', [0.5 0 0], 'Parent', h.signalAx2)
                legend(h.signalAx2, {'signal of the new ROI'})
                h.signalAx2.YLim = [round(min(signal),-2)-100, round(max(signal),-2)+100];
                h.signalAx2.YLim = [round(min(signal),-2)-100, round(max(signal),-2)+100];
            end
            h.signalAx2.YDir = 'normal';
            % Info about the new neuron
            str{1} = sprintf('threshold value: %.3f', threshold);
            str{2} = sprintf('maximum value:  %.3f', maxval);
            % Do not allow any sized ROI to be saved. ROI needs to be
            % bigger then 35 pixels
            if newRoi.A < 10
                str{3} = sprintf('size: %d px <- VERY LOW!', newRoi.A);
                switches.creationAllow = false; % do not allow this ROI to be saved
            elseif newRoi.A < 20
                str{3} = sprintf('size: %d px <- LOW!', newRoi.A);
                switches.creationAllow = false;
            else
                str{3} = sprintf('size: %d px', newRoi.A);
            end
            roundedness = perimarea(con.x, con.y);
            if roundedness<0.4
                str{4} = sprintf('roundedness=%.2f <- LOW!', roundedness);
            else
                str{4} = sprintf('roundedness=%.2f', roundedness);
            end
            h.creationRoiText.String =  strjoin(str, '\n');
            h.creationRoiText.FontWeight =  'normal';
            
            % Saving new neuron to data, where it will wait for
            % confirmation to be applied to actual dataset by the user
            data.newRoi = newRoi;
            data.newMask = MaskNew;
            setappdata(h.hGUI, 'data', data)
            
        else % No contour found around the clicked point
            str = cell(4,1);
            str{1} = sprintf('threshold value: %.3f', threshold);
            str{2} = 'did not find a contour around click!';
            str{3} = 'try changing the threshold';
            str{4} = 'or click more accurately on a neuron';
            h.creationRoiText.String =  strjoin(str, '\n');
            h.creationRoiText.FontWeight =  'normal';
            
            switches.creationAllow = false; % After this fail, do not allow creation of ROI
        end
        
        h.createRoiAx1.NextPlot = 'replacechildren';
        bins = 50; % number of bins
        piece = img(piecex, piecey); % Get a 'voxel' of the data
        piece(piece<=vals(2)) = NaN;
        histogram(piece, bins,'EdgeAlpha',0, 'Parent', h.createRoiAx1)
        h.createRoiAx1.NextPlot = 'add';
        h.createRoiAx1.YLim = h.createRoiAx1.YLim;% Set ylimits fixed
        line([threshold, threshold], h.createRoiAx1.YLim, 'color','k',...
            'linewidth', 2,'Parent', h.createRoiAx1)
        legend(h.createRoiAx1, {'spectral density values', 'threshold'},...
            'Location','NorthOutside','FontSize',8)
        
    else % User clicked on an existing ROI
        str = cell(2,1);
        str{1} = '';
        str{2} = 'There is already a ROI there!';
        str{3} = 'try not clicking on existing ROIs';
        h.creationRoiText.String =  strjoin(str, '\n');
        h.creationRoiText.FontWeight =  'bold';
        
        switches.creationAllow = false; % After this fail, do not allow creation of ROI
    end
    
    % Change color of the apply button based on whether saving is allowed
    if switches.creationAllow
        h.applyNewRoi.ForegroundColor = [0 1 0]; % green if allowed to save
    else
        h.applyNewRoi.ForegroundColor = [1 0 0]; % red if not allowed to save
    end
    setappdata(h.hGUI, 'switches', switches)
end


% --- Executes on button press in applyNewRoi.
function applyNewRoi_Callback(~, ~, h) %#ok<DEFNU>
    % Apply the new ROI
    switches = getappdata(h.hGUI, 'switches');
    data = getappdata(h.hGUI, 'data');
    idx = getappdata(h.hGUI, 'idx');
    
    if switches.creationAllow
        
        h.applyNewRoi = TurnOn(h.applyNewRoi);
        h.applyNewRoi.String = {'Applying new ROI'};
        % Update indexes
        idx.Size  = [false,idx.Size];
        idx.Thres = [false,idx.Thres];
        idx.ThresCor = [false,idx.ThresCor];
        idx.Round = [false, idx.Round];
        idx.Black = [false,idx.Black];
        idx.White = [true, idx.White]; % immediately whitelist

        % Update the data
        data.Mask = data.newMask;
        data.PP.A = [data.newRoi.A, data.PP.A];
        data.PP.P = [data.newRoi.P, data.PP.P];
        data.PP.Con = [data.newRoi.Con, data.PP.Con];
        data.PP.creationMethod = [{'manual click'}; data.PP.creationMethod];
        data.PP.Cnt = data.PP.Cnt + 1;
        
        % Update the spectral profiles and peak frequency
        nSpec = size(data.imgStackT,3);
        ROIi = data.Mask == 1;
        ROIi3D = repmat(ROIi, [1,1,nSpec]);
        specProfile = reshape(data.imgStackT(ROIi3D), [sum(ROIi(:)), nSpec]);
        data.PP.SpecProfile = [mean(log(specProfile))' data.PP.SpecProfile];
        [~, peakFreq] = max(data.PP.SpecProfile);
        data.PP.peakFreq = data.Sax(peakFreq);
        
        switches.creationAllow = false;
        
        setappdata(h.hGUI, 'data', data);
        setappdata(h.hGUI, 'switches', switches);
        setappdata(h.hGUI, 'idx', idx);

        UpdateMainImg(h)
        
        try
            pause(2)
            h.applyNewRoi = TurnOff(h.applyNewRoi);
            h.applyNewRoi.String = {'Apply'};
        catch
            fprintf('Figure got closed quickly after creation a new ROI\n')
        end
    else % If applying the ROI is not allowed tell the user
        
        h.applyNewRoi.String = {'bad ROI, not saving'};
        try
            pause(2)
            h.applyNewRoi.String = {'Apply'};
        catch
            fprintf('Figure got closed suddenly\n')
        end
    end
end

%% ROI manual creation functions

function ManualRoi(pos, button, h)
    % Add a new part of the contours
    data = getappdata(h.hGUI, 'data');
    switches = getappdata(h.hGUI, 'switches');
    
    nCords = length(data.manRoi.Con.x); % number of coordinates of contour
    
    if button == 1 % left mouse button: add position to the contour
        data.manRoi.Con.x(end+1) = pos(1);
        data.manRoi.Con.y(end+1) = pos(2);        
        nCords = nCords + 1;
    else % other mouse button: remove last positition of contour
        if nCords > 0
            data.manRoi.Con.x(end) = [];
            data.manRoi.Con.y(end) = [];
        nCords = nCords - 1;
        end
    end
    
    switches.creationManAllow = false;
    
    % Create a mask that includes this new ROI
    if nCords > 2
        x = data.manRoi.Con.x;
        x = [x x(1)]; % close the contour
        y = data.manRoi.Con.y;
        y = [y y(1)]; % close the contour
        
        delete(h.manRoi)
        h.manRoi = plot(x, y, '-co', 'Parent',h.mainAx, 'Hittest','off');
        
        manMask = poly2mask(x, y, data.dim(3), data.dim(2));
        
        npixels = sum(manMask(:) > 0);
      
        data.manRoi.Mask = manMask;
        data.manRoi.A = npixels;
        
        if npixels>0 
            
            % find coordinates and maximum value inside the new ROI
            imgTemp = data.BImg;
            imgTemp(manMask==0) = min(data.BImg(:)); 
            [maxVal, maxy] = max(imgTemp);
            [maxVal, maxx] = max(maxVal);
            maxy = maxy(maxx);
            [maxCord] = sub2ind(size(manMask'), maxx, maxy);
        
            [~,~,contourPixels] = BufferMask(double(manMask),2);
            
            data.manRoi.P = zeros(4,1);
            data.manRoi.P(1) = maxx; % x coordinate of brightest point
            data.manRoi.P(2) = maxy; % y coordinate of brightest point
            data.manRoi.P(3) = mean(data.BImg(contourPixels==1));
            data.manRoi.P(4) = maxVal;  % maximum value of BImg in ROI
            str = cell(2,1); % Print text for new manual mask
            
            % Plot image of new manual mask
            h.manRoiAx.NextPlot = 'replacechildren';
            im = CreateRGB({data.BImg, manMask, contourPixels}, 'g r b'); % image of new contour over the BImg
            imagesc(im, 'Parent',h.manRoiAx, 'hittest','off');
            h.manRoiAx.NextPlot = 'add';
            plot(maxx, maxy, 'wx', 'Parent', h.manRoiAx)
            h.manRoiAx.XLim = [min(x)-20, max(x)+20];
            h.manRoiAx.YLim = [min(y)-20, max(y)+20];            
            h.manRoiAx.YDir = 'reverse';
            
            % If the correlation values are saved (in PP.P(5,:)) calculate
            % for this new ROI as well
            sbxt = getappdata(h.hGUI, 'sbxt');
            signal = sbxt.Data.y(:, manMask');
            signaldecim = zeros(ceil(size(signal,1)/floor(data.freq)), size(signal,2));
            signalCenter = mean(sbxt.Data.y(:, maxCord-2:maxCord+2),2);
            %signals for all pixels in ROI decimated by freq
            for i = 1:size(signal,2) 
                signaldecim(:,i) = decimate(double(signal(:,i)), floor(data.freq));
            end
            signalCenter = decimate(double(signalCenter), floor(data.freq));

            % Plot the signal of the found ROI
            h.signalAx2.NextPlot = 'replacechildren';
            plot(signalCenter, 'color', [0 0 0.5 0.4], 'Parent', h.signalAx2)
            h.signalAx2.NextPlot = 'add';
            plot(mean(signaldecim,2), 'color', [0.5 0 0], 'Parent', h.signalAx2)
            legend(h.signalAx2, {'signal of the brightest pixels','signal of the new ROI'})
            h.signalAx2.YLim = [round(min([mean(signal,2); signalCenter]),-2)-100,...
                                round(max([mean(signal,2); signalCenter]),-2)+100];
            h.signalAx2.YDir = 'normal';
            
            if npixels > 35 % MINIMUM NUMBER OF PIXELS FOR ROI: ALLOW CREATION
                switches.creationManAllow = true;
                str{1} = sprintf('contour size = %d pixels', npixels);
            else
                str{1} = sprintf('contour size = %d pixels <- Too small to apply', npixels);
            end
            roundedness = perimarea(x, y);
            if roundedness<0.4
                str{2} = sprintf('roundedness=%.2f <- LOW!', roundedness);
            else
                str{2} = sprintf('roundedness=%.2f', roundedness);
            end
            str{3} = sprintf('%d points for contour', nCords);
            str{4} = 'red = new mask';
            str{5} = 'green = Background image';
            str{6} = 'blue = contourpixels';
            
            if size(data.PP.P, 1)==5 % or inner correlation value
                corVal = mean(corr(signalCenter,signaldecim,'rows','pairwise'));
                data.manRoi.P = [data.manRoi.P; corVal];
                str = [str(1:3); sprintf('mean inner correlation: %.3f', corVal); str(4:end) ];
            end
            h.manRoiText.String =  strjoin(str, '\n');
            
        else
            h.manRoiText.String = 'zero pixels inside contour';
        end
    else % less then 2 coordinates are clicked. no ROI possible yet
        h.manRoiAx.NextPlot = 'replacechildren';
        im = zeros(2,2,3); % black image
        imagesc(im, 'Parent',h.manRoiAx);
        h.manRoiAx.XLim = [1 2];
        h.manRoiAx.YLim = [1 2];
        
        h.manRoiText.String = sprintf('Too few points to create ROI: %d points', nCords);
        
        delete(h.manRoi)
        h.manRoi = plot(data.manRoi.Con.x, data.manRoi.Con.y, '-co',...
            'Parent',h.mainAx, 'Hittest','off');
    end
    
    guidata(h.hGUI,h); % Save the handles
    setappdata(h.hGUI, 'data', data);
    setappdata(h.hGUI, 'switches', switches);
end

% --- Executes on button press in applyManRoi.
function applyManRoi_Callback(~, ~, h) %#ok<DEFNU>
    switches = getappdata(h.hGUI, 'switches');
    data = getappdata(h.hGUI, 'data');
    idx = getappdata(h.hGUI, 'idx');
    
    if switches.creationManAllow
        % turn on the button as a response
        h.applyManRoi = TurnOn(h.applyManRoi);
        h.applyManRoi.String = {'Applying new ROI'};
        
        % Check if there are no gaps in between ROI numbers
        bad = find(diff(unique(data.Mask))>1);
        if ~isempty(bad)
            warning('ROI numbers have gaps! BAD!!! (applyManRoi_Callback)')
        end
        
        % Update indexes
        idx.Size  = [false,idx.Size];
        idx.Thres = [false,idx.Thres];
        idx.ThresCor = [false,idx.ThresCor];
        idx.Round = [false, idx.Round];
        idx.Black = [false,idx.Black];
        idx.White = [true, idx.White]; % immediately whitelist

        % close the contour
        data.manRoi.Con.x = [data.manRoi.Con.x, data.manRoi.Con.x(1)];
        data.manRoi.Con.y = [data.manRoi.Con.y, data.manRoi.Con.y(1)];
        
        % Update the data
        data.Mask(data.Mask>0) = data.Mask(data.Mask>0) + 1;
        data.Mask(data.manRoi.Mask) = 1;
        data.PP.A = [data.manRoi.A, data.PP.A];
        data.PP.P = [data.manRoi.P, data.PP.P];
        data.PP.Con = [data.manRoi.Con, data.PP.Con];
        data.creationMethod = [{'manual draw'}; data.PP.creationMethod];
        data.PP.Cnt = data.PP.Cnt + 1;
        
        % Update the spectral profiles and peak frequency
        nSpec = size(data.imgStackT,3);
        ROIi = data.Mask == 1;
        ROIi3D = repmat(ROIi, [1,1,nSpec]);
        specProfile = reshape(data.imgStackT(ROIi3D), [sum(ROIi(:)), nSpec]);
        data.PP.SpecProfile = [mean(log(specProfile))' data.PP.SpecProfile];
        [~, peakFreq] = max(data.PP.SpecProfile);
        data.PP.peakFreq = data.Sax(peakFreq);
        
        bad = find(diff(unique(data.Mask))>1);
        if isempty(bad) % Save the updated handles

            % reset manroi data
            data.manRoi = struct('Con',struct('x',[],'y',[]),'mask',[],'A',0);
            switches.creationManAllow = false;

            setappdata(h.hGUI, 'data', data);
            setappdata(h.hGUI, 'switches', switches);
            setappdata(h.hGUI, 'idx', idx);

            UpdateMainImg(h)
        else
            str = sprintf('NOT APPLYING manual ROI!!! %d ROIs would have been deleted!', length(bad));
            h.manRoiText.String = str;
            h.applyManRoi.ForegroundColor = [1 0 0];
            h.applyManRoi.String = 'NO, I will not overwrite ROIs';
        end
        
        
        % reset text of button
        try
            pause(2)
            h.applyManRoi = TurnOff(h.applyManRoi);
            h.applyManRoi.String = {'Apply'};
        catch
            fprintf('Figure got closed quickly after creation a new ROI\n')
        end
    else
        fprintf('not allowed to save manually created ROI atm.\n')
    end
end


%% Roi Clustering functions %%
% --- Executes on slider movement.
function RoiSplitter(pos, h)
    % Select a ROI for splitting and stuff
    switches = getappdata(h.hGUI, 'switches');
    data = getappdata(h.hGUI, 'data');
    
    x = round(pos(1));
    y = round(pos(2));
    
    i = data.Mask(y, x);
    if i > 0
        roi = i;
        try
            [data.corners, ~] = GetRoiCorners(data.Mask, data.PP, 2, roi);
        catch
            fprintf('Had to take no buffer around mask for roi corners.\n')
            [data.corners, ~] = GetRoiCorners(data.Mask, data.PP, 0, roi);
        end
        
        switches.splitStarted = true;
        switches.splitRoi = roi;
        setappdata(h.hGUI, 'switches', switches)
        setappdata(h.hGUI, 'data', data)
        ClusterRoi(h)
    end
end


function ClusterRoi(h)
    % This function does the actual ROI splitting 
    % Get correlation of corners with other signals
    data = getappdata(h.hGUI, 'data');
    switches = getappdata(h.hGUI, 'switches');
    sbxt = getappdata(h.hGUI, 'sbxt');
    
    [corrImg, ~, corrs, corrsSub, pos] =...
        CorrRoiCorners(data.corners, switches.splitRoi, data.Mask, sbxt, data.dim);
    
    % Clustering with k-means
    nCluster = switches.nCluster;
    % normalization of coordinates
    info1 = [pos.x./(std(pos.x)), pos.y./(std(pos.y))]; % pixel coordinates
    info2 = corrs{1}./mean(std(corrs{1})); % pixel correlation with 4 start points
    info3 = corrsSub{1}./mean(std(corrsSub{1}));
    clIdx = kmeans([info1./2 info2 info3], nCluster); % splitting into clusters
    
    xpos = [min(pos.x)-4, max(pos.x)+4];
    ypos = [pos.y(1)-4, pos.y(end)+4];
    
    % correlation image
    strC = 'r g b gb';
    newMaskc = squeeze(squeeze(num2cell(corrImg, [1 2])));
    corrImg = permute(CreateRGB(newMaskc(1:4), strC), [2 1 3]);
    h.clusterPlot1.NextPlot = 'replacechildren';
    imagesc(corrImg, 'Parent',h.clusterPlot1, 'hittest','off');
    h.clusterPlot1.XLim = xpos;
    h.clusterPlot1.YLim = ypos;
    h.clusterPlot1.YDir = 'reverse';
    
    RGB = cell(1, 5);
    C = [];
    for k = 1:nCluster
        RGB{k} = zeros(size(newMaskc{1}));
        RGB{k}(pos.d1(clIdx==k)) = 1;
        
        % Calculate contour around each cluster
        Ctemp = contourc(RGB{k}, [1 1]);
        Ctemp = getcontourlines(Ctemp);
        maxidx = 1;
        if length(Ctemp)>1 % Only keep the largest contour.
            maxsize = 1;
            for i = 1:length(Ctemp)
                if maxsize<length(Ctemp(i).x)
                    maxsize = length(Ctemp(i).x);
                    maxidx = i;
                end
            end
        end
        C = [C, Ctemp(maxidx)];
        
        % Use contour to reject small island clusters
        in = inpolygon(pos.x, pos.y, Ctemp(maxidx).y, Ctemp(maxidx).x);
        RGB{k}(pos.d1(~in & clIdx==k)) = 0.5; % darken the image
        clIdx(~in & clIdx==k) = 0;
    end
    
    % Plot Mask of clusters
    [RGB, cval] = CreateRGB(RGB, [strC, ' rg']);
    
    h.clusterPlot2.NextPlot = 'replacechildren';
    imagesc(permute(RGB,[2 1 3]),'Parent',h.clusterPlot2, 'hittest','off');
    h.clusterPlot2.NextPlot = 'add';
    for i = 1:nCluster
        plot(C(i).y, C(i).x, 'color', 'w','Parent',h.clusterPlot2, 'hittest','off');
    end
    h.clusterPlot2.XLim = xpos;
    h.clusterPlot2.YLim = ypos;
    h.clusterPlot2.YDir = 'reverse';
    
    % plot fluoresence signals of the clusters
    clusterSignal = zeros(nCluster, length(data.xas));
    for i = 1:nCluster
        clusterSignal(i,:) = mean(sbxt.Data.y(:, pos.d1(clIdx==i)),2);
    end
    h.signalAx2.NextPlot = 'replacechildren';
    h.clusterSignal = plot(data.xas, clusterSignal, 'parent', h.signalAx2);
    for i = 1:nCluster
        h.clusterSignal(i).Color = cval(i,:);
    end
    h.signalAx2.YLim = [min(2000, round(min(clusterSignal(:)),-2)-100),...
                        max(6000, round(max(clusterSignal(:)),-2)+100)];
    h.signalAx2.XLim = [1, data.xas(end)];
    h.signalAx2.YDir = 'normal';
    
    data.pos = pos.d1;
    data.clIdx = clIdx;
    data.clCon = C;
    setappdata(h.hGUI, 'data', data)
end


function AddPoint(pos, h)
    % Add a reference point from which to start a correlation, to inform
    % kmeans clustering how to cluster a ROI
    data = getappdata(h.hGUI, 'data');
    pos = round(pos);
    
    pos1D = (pos(2)-1)*data.dim(2)+pos(1);
    if ismember(pos1D, data.pos)
        data.corners = [pos', data.corners];
        setappdata(h.hGUI, 'data', data)
        ClusterRoi(h)
    end
end


function DeleteCluster(pos, h)
    % Remove a cluster so it gets rejected immediately after applying clustering
    
    data = getappdata(h.hGUI, 'data');
    pos = round(pos);
    
    pos1D = (pos(2)-1)*data.dim(2)+pos(1);
    if ismember(pos1D, data.pos)
        toDel = data.clIdx(data.pos==pos1D);
        if toDel ~= 0
            data.clIdx(data.clIdx==toDel) = 0; % Set the clicked cluster to 0

            % Update the mask
            img = permute(h.clusterPlot2.Children(end).CData, [2 1 3]);
            for i = 1:3
                imgi = img(:,:,i);
                imgi(data.pos(data.clIdx==0)) = 0.1; % Change the color to dark
                img(:,:,i) = imgi;
            end
            % Change the plots.
            h.signalAx2.Children(end-toDel+1).Color = ...
                [h.signalAx2.Children(end-toDel+1).Color, 0.1];
            h.clusterPlot2.Children(end-toDel).Color = [0.5 0.5 0.5];
            h.clusterPlot2.Children(end).CData = permute(img, [2 1 3]);
            setappdata(h.hGUI, 'data', data)
        end
    end
end


function nClusterSlider_Callback(hObject, ~, h) %#ok
    % Changes the number of cluster kmeans finds in ROIs
    hObject.Value = round(hObject.Value);
    switches = getappdata(h.hGUI, 'switches');
    switches.nCluster = hObject.Value;
    h.nClusterTitle.String = sprintf('number of clusters: %d', hObject.Value);
    setappdata(h.hGUI, 'switches', switches)
    ClusterRoi(h)
end


% --- Executeson button press in applyClustering.
function applyClustering_Callback(~, ~, h) %#ok
    idx = getappdata(h.hGUI, 'idx');
    data = getappdata(h.hGUI, 'data');
    if size(data.PP.P,1)==5 % if pixel correlations are in the PP.P variable
        sbxt = getappdata(h.hGUI, 'sbxt'); % load the raw data
    end
    PP = data.PP;
    Mask = data.Mask;
    C = data.clCon;
    switches = getappdata(h.hGUI, 'switches');
    roi = switches.splitRoi;
    
    if switches.splitStarted
        switches.splitStarted = false;
        h.applyClustering.String = 'applying...';
        h.applyCluster = TurnOn(h.applyClustering);
        n = unique(data.clIdx); % which clusters to be saved
        n(n==0) = [];
        nnew = length(n);

        % The old roi gets deleted
        threshold = PP.P(3,roi);
        PP.Con(roi) = [];
        PP.A(roi) = [];
        PP.P(:,roi) = [];
        PP.SpecProfile(:,roi) = [];
        PP.peakFreq(roi) = [];
        PP.creationMethod(roi) = [];
        Mask(Mask==roi) = 0;
        Mask(Mask>roi) = Mask(Mask>roi) - 1; % Close the gap that was just created
        Mask(Mask>0) = Mask(Mask>0) + nnew; % Make room for the new ROIs
        PP.Cnt = PP.Cnt - 1 + nnew;

        % New ones get inserted
        for i = nnew:-1:1
            
            Mask = Mask';
            Mask(data.pos(data.clIdx==n(i))) = i;
            Mask = Mask';
            
            newCon = struct('x', C(n(i)).y, 'y', C(n(i)).x);
            PP.Con = [newCon, PP.Con];
            
            maxval = max(data.BImg(Mask==i));
            xmass = round(mean(C(n(i)).y)); % centre of mass of this new ROI
            ymass = round(mean(C(n(i)).x));
            if size(PP.P, 1)==5
                % Calculate the average correlation of the signal of each
                % pixel in the ROI with the signal of the 3 pixels around
                % the brightest pixels of the ROI per cluster
                cord1D = sub2ind(size(data.Mask'),xmass, ymass);
                signalCenter = mean(sbxt.Data.y(:,cord1D-1:cord1D+1),2);
                signalCenter = decimate(double(signalCenter), floor(data.freq));
                signalOthers = sbxt.Data.y(:,data.pos(data.clIdx==n(i)));
                signalOthersdecim = zeros(ceil(size(signalOthers,1)/floor(data.freq)), size(signalOthers,2));
                for j = 1:size(signalOthers,2)
                    signalOthersdecim(:,j) = decimate(double(signalOthers(:,j)),floor(data.freq));
                end
                corVal = mean(corr(signalCenter, signalOthersdecim, 'rows', 'pairwise'));
                PP.P = [[xmass; ymass; maxval; threshold; corVal], PP.P];
            else
                PP.P = [[xmass; ymass; maxval; threshold], PP.P];
            end
            PP.A = [sum(Mask(:)==i), PP.A];
            PP.creationMethod = [{'splitted'}; PP.creationMethod];
            
            % Calculate spectral profiles
            nSpec = size(data.imgStackT,3);
            ROIi = Mask == i;
            ROIi3D = repmat(ROIi, [1,1,nSpec]);
            specProfile = reshape(data.imgStackT(ROIi3D), [sum(ROIi(:)), nSpec]);
            PP.SpecProfile = [mean(log(specProfile))' PP.SpecProfile];
            [~, peakFreq] = max(PP.SpecProfile);
            PP.peakFreq = [data.Sax(peakFreq), PP.peakFreq];
        end

        % Update the rejection lists
        idx.Size(roi)  = [];
        idx.Thres(roi) = [];
        idx.ThresCor(roi) = [];
        idx.Round(roi) = [];
        idx.Black(roi) = [];
        idx.White(roi) = [];
        idx.Size  = [false(1,nnew),idx.Size];
        idx.Thres = [false(1,nnew),idx.Thres];
        idx.ThresCor = [false(1,nnew),idx.ThresCor];
        idx.Round = [false(1,nnew),idx.Round];
        idx.Black = [false(1,nnew),idx.Black];
        idx.White = [true(1,nnew), idx.White]; % immediately whitelist

        % Save the clustered ROIs
        data.PP = PP;
        data.Mask = Mask;
        setappdata(h.hGUI, 'data', data)
        setappdata(h.hGUI, 'switches', switches)
        setappdata(h.hGUI, 'idx', idx)
        
        % Redraw the roi contours
        h.mainAx.NextPlot = 'replacechildren';
        h.im = imagesc(h.im.CData, 'Parent',h.mainAx, 'HitTest','off');
        DrawRois(h)
        h = guidata(h.hGUI); % Update the handles, were edited by DrawRois
        UpdateRois(h)
        h.applyClustering.String = 'clustering applied!';
    else
        % Clustering not started/ just finished-> Don't change ROIs!
        h.applyClustering.String = 'not clustering atm!';
    end
    
    try 
        pause(2)
        h.applyClustering.String = 'Apply Clustering';
        h.applyClustering = TurnOff(h.applyClustering);
    catch
        fprintf('you closed the figure quickly after applying clustering\n')
    end
end


%% Data selection functions %%

% --- Executes on slider movement.
function numSelectSlider_Callback(hObject, ~, h)
    % Change the number of signal traces that will be visible
    switches = getappdata(h.hGUI, 'switches');
    colors   = getappdata(h.hGUI, 'colors');
    
    hObject.Value = round(hObject.Value); % Number of max data tracers
    
    if switches.selNumCur > hObject.Value % Current datatrace to draw
        switches.selNumCur = 1;
    end
    
    h.signalAx1.NextPlot = 'add';
    h.signalAx1.YDir = 'normal';
    h.mainAx.NextPlot = 'add';
    if isfield(h,'signalPlt')
        fprintf('deleting existing lines & handles\n')
        delete(h.signalPlt) % delete the lines
        delete(h.selectedPlt) % delete the boxes
    end

    % Create handles for the data lines
    h.signalPlt = plot(zeros(2, hObject.Value), zeros(2, hObject.Value),...
        'Parent', h.signalAx1); % replace or create handles and lines
        
    % Create handles for the selection boxes
    h.selectedPlt = plot(zeros(2, hObject.Value), zeros(2, hObject.Value),...
        'Parent', h.mainAx);
        
    switches.selNumMax = hObject.Value;
    
    colors.summer = summer(switches.selNumMax+1);
    % Change green value to make it look better
    colors.summer(:,2) = colors.summer(:,2).*1.25-0.25; 
    
    h.numSelectTitle.String = sprintf('number of selections: %d', hObject.Value);
    guidata(h.hGUI, h)
    setappdata(h.hGUI, 'switches', switches)
    setappdata(h.hGUI, 'colors', colors)
end 


% --- Executes on slider movement.
function selSize_Callback(hObject, ~, h) %#ok
    % Change the amount of pixels a dataselection will be wide/2
    switches = getappdata(h.hGUI, 'switches');
    switches.selSize = round(hObject.Value);
    hObject.Value = switches.selSize;
    setappdata(h.hGUI, 'switches', switches)
    
    % Update slider title text
    val = hObject.Value*2;
    h.selSizeTitle.String = {sprintf('selection size : %2.0fx%2.0f=%4.0f pixels', val, val, val*val)};
end


function [xmesh, ymesh, signal, varargout] = SelectData(pos, selSize, h, modus)
    % Selects data from the sbx file using coordinates the main axes image
    sbxt = getappdata(h.hGUI, 'sbxt');
    data = getappdata(h.hGUI, 'data');
    switches = getappdata(h.hGUI, 'switches');
    dim = data.dim;
    x = round(pos(1));
    y = round(pos(2));
    
    % Make sure the selection window does not fall outside of the BImg
    if x-selSize < 1
        x = selSize+1;
    elseif x+selSize > dim(2)
        x = dim(2)-selSize-1;
    end
    if y-selSize < 1
        y = selSize+1;
    elseif y+selSize > dim(3)
        y = dim(3)-selSize-1;
    end
    
    % Get 2D indexes of data
    xmesh = meshgrid(x-selSize:x+selSize);
    ymesh = meshgrid(y-selSize:y+selSize)';

    % Convert 2D indexes to 1D indexes
    ind = xmesh + (ymesh-1).*dim(2);
    ind = reshape(ind, [numel(ind),1]);

    % Select the data
    if strcmp(modus, 'mean')
        signal = mean(sbxt.Data.y(:,ind),2);
    else
        signal = sbxt.Data.y(:,ind);
        
        % Select nine pixels as center of data (to reduce noise)
        x2 = meshgrid(x-1:x+1);
        y2 = meshgrid(y-1:y+1)';
        ind = x2 + (y2-1).*dim(2);
        ind = reshape(ind, [numel(ind),1]);
        varargout{1} = mean(sbxt.Data.y(:, ind),2);
    end
end


%% Data visualisation functions %%
% --- Executes on slider movement.
function AlphaSlider_Callback(hObject, ~, h) %#ok
    % Determines the alpha value for the ROI contour lines
    switches = getappdata(h.hGUI, 'switches');
    switches.alph = hObject.Value; % Update the alpha value
    h.alphaSliderTitle.String = sprintf('Contour Alpha: %.2f',hObject.Value);
    setappdata(h.hGUI, 'switches', switches)
    UpdateRois(h) % Immediately update the contours
end


% --- Executes on selection change in backGrdView.
function backGrdView_Callback(hObject, ~, handles) %#ok
    backGrdView(hObject.Value, handles)
end


function BackgrdCLimStart(h)
    % Starts a nice histogram and the caxis lines at the selection of new
    % background view
    switches = getappdata(h.hGUI, 'switches');
    
    lineW = 2;
    lineC = 'k';
    lineVis = '-';
    
    h.backgrdCLimAx.NextPlot = 'replacechildren';
    img = h.im.CData(:); % linearized background image
    if isfield(h, 'backgrdHist')
        delete(h.backgrdHist)
    end
    h.backgrdHist = histogram(h.backgrdCLimAx, img, 100,...
        'EdgeAlpha', 0, 'HitTest', 'off');
    
    if isfield(h, 'backgrdCLimLine1')
        delete(h.backgrdCLimLine1)
        delete(h.backgrdCLimLine2)
    end
    backgrdCLim = prctile(img, [0.01 99.99]);
    
    ylims = get(h.backgrdCLimAx, 'ylim');
    
    h.backgrdCLimAx.NextPlot = 'add';
    h.backgrdCLimLine1 = plot(backgrdCLim([1 1]), ylims,...
        'linewidth',lineW, 'color', lineC, 'linestyle', lineVis,...
        'Parent', h.backgrdCLimAx, 'HitTest','off');
    h.backgrdCLimLine2 = plot(backgrdCLim([2 2]), ylims,...
        'linewidth',lineW, 'color', lineC, 'linestyle', lineVis,...
        'Parent', h.backgrdCLimAx, 'HitTest','off');
    h.backgrdCLimAx.YDir = 'normal';
    h.backgrdCLimAx.YTick = []; % Turn off y axis ticks
    h.backgrdCLimAx.XTick = []; % Turn off x ticks
    
    switches.backgrdCLim = backgrdCLim;
    
    % Save the edited data
    setappdata(h.hGUI, 'switches', switches)
    guidata(h.hGUI, h)
    
    backgrdCLim_Apply(backgrdCLim(2),h) % Update the color axis correctly
end

function backGrdView(selected, h)
    % Changes the background image of mainAx
    
%     currentView = [h.mainAx.XLim, h.mainAx.YLim];
    switches = getappdata(h.hGUI, 'switches');
    % If the background image should be changed
    if selected ~= switches.viewToggle
        switches.viewToggle = selected;
        data = getappdata(h.hGUI, 'data');
        h.mainAx.NextPlot = 'add';
    
    
        % Legend text off
        h.legendText2.Visible = 'off';
        h.legendText3.Visible = 'off';
        h.legendText4.Visible = 'off';
        h.legendText5.Visible = 'off';
        h.legendText6.Visible = 'off';
        
        switch selected
            case 1 % Spectral image
                h.im.CData = data.BImg;
                colormap(gray)
                
            case 2 % Spectral image, colored by frequency
                
                % Cheat view toggle so you don't have to switch to other
                % view before you can select view 2 again
                switches.viewToggle = 999; 
                
                s = data.Sax; % spectral frequency axis
                % Get which frequencies to show
                question = [sprintf('which frequency (range) to select? min=%.2fHz, max=%.2fHz\n',...
                                    s(1), s(end)),...
                            'If you enter nothing all frequencies will be selected',...
                            sprintf('\nlowest selected frequency will be red, highest will be blue\n'),...
                            'colormap will be jet'];
                answer = inputdlg(question, 'spectral colored by intensity at different frequencies');
                if isempty(answer)  % cancel was pressed
                    sselect = 1:length(s); % spectral frequency selection
                elseif isempty(answer{:}) % ok was pressed without input
                    sselect = 1:length(s);
                else
                    stoselect = sort(str2num(answer{:})); % which frequency range to select
                    if isempty(stoselect) % letters were given so now stoselect is empty
                        stoselect = s([1 end]); % select all frequencies
                    elseif length(stoselect) == 1
                        % If only one frequency was given, make sure that a
                        % frequency can be used by using a range around the given frequency
                        fstep = mean(diff(s)); % difference in frequencies
                        stoselect = [stoselect(1)-fstep*2, stoselect(1)+fstep*2];
                    end
                    sselect = find(s>=stoselect(1) & s<=stoselect(2));
                end
                
                nsselect = length(sselect);
                if nsselect == 0 % no frequency was selected because of bad input range, select lowest frequency
                    sselect  = 1:length(s);
                    nsselect = length(sselect);
                end
                colors = flipud(jet(nsselect));
                
                h.im.CData = CreateRGB2(num2cell(log1p(data.imgStackT), [1 2]), sselect, colors, true, true);
                
            case 3 % SpatialCorr
                h.im.CData = data.SpatialCorr;
                colormap(jet)
                
            case 4 % ROI corners signal correlations image
                if switches.roiCorrIm == false % image needs to be created
                    sbxt = getappdata(h.hGUI, 'sbxt');
                    
                    try % If ROIs are too small they get deleted when trying to remove
                        % the outer edge of the ROI
                        [piece, check] = GetRoiCorners(data.Mask, data.PP, 2, 20, 'size');
                    catch
                        fprintf('\nprobably some ROIs are too small, adjusting buffer zone\n')
                        [piece, check] = GetRoiCorners(data.Mask, data.PP, 0, 20, 'size');
                    end

                    [newMask, ~, ~, ~, ~] = CorrRoiCorners(piece, check, data.Mask, sbxt, data.dim);

                    % Plotting
                    newMaskc = squeeze(num2cell(newMask, [1 2]));
                    [data.roiCorrImg, ~] = CreateRGB(newMaskc, 'r g b gb');                 
                    data.roiCorrImg = permute(data.roiCorrImg, [2 1 3]);
                    h.im.CData = data.roiCorrImg;
                    
                    switches.roiCorrIm = true;
                    setappdata(h.hGUI, 'data', data)
                 
                else % Image is already made
                    h.im.CData = data.roiCorrImg;
                end
                
                colormap(gray)
                
                % info text about colors
                h.legendText2.String = 'left';
                h.legendText2.ForegroundColor = 'r'; % red
                h.legendText2.Visible = 'on';
                h.legendText3.String = 'right';
                h.legendText3.ForegroundColor = 'g'; % green
                h.legendText3.Visible = 'on';
                h.legendText4.String = 'top';
                h.legendText4.ForegroundColor = 'b'; % blue
                h.legendText4.Visible = 'on';
                h.legendText5.String = 'bottom';
                h.legendText5.ForegroundColor = 'c'; % cyan = green + blue
                h.legendText5.Visible = 'on';
                h.legendText6.String = 'corner correlations good overlap';
                h.legendText6.ForegroundColor = 'w'; % white
                h.legendText6.Visible = 'on';
                
            case 5 % Mask view
                picMask = data.Mask;
                picMask(picMask==0) = data.PP.Cnt+100;
                h.im.CData = picMask;
                colormap(flipud(gray))
                
            case 6 % Chronic & current backgroudn overlay view
                if ~switches.chronic % Load the chronic recording
                    importChronicFile(h)
                    data = getappdata(h.hGUI, 'data');
                    switches = getappdata(h.hGUI, 'switches');
                end
                h.im.CData = CreateRGB({data.BImg, data.chronicImg}, 'g r');
                
            case 7 % Chronic view only
                if ~switches.chronic % Load the chronic recording
                    importChronicFile(h)
                    data = getappdata(h.hGUI, 'data');
                    switches = getappdata(h.hGUI, 'switches');
                end
                h.im.CData = CreateRGB({data.chronicImg}, 'r');
                
            case 8 % Colors are based on what the peak frequency is
                picMask  = data.Mask;
                img = zeros(size(picMask));
                for i = 1:data.PP.Cnt
                    img(picMask==i) = data.PP.peakFreq(i);
                end
                h.im.CData = img;
                colors = [0 0 0; jet(255)]; % jet with black
                colormap(colors)
                
            case 9 % Import variable from workspace
                
                % cheat view toggle so you don't have to switch to other view before you can import another time
                switches.viewToggle = 999; 
                
                infoStr = ['variable should result in a matrix of doubles.',...
                           'first two dimensions should same as background image, may be 3D. ',...
                           'three examples: BImgA   or   permute(log1p(SPic), [2 1 3])   or   cat(3, BImgAverage, BImgMax)'];
                varName = inputdlg(infoStr, 'Import a variable from the workspace');
                
                % quit if nothing was typed
                if isempty(varName) 
                    return
                end
                
                try
                    impImg = evalin('base',varName{:}); % load the requested variable
                    originDim = size(data.BImg);
                    impDim = size(impImg);
                    if isa(impImg, 'double')...
                            && impDim(1) == originDim(1)...
                            && impDim(2) == originDim(2)
                        
                        % Create a color image if necessary
                        if length(impDim) == 3
                            if impDim(3)==2
                                % 2 colors : red, green
                                colors = [1 0 0; 0 1 0]; 
                            elseif impDim(3)==3
                                % 3 colors : red, green, brighter blue
                                colors = [1 0 0; 0 1 0; 0.1 0.1 1];
                            else
                                % more colors: jet colormap
                                colors = flip(jet(impDim(3))); 
                            end
                            impImgCell = squeeze(num2cell(impImg, [1 2]));
                            answer = questdlg('normalize color image for every data slice and color?',...
                                'normalize?',...
                                'yes','no','no');
                            switch answer
                                case 'yes'
                                    norma = true;
                                case 'no'
                                    norma = false;
                            end
                            impImg = CreateRGB2(impImgCell, 1:impDim(3), colors, norma, norma);
                        elseif length(impDim) > 3
                            msgbox('variables with more than 3 dimensions not supported', 'error')
                        end
                        
                        % Apply the chosen variable
                        h.im.CData = impImg;
                    end
                catch ME
                    str = sprintf('unable to load requested variable!: %s', varName{:});
                    msgbox(str, 'error')
                    throw(ME)
                end
        end
        h.im.HitTest = 'off';
        caxis([min(h.im.CData(:)), max(h.im.CData(:))])
        h.mainAx.YDir = 'reverse';
        setappdata(h.hGUI, 'switches', switches); % Save current backgroundview
        guidata(h.hGUI, h); % Save the handle of im

        % Start up the clim selecter axis
        BackgrdCLimStart(h)
    end
end

function importChronicFile(h)
    % Load a chronic file, and register the averaged background images to
    % the current background image
    data = getappdata(h.hGUI, 'data');
    switches = getappdata(h.hGUI, 'switches');
    [f,p] = uigetfile('','load the chronic file for this dataset');

    chronicImg = load([p,f],'BImg2');
    chronicImg = chronicImg.BImg2;
    dims = size(chronicImg{1});
    % Average the other background images to one nice image
    chronicImg = mean(reshape(cell2mat(chronicImg),[dims(1),dims(2),length(chronicImg)]),3)';

    % Register the chronic dataset to this dataset
    waBa = waitbar(1/3, 'translating chronic image'); % Waitbar
    buffer = 100;
    reference = chronicImg(buffer:size(data.BImg,1)-buffer, buffer:size(data.BImg,2)-buffer);
    reference = reference-mean(reference(:));
    base = data.BImg-mean(data.BImg(:));
    correl = xcorr2(base,reference);
    [~,snd] = max(correl(:));
    [ij, ji] = ind2sub(size(correl),snd);
    x = size(base,1) - ij - buffer;
    y = size(base,2) - ji - buffer;
    % Cutting or padding the chornic image if necessary
    if x>=1 % Cut top
        chronicImg = chronicImg(x:end, :);
    else % Padd top % Maybe add a +1; code untested because condition doesn't happen often
        chronicImg(-x+1:-x+size(chronicImg,1),:) = chronicImg;
    end
    if y>=1 % Cut left beginning
        chronicImg = chronicImg(:, y:end);
    else % Padd left beginning
        chronicImg(:,-y+1:-y+size(chronicImg,2)) = chronicImg;
    end
    if size(chronicImg,1) > size(data.BImg,1) % Cut bottom
        chronicImg(size(data.BImg,1)+1:end,:) = [];
    else % padd bottom
        chronicImg(size(data.BImg,1),1) = 0;
    end
    if size(chronicImg,2) > size(data.BImg,2) % Cut right side
        chronicImg(:,size(data.BImg,2)+1:end) = [];
    else % padd right side
        chronicImg(1,size(data.BImg,2)) = 0;
    end

    % Rotation by checking rotations
    waitbar(2/3, waBa, 'rotating chronic image'); % Waitbar
    rotation = -1:0.1:1; % rotations to apply: counterclockwise to clockwise degrees
    simil  = zeros(1, length(rotation));
    % correl = zeros(nfiles, length(rotation));
    buffer = 30; % buffer image edges because rotation makes part fall off
    BImgRef = data.BImg(buffer:end-buffer, buffer:end-buffer);
    for j = 1:length(rotation)
        BImgRot = imrotate(chronicImg, rotation(j),'bicubic','crop');
        BImgRot = BImgRot(buffer:end-buffer, buffer:end-buffer);
        simil(j)= corr2(BImgRot, BImgRef);
    end

    % Best rotations
    [~, rotIdx] = max(simil);
    rotAng = rotation(rotIdx);

    chronicImg = imrotate(chronicImg, rotAng,'bicubic', 'crop');

    % Save the chronic, registered background image
    data.chronicImg = chronicImg;
    switches.chronic = true;
    close(waBa)

    setappdata(h.hGUI, 'data', data)
    setappdata(h.hGUI, 'switches', switches)
end

function UpdateMainImg(h)
    % Update plot
    h.mainAx.NextPlot = 'replacechildren';
    h.im = imagesc(h.im.CData, 'Parent',h.mainAx, 'HitTest','off');
    % Redraw the roi contours
    DrawRois(h)
    h = guidata(h.hGUI); % Update the handles, was edited by DrawRois
    UpdateRois(h)
    RejectInfoFunc(h)
end


function DrawRois(h)
    % Draws the ROI contours and centres.
    
    data = getappdata(h.hGUI, 'data');
    PP = data.PP;
    
    h.Con = gobjects(size(PP.Con,2),1);
    h.PPP = gobjects(size(PP.P,2),1);
    
    h.mainAx.NextPlot = 'add';
    for i = 1:PP.Cnt
        h.Con(i) = plot(PP.Con(i).x, PP.Con(i).y,...
            'Parent',h.mainAx,'hitTest','off');
%         hCon(i).ZData = ones(1,length(PP.Con(i).x))./10;
        h.PPP(i) = plot(PP.P(1,i), PP.P(2,i), 'x', 'markersize',4,...
            'Parent',h.mainAx,'hitTest','off');
    end
    
    h.manRoi = gobjects(1);
    
    guidata(h.hGUI,h); % Save the handles
end


function UpdateRois(h)
    % Update colors and alpha of the ROI contours
    
    idx = getappdata(h.hGUI, 'idx');
    
    % get alpha value
    switches = getappdata(h.hGUI, 'switches');
    alph = switches.alph;
    
    % colors
    cyan    = [0, 1, 1];
    red     = [1, 0, 0];
    green   = [0, 1, 0];
    magenta = [1, 0, 1];
    for j = 1:length(h.Con)
        set(h.Con(j), 'Color',[cyan, alph], 'linewidth',0.5);
        set(h.PPP(j), 'Color', cyan);
        if h.thresToggle.Value
            if idx.Thres(j) == 1
                set(h.Con(j), 'Color', [red, alph]);
            end
        else % The correlation threshold should be used
            if idx.ThresCor(j) == 1
                set(h.Con(j), 'Color', [red, alph]);
            end
        end
        if idx.Round(j) == 1
            set(h.PPP(j), 'color', magenta)
        end
        if idx.Size(j) == 1
            set(h.PPP(j), 'Color', red);
        end
        if idx.White(j) == 1
            set(h.Con(j), 'Color', [green, (alph/2)+0.5],'linewidth',1.2)
        end
        if idx.Black(j) == 1
            set(h.Con(j), 'Color', [magenta, (alph/2)+0.5], 'linewidth',1.2)
        end
    end
end


function DrawData(xm, ym, signal, h)
    % Draw raw 2P data after it is retreived by SelectData
    
    switches = getappdata(h.hGUI, 'switches');
    data   = getappdata(h.hGUI, 'data');
    colors = getappdata(h.hGUI, 'colors');
    
    num = switches.selNumCur; % current trace to draw
    
    % Delete the previous line
    delete(h.signalPlt(num))
    delete(h.selectedPlt(num))
    
    % plot the signal
    h.signalPlt(num) = plot(data.xas, signal,'color',colors.summer(num,:),...
        'Parent',h.signalAx1);
    
    % Save extreme values for calculation of nice YLims
    switches.ylims(1,num) = min(signal);
    switches.ylims(2,num) = max(signal);
    
    % Draw the box area from which the data is sampled
    x = [xm(1), xm(1), xm(end), xm(end), xm(1)]-0.5;
    y = [ym(1),ym(end),ym(end),ym(1), ym(1)]-0.5;
    h.selectedPlt(num) = plot(x,y,'color', colors.summer(num,:),...
        'linewidth',2,'hittest','off','Parent',h.mainAx);
    
    % And set the nice limits for the data
    h.signalAx1.XLim = [data.xas(1) data.xas(end)];
    h.signalAx1.YLim =  [round(min(switches.ylims(1,:)),-2)-100,...
                        round(max(switches.ylims(2,:)),-2)+100];
    
    switches.selNumCur = mod(num, switches.selNumMax) + 1;
    setappdata(h.hGUI, 'switches', switches)
    guidata(h.hGUI, h)
end


function LocalCorr(xm, ym, signal, centersignal, h)
    % Calculate correlations between the centersignal and all rows of the
    % signal, transform those correlations into the 2D image and overlay it
    % on the background image
    
    colors = getappdata(h.hGUI, 'colors');
    
    correl = corr(double(centersignal), double(signal),'rows','pairwise');
    correl = round(correl*128)+128;
    correl = colors.parula(correl,:);
    correl = reshape(correl, [size(xm),3]);
    
    if isfield(h, 'correl')
        delete(h.correl)
    end
    h.correl = imagesc(xm(1,:), ym(1,:), correl,...
        'Parent', h.mainAx, 'hittest', 'off');
    
    % Save the new h.correl handle
    guidata(h.hGUI, h)
end


%% Axes click callback
% --- Executes on mouse press over axes background.
function mainAx_ButtonDownFcn(~, eventdata, h) %#ok<DEFNU>
    % Gives the position of a click on the main axes to a function if a
    % function that can use clicks of the main axes is turned on
    position = eventdata.IntersectionPoint([1 2]);
    switches = getappdata(h.hGUI, 'switches');
    
    if switches.currentMajor == 2 % The ROI split panel
        % ROI splitting panel is active -> be ready for ROI selection
        RoiSplitter(position, h)
        
    elseif switches.currentMajor == 3 % The ROI creation panel
        % The ROI creation panel is active -> be ready for ROI selection
        switches.creationPos = round(position);
        setappdata(h.hGUI, 'switches', switches)
        CreateRoi(h,'new')
        
    elseif switches.currentMajor == 5 % The Manual ROI creation panel
        % The manual ROI creation panel is active -> be ready for Contour
        % selection
        ManualRoi(position, eventdata.Button, h)
        
    elseif strcmp(switches.currentMinor, 'whiteButton')
        % select ROI for white listing
        SelectROI(position, h, 'white')
        
    elseif strcmp(switches.currentMinor, 'blackButton')
        % select ROI for black listing
        SelectROI(position, h, 'black')
        
    elseif strcmp(switches.currentMinor, 'sigSelectButton')
        % Select 2P signal from mouse click surrounding selsize x selsize pixels
        [xmesh, ymesh, signal] = ...
            SelectData(position, switches.selSize, h, 'mean');
        % Plot the select data
        DrawData(xmesh, ymesh, signal, h)
        h = guidata(h.hGUI);
        guidata(h.hGUI, h)
        
    elseif strcmp(switches.currentMinor, 'corrButton')
        % Local correlation visualisation
        
        [xmesh, ymesh, signal, centersignal] = SelectData(position, 16, h, '');
        
        % subsampling for long signals
        time = 1:length(centersignal);
        if length(centersignal) > 1000
            step = round(length(centersignal)./1000);
            time = 1:step:length(centersignal);
        end
        signal = signal(time,:);
        centersignal = centersignal(time);
        
        LocalCorr(xmesh, ymesh, signal, centersignal, h)
        
    elseif strcmp(switches.currentMinor, 'plotROIsButton')
        % Plot the ROI signal that is clicked
        SelectROI(position, h, 'nothing')
    end
end


% --- Executes on mouse press over axes background.
function clusterPlot_ButtonDownFcn(~, eventdata, h) %#ok
    % Both cluster plots, which show data for ROI splitting use this
    % callback
    
    position = eventdata.IntersectionPoint([1 2]);
    switches = getappdata(h.hGUI, 'switches');
    
    if strcmp(switches.currentMinor, 'addPointButton')
        if switches.splitStarted
            AddPoint(position, h)
        end
    elseif strcmp(switches.currentMinor, 'delClusterButton')
        if switches.splitStarted
            DeleteCluster(position, h)
        end
    end
end


% --- Executes on mouse press over axes background.
function backgrdCLimAx_ButtonDownFcn(~, eventData, h) %#ok<DEFNU>
    % The plot in Show Data tab uses this callback to set the color axis
    % limits of the main background image
    x = eventData.IntersectionPoint(1);
    
    backgrdCLim_Apply(x, h)
end

function backgrdCLim_Apply(x, h)
    % Change the color axis or color values of the main image
    
    switches = getappdata(h.hGUI, 'switches');
    
    xMin = switches.backgrdCLim(1);
    xMax = switches.backgrdCLim(2);
    
    % Calculate to which line the click was closest
    difference = abs([xMin-x, xMax-x]);
    [~, idx] = min(difference);
    
    % Edit the line that was closest to the click
    if idx == 1
        switches.backgrdCLim(1) = x;
        h.backgrdCLimLine1.XData = [x x];
    else
        switches.backgrdCLim(2) = x;
        h.backgrdCLimLine2.XData = [x x];
    end
    
    % Set the color limits
    if length(size(h.im.CData))==2 % if the background image is 2D use Caxis
        if diff(switches.backgrdCLim) > 0 % if CLimits are valied apply them
            set(h.mainAx, 'clim',switches.backgrdCLim)
        else % In case CLimits are not increasing, force them
            set(h.mainAx, 'clim', [switches.backgrdCLim(1), switches.backgrdCLim(1)+0.1])
        end
    else % 3D CData requires direct mapping of the color values
        % reset colors to original, which should always be range 0-1
        img = (h.im.CData - min(h.im.CData(:))) / range(h.im.CData(:));
        % set colors to requested
        img = img / switches.backgrdCLim(2) - switches.backgrdCLim(1);
        h.im.CData = img;
    end
    
    setappdata(h.hGUI, 'switches', switches)
end

%% Save changes function
% --- Executes on button press in saveButton.
function saveButton_Callback(~, ~, h) %#ok
    % Save the new ROIs to workspace
    data = getappdata(h.hGUI, 'data');
    switches = getappdata(h.hGUI, 'switches');
    PP = data.PP;
    Mask = data.Mask;
    BImg = data.BImg;
    SpatialCorr = data.SpatialCorr;
    SpatialCorr(Mask==0) = 0;
    
    if isequal(data.PP.A, switches.originalROIs)
        h.saveButton.String = {'No change'};
        pause(2)
        try
            h.saveButton.String = {'Save'};
        catch
            fprintf('GUI is closed probably.\n')
        end
    else
        h.saveButton = TurnOn(h.saveButton);
        
        % If at least one of the sliders is edited, save their settings
        saveSliders = ~all(isnan(table2array(switches.sliderSettings)));
        
        if isfield(data, 'SPSIGfile')
            h.saveButton.String = {'to SPSIG file'};
            pause(0.01)
            save(data.SPSIGfile, 'PP', 'Mask', 'BImg', 'SpatialCorr', '-append')
            if saveSliders
                rrSliderSettings = switches.sliderSettings;
                save(data.SPSIGfile, 'rrSliderSettings', '-append')
            end
            % TODO: update the SpatialCorr!!!
            disp('Changes saved to SPSIG file & workspace')
        else
            h.saveButton.String = {'to workspace'};
        end
        
        assignin('base', 'Mask', Mask);
        assignin('base', 'PP', PP);
        assignin('base', 'BImg', BImg);
        assignin('base', 'SpatialCorr', SpatialCorr);
        if saveSliders
            rrSliderSettings = switches.sliderSettings;
            assignin('base', 'rrSliderSettings', rrSliderSettings);
        end
        
        switches.originalROIs = data.PP.A;
        setappdata(h.hGUI,'switches',switches);
        
        try
            pause(1)
            h.saveButton.String = {'Saved'};
            pause(2)
            h.saveButton.String = {'Save'};
            h.saveButton = TurnOff(h.saveButton);
        catch
            fprintf('You closed fast after saving. Saving successfull\n')
        end
    end
end
