function Showsbx(varargin)

global info Pos clim Chan %bRun SliderChan
persistent Crop

if nargin == 1
    strfp = varargin{1};
    [Pn,  Fn, ~] = fileparts(strfp);
    filename{1} = Fn;
elseif nargin == 3   %callback
    Pn = varargin{3};
    [Fn , Pn] = uigetfile([Pn '*.sbx']);
    if ~ischar(Fn)
        return;
    end
    filename = strsplit(Fn, '.');
    strfp = [Pn filename{1}];
else
    [Fn , Pn] = uigetfile('*.sbx');
    if ~ischar(Fn)
        return;
    end
    filename = strsplit(Fn, '.');
    strfp = [Pn filename{1}];
end

hmb = msgbox('Opening sbx file........');
Img = sbxread(strfp, 0,1);
close(hmb)

Pos = 1; %frame position
bRun = 0; %play sbx movie
Chan = 1; %which channel
Slice = 0; %which level
bsplit = 0; %don't show slices
bVers = 0; %gpu or master
Splitnm = 1;
bUpCrop = 0;%has the crop been changed after loading the stack
bFilt = 1; %Median filter images (Yes(1)/No(0))

%number of splits
if ~isfield(info, 'Slices')
    if isfield(info, 'otparam') && ~isempty(info.otparam)
        Splitnm = info.otparam(3);
    else
        answ = inputdlg('How many Slices?', 'Slice info!!...', [1 40],  {'1'});
        Splitnm = str2double(answ{1});
    end
    
    info.Slices = Splitnm;
    info.Slice = Slice;
    info.bsplit = bsplit;
else
    Splitnm = info.Slices;
    if isfield(info, 'Slice')
        Slice = info.Slice;
    else
        info.Slice = 0;
    end
    if isfield(info, 'bsplit') && Splitnm > 1
        bsplit = info.bsplit;
    else
        info.bsplit = 0;
    end
end

framenm = info.max_idx; %total number of frames

if isfield(info, 'crop')
    %if imported sbx already has been cropped
    Crop = info.crop;
    % Set crop to position that doesn't cause crashes
    if isfield(info, 'simon') || (isfield(info, 'scanbox_version') &&  info.scanbox_version == 2.5 )
        dwidth = info.Shape(info.Perm(2));
        dheight = info.Shape(info.Perm(1));
    else
        dwidth = info.sz(2);
        dheight = info.sz(1);
    end
    if Crop.x(end)>dwidth
        Crop.x = Crop.x - (Crop.x(end)-dwidth);
    end
    if Crop.y(end)>dheight
        Crop.y = Crop.y - (Crop.y(end)-dheight);
    end
    if length(Crop.x)>dwidth
        fprintf('Adjusting crop.x settings to fit data\n')
        Crop.x = 1:dwidth;
    end
    if length(Crop.y)>dheight
        fprintf('Adjusting crop.y settings to fit data\n')
        Crop.y = 1:dheight;
    end
end
if ~isempty(Crop)
    %by default use previous crop parameters
    info.crop = Crop;
    info.d1 = length(Crop.x);
    info.d2 = length(Crop.y);
end

if size(Img,1) < 3
    bVers = 1;
    Img = permute(Img, [2 3 1]);
end

info.bVers = bVers;
lnChan = size(Img,3); %number of channels
info.strfp = strfp;

%determine pixel value range
clim = prctile(unique(Img(:)), [2 99.8]);
%clim = prctile(Img(:), [5 95]);

f  = figure('ToolBar','none');

scrnSz = get(0,'ScreenSize');
desiredSz = [min(725, scrnSz(3)), min(500, scrnSz(4))];
p1 = round(scrnSz(3)/2)-round(desiredSz(1)/2); % Set figure around center of screen
p2 = round(scrnSz(4)/2)-round(desiredSz(2)/2);
set(f, 'Position', [p1 p2 desiredSz(1) desiredSz(2)-25]);
set(f, 'Name', filename{1});
helpabout = findall(f, 'tag', 'figMenuHelp');

if verLessThan('matlab', '9.3') %matlab R2017b
    uimenu(helpabout, 'Label', '&Help Showsbx', 'Callback', @helpShowsbx);
else
    uimenu(helpabout, 'Text', '&Help Showsbx', 'MenuSelectedFcn', @helpShowsbx);
end

%_______ images for the toolbar
t = uitoolbar(f);
iconcrop = [zeros(3,1); ones(10,1); zeros(3,1)];
iconcrop = cat(2, zeros(16,3), repmat(iconcrop, 1, 10), zeros(16,3));
iconcrop = repmat(iconcrop, 1, 1, 3);

iconsplit = [zeros(4,16); ones(2,16); zeros(4,16); ones(2,16); zeros(4,16)];
iconsplit = repmat(iconsplit, 1, 1, 3);

iconsig = ones(16,16)*0.8;
for i1 = 1:16
    yicon = round((sin(i1/8*pi)+1)*7+1);
    iconsig(yicon:yicon+1, i1) = 0;
    iconsig(8, i1) = 0;
end
iconsig = repmat(iconsig, 1, 1, 3);

%________buttons in toolbar
p = uipushtool(t,'TooltipString','Set Crop',...
    'ClickedCallback',...
    @SetCrop);
p.CData = iconcrop;

p = uipushtool(t,'TooltipString','Show Average of multiple frames',...
    'ClickedCallback',...
    @sbxaverage);

p.CData = imread('sum.png');

p = uipushtool(t,'TooltipString','Split the sbx file into multiple files (for different depths)',...
    'ClickedCallback',...
    {@sbxsplit, info, strfp});

p.CData = iconsplit;

if isempty(Crop)
    tooltipString = 'Align: Starts NoRMCorre motion correction! Remember to set crop first to remove artefacts';
else
    tooltipString = 'Align: Starts NoRMCorre motion correction! Crop present';
end

p = uipushtool(t,'TooltipString', tooltipString,...
    'ClickedCallback',...
    @sbxalign);

p.CData = imread('a.png');

p = uipushtool(t,'TooltipString','Export signal to workspace: draw rectangle, double click in rectangle to confirm',...
    'ClickedCallback',...
    @sbxsigout);
p.CData = iconsig;

p = uipushtool(t,'TooltipString', ...
    'Export mp4', 'ClickedCallback',...
    @sbxmp4out);

p.CData = imread('iconvideo.png');

%_______


im = imagesc(medfilt2(Img(:,:,Chan)), clim); colormap gray

if isstruct(Crop)
    hCropRectangle = rectangle('Position', [Crop.x(1), Crop.y(1), length(Crop.x), length(Crop.y)],...
                               'EdgeColor', 'y', 'LineWidth', 1);
else
    hCropRectangle = gobjects(1);
end
    
%_____ buttons inside the figure itself
% movie position slider
Slider = uicontrol(f,'Style','slider', 'String', 'Slider 1',...
    'Min',1,'Max', framenm, 'SliderStep', [1/(framenm-1) 10/(framenm-1)], ...
    'Value',1, 'Position', [80 10 500 20], 'Callback', @SetPos);

uicontrol(f,'Style', 'pushbutton', 'String', 'Run', 'Position', [10 10 60 20],...
    'Callback', @PlaySbx);

uicontrol(f,'Style', 'pushbutton', 'String', 'Chan', 'Position', [590 10 50 20],...
    'BackGroundColor', 'green', 'Callback', @ChangeChannel);

uicontrol(f,'Style', 'checkbox', 'String', 'median filter', 'Position', [10 40 80 20],...
    'Value', 1, 'Callback', @ToggleFilter);

if Splitnm > 1
    uicontrol(f,'Style', 'listbox', 'String', {int2str((1:Splitnm)')}, 'Position', [20 70 20 Splitnm*15],...
        'Min', 1, 'Max', Splitnm, 'Value', Slice+1, 'Callback', @Setsplit);
    
    uicontrol(f,'Style', 'checkbox', 'String', 'split', 'Position', [10 130 50 20],...
        'Value', bsplit, 'Callback', @Showsplit);
end

cm = uicontextmenu;
im.UIContextMenu = cm;
uimenu(cm,'Label','Change colormap scale','Callback', @ColormapScale);
uimenu(cm,'Label','Save current image as BImg in workspace','Callback', @BImgToWorkspace);
uimenu(cm,'Label','Save current image as tiff','Callback', @savetif);
uimenu(cm,'Label','Save sbx file with Deleted Frame','Callback', @deleteframe);
uimenu(cm,'Label','Save cropped sbx file','Callback', @cropfile);

    function SetPos(source, ~)
        % Jump to specific position/frame in sbx file
        %  global Pos bVers Chan
        Pos = round(source.Value);
        if Pos > framenm
            Slider.Value = framenm;
        else
            Slider.Value = Pos;
        end
        DrawFrame
        title(im.Parent, num2str(Pos))
        drawnow
    end

    function ToggleFilter(~, ~)
        % Toggle median filtering
        bFilt = bFilt ~= 1;
        DrawFrame % redraw the current frame
        drawnow
    end

    function PlaySbx(~, ~) % play movie
        bRun = bRun ~= 1;
        while Pos < framenm && bRun
            DrawFrame
            
            Slider.Value = Slider.Value +1;
            title(im.Parent, num2str(Pos)) % Useful for debugging
            drawnow, pause(0.01)
            Pos = Slider.Value;
        end
    end

    function DrawFrame
        % Reads and draws a frame from the sbx file into the figure
        if bsplit
            Img = sbxread(strfp, (Pos-1)*Splitnm + Slice,1);
        else
            Img = sbxread(strfp, (Pos-1),1);
        end
        if bVers
            Img = permute(Img, [2 3 1]);
        end
        if bFilt
            Img =  medfilt2(Img(:,:,Chan));
        else
            Img = Img(:,:,Chan);
        end
        
        if ishandle(im)
            im.CData = Img;
        else
            figure(f)
            im = imagesc(Img, clim); colormap gray
        end
    end


    function ColormapScale(~,~)
        % Change the colormap scale
        figure(f)
        strClim = inputdlg({'Cutoff low', 'Cutoff high'} , 'Change Colormap scale' , 1,...
                            {num2str(clim(1)), num2str(clim(2))});
        if ~isempty(strClim)
            clim(1) = str2double(strClim{1});
            clim(2) = str2double(strClim{2});
            caxis(clim)
        end
    end


    function Showsplit(source, ~)
        % Toggle to ignore splits or to remain in same split
        bsplit = source.Value;
        info.bsplit = bsplit;
        if bsplit
            framenm = floor(info.max_idx/Splitnm);
        else
            framenm = info.max_idx;
        end
        if framenm < Slider.Value
            Slider.Value = framenm;
        end
        Slider.Max = framenm;
        Slider.SliderStep = [1/(framenm-1) 10/(framenm-1)];
    end


    function Setsplit(source, ~)
        % Set which depth to show
        if bsplit
            Slice = source.Value - 1;
            info.Slice = Slice;
            DrawFrame
            drawnow
        end
    end
    

    function sbxaverage(~,~)
        % Average some number frames to decrease noise a lot
        NmFrames = 100;
        prompt = {'Number of frames to average'};
        dlg_title = 'Average frames';
        num_lines = [1 30];
        def = {'100'};
        answer = inputdlg(prompt,dlg_title, num_lines,def);
        if ~isempty(answer)
            NmFrames = round(str2double(answer{1}));
        end

        if info.bsplit %if splitting has been set on
            Img = sbxread(info.strfp, (Pos-1)*info.Slices, NmFrames*info.Slices);
            if bVers
                Img = squeeze(Img(Chan,:,:, 1+info.Slice:info.Slices:NmFrames*info.Slices));
            else
                Img = squeeze(Img(:,:,Chan, 1+info.Slice:info.Slices:NmFrames*info.Slices));
            end
        else
            Img = sbxread(info.strfp, (Pos-1), NmFrames);
            if bVers
                Img = squeeze(Img(Chan,:,:,:));
            else
                Img = squeeze(Img(:,:,Chan,:));
            end
        end

        Img = uint16(mean(Img,3));
        if ishandle(im)
            im.CData = Img;
        else
            figure(f)
            im = imagesc(Img, clim); colormap gray
        end
    end


    function BImgToWorkspace(~,~)
        % Save the Img to workspace
        BImg = Img - double(clim(1));
        BImg = BImg./double(clim(2)-clim(1));
        assignin('base', 'BImg', BImg)
    end


    function savetif(~,~)
        % Saves the current image as a tif file
        clim = double(clim);
        BImg = double(Img) - clim(1);
        BImg = BImg./(clim(2)-clim(1)).*65535;

        [fn, pn] = uiputfile([info.strfp '_av.tiff']);
        if ischar(fn)
            imwrite(uint16(BImg),[pn fn])
        end
    end

    
    function ChangeChannel(source, ~)
        % Change the color channel
        if lnChan > 1
            if Chan == 1
                Chan = 2;
                set(source, 'BackGroundColor', 'red');
            else
                Chan = 1;
                set(source, 'BackGroundColor', 'green');
            end
        end
        
    end


    function deleteframe(~,~)
        % Save a new sbx file with a frame deleted
        pn = uigetdir('Save file location');
        
        if ischar(pn)
            fn = [pn '\' filename{1} '_copy.sbx'];
            file = fopen(fn, 'w');
            disp('Copying file.........')

            wb = waitbar(0);
            frewind(info.fid);
            for i = 1:(Pos-1)
                data = fread(info.fid,info.nsamples/2,'uint16=>uint16');
                fwrite(file, data, 'uint16');
                waitbar(i/info.max_idx, wb)
            end

            %frame to be skipped
            fread(info.fid,info.nsamples/2,'uint16=>uint16');

            for i = 1:(info.max_idx-Pos)
                data = fread(info.fid, info.nsamples/2,'uint16=>uint16' );
                fwrite(file, data, 'uint16');
                waitbar((i+Pos)/info.max_idx, wb)
            end
            close(wb)
            info.max_idx = info.max_idx -1;

            %also update the frame numbers in the info structure
            if isfield(info, 'frame')
                Mx = max(info.frame);
                ix = find(info.frame == Mx);
                %take care of bug in counting frames
                if ix < length(info.frame)
                    info.frame(ix+1:end) = info.frame(ix+1:end) + 65535;
                end
                %update frame positions now that one frame has been removed
                info.frame(info.frame > Pos) = info.frame(info.frame > Pos) - 1;
            end

            save(fn(1:end-4), 'info');
            fclose(file);
            disp('Done!!!!!!')
        else
            fprintf('Not saving sbx file with a deleted frame because canceled\n')
        end
    end


    function SetCrop(~,~)
        % Set a crop which gets used by other functions like align and cropsbx
        % saved crop from earlier file
        if isfield(info, 'simon') || (isfield(info, 'scanbox_version') &&  info.scanbox_version == 2.5 )
            IWidth = info.Shape(info.Perm(2));
            IHeight = info.Shape(info.Perm(1));
        else
            IWidth = info.sz(2);
            IHeight = info.sz(1);
        end
        
        try
            if ~isempty(Crop) && Crop.x(end) <= IWidth && Crop.y(end) <= IHeight
                org(1) = Crop.x(1);
                org(2) = Crop.y(1);
                wide = length(Crop.x);
                high = length(Crop.y);
                
                info.crop = Crop;
                info.d1 = wide;
                info.d2 = high;
                
                if isgraphics(hCropRectangle)
                    delete(hCropRectangle)
                end
                hCropRectangle = rectangle('Position',[org(1), org(2), wide, high],...
                                           'EdgeColor', 'y');
                                       
                button = questdlg('Renew/Change Crop','Crop', 'Cancel');
                if strcmp(button, 'Cancel') || strcmp(button, 'No')
                    return;
                else
                    bUpCrop = 1;
                end
            end
        catch
            %just ignore if this fails
        end
        
        [x, y] = getpts;
        
        if length(x)>1 && isnumeric(x)
            x = sort(round(x));
            y = sort(round(y));
            if x(1) < 1
                x(1) = 1;
            end
            if x(2) > IWidth
                x(2) = IWidth;
            end
            if y(1) < 1
                y(1) = 1;
            end
            if y(2) > IHeight
                y(2) = IHeight;
            end
            
            w = x(2)-x(1)+1;
            h = y(2)-y(1)+1;
            wide = floor(w/8)*8;
            high = floor(h/8)*8;
            org = round((w-wide)/2+x(1));
            org(2) = round((h-high)/2+y(1));
            
            if isgraphics(hCropRectangle)
                delete(hCropRectangle)
            end
            hCropRectangle = rectangle('Position',[org(1), org(2), wide, high],...
                                       'EdgeColor', 'y');
            
            info.crop.x = (org(1):org(1)+wide-1);
            info.crop.y = (org(2):org(2)+high-1);
            info.d1 = wide;
            info.d2 = high;
            Crop = info.crop;
        end
        
    end

    function cropfile(~,~)
        % Save a new sbx file that's a cropped version of this sbx file
        if ~isempty(Crop) % crop present
            pn = uigetdir('Save file location');
            fn = [pn '\' filename{1} '_cropped.sbx'];

            file = fopen(fn, 'w');
            disp('Copying file.........')

            simon = 0;
            if isfield(info, 'simon')
                simon = info.simon;
            end

            wb = waitbar(0);
            frewind(info.fid);
            for i = 1:info.max_idx
                Y_temp = fread(info.fid, info.nsamples/2,'uint16=>uint16' );
                Y_temp = reshape(Y_temp,info.Shape);
                if ~simon
                    Y_temp = intmax('uint16')- Y_temp;
                end
                if info.scanbox_version == 2.5 || simon
                    if info.nchan == 2
                        Y_temp = Y_temp(Crop.x,Crop.y,:);
                        dim = size(Y_temp);
                        Yt = single(reshape(Y_temp, dim(1), dim(2)*2));
                    else
                        Yt = single(Y_temp(Crop.x,Crop.y,1));
                    end
                else
                    if info.nchan == 2
                        Y_temp = Y_temp(:,Crop.x,Crop.y);
                        Y_temp = permute(Y_temp, [2 3 1]);
                        dim = size(Y_temp);
                        Yt = single(reshape(Y_temp, dim(1), dim(2)*2));
                    else
                        Yt = single(squeeze(Y_temp(1,Crop.x,Crop.y)));
                    end
                end
                fwrite(file, Yt, 'uint16');
                waitbar(i/info.max_idx, wb)
            end
            fclose(file);
            close(wb)

            info.sz = [length(Crop.y),length(Crop.x)];
            info.Shape = [length(Crop.x),length(Crop.y), info.nchan];
            info.Perm = [2 1 3 4];
            info.simon = 1;

            save(fn(1:end-4), 'info');
            disp('Done!!!!!!')
        else % No crop
            msgbox('No crop yet, select a crop first via the button shown in this message',...
                'no crop', 'custom', iconcrop)
        end
    end

    function sbxalign(~,~)
        % Starts motion correction!
        f  = figure('ToolBar','none',...
                    'Position', [100 750 250 180], 'name', 'Align'); %'WindowStyle', 'modal',
        if info.nchan == 2
            ChCntr = uicontrol(f,'Style', 'listbox', 'String', {'Green', 'Both R & G'},...
                'Position', [130 90 100 50],  'Value', 1);
            uicontrol(f,'Style', 'text', 'String', 'Which channels to align',...
                'Position', [20 100 100 40],  'Value', 1);
        end
        
        method = [];
        bg = uibuttongroup('Position',[0 0.78 1 0.2], 'SelectionChangedFcn', @bselection);
        uicontrol(bg,'Style', 'radiobutton', 'String','Rigid', 'Position',[60 10 100 20]);
        uicontrol(bg,'Style', 'radiobutton', 'String','Nonrigid', 'Position',[130 10 100 20]);
        
        function bselection(~,event)
            method = event.NewValue.String;
            % disp(method)
        end
        
        if ~(isfield(info, 'firmware') && str2double(info.firmware) >= 4.1)
            SkipCtrl = uicontrol(f,'Style', 'edit', 'String', '65537',...
                'Position', [130 60 100 20]);
            uicontrol(f,'Style', 'text', 'String', 'Skip frame (or -1)',...
                'Position', [20 60 100 20]);
        end
        
        uicontrol('Position',[130 10 100 30],'String','Continue',...
            'Callback','uiresume(gcbf)');
        uiwait(gcf)
        
        if ~isempty(method)
            info.AlignMethod = method;
        end
        if info.nchan == 2
            info.SelChan = ChCntr.Value; %green or both
        else
            info.SelChan = 1;
        end
        if ~(isfield(info, 'firmware') && str2double(info.firmware) >= 4.1)
            info.Skipframe = str2double(SkipCtrl.String) - 1; %frame indices are 0 based
        end
        close(f)
        
        simonalign3();
    end

    function sbxmp4out(~,~)
        % res = questdlg('Save image sequence to avi file', 'Avi', 'Yes', 'No', 'Yes');
        strFileOut = [strfp '.mp4'];
        [newFn, newPn] = uiputfile('.mp4','Save image sequence to avi file', strFileOut);
        
        if ischar(newPn) && ischar(newFn)
            
            if isfield(info, 'Freq')
                Framerate = info.Freq;
            elseif isfield(info, 'Tframe') %frame time
                Framerate = round(1/info.Tframe);
            else
                if info.scanmode
                    Framerate = round(info.resfreq/512); %unidirectional
                else
                    Framerate = round(info.resfreq/256); %bidirectional
                end
            end
            
            if bsplit
                [~, fname, ~] = fileparts(newFn);
                newFn = [fname '_Split' num2str(Splitnm) '.mp4' ];
                strFileOut = fullfile(newPn, newFn);
                mssg1 = ['Saving split ' num2str(Slice + 1) ' from '...
                          num2str(Splitnm) ' sequences to :' newFn];
                Framerate = round(Framerate/Splitnm);
            else
                strFileOut = fullfile(newPn, newFn);
                mssg1 = ['Saving image sequence to :' newFn];
            end
            mssg2 = ['Framerate = ' num2str(Framerate)];
            
            if bUpCrop
                Crop = info.crop;
                [~, ic] = min(info.Shape);
                if ic == 1
                    xi = Chan;
                    yi = Crop.y;
                    zi = Crop.x;
                else
                    zi = Chan;
                    xi = Crop.y;
                    yi = Crop.x;
                end
            end
            if length(xi) < 62 || length(yi) < 62
                warning('MP4 may be too small to write!, select bigger crop')
            end
            
            outputVideo = VideoWriter(strFileOut, 'MPEG-4');
            outputVideo.FrameRate = Framerate;
            outputVideo.Quality = 100;
            mssg3 = ['Compression: ' outputVideo.VideoCompressionMethod];
            
            open(outputVideo)
            Frame = 1;
            Scl = single(clim);
            Scl(2) = Scl(2)-Scl(1); %scale imagedata between 0 -1.0
            %The videodata is cropped, run through a median filter and scaled
            
            hw = waitbar(0, {mssg1, mssg2, mssg3});
            while Frame < framenm
                if bsplit
                    Im = sbxread(strfp, (Frame-1)*Splitnm + Slice,1);
                else
                    Im = sbxread(strfp, (Frame-1),1);
                end
                if bUpCrop
                    Im = squeeze(Im(xi, yi, zi));
                else
                    Im = squeeze(Im);
                end
                Im = (single(medfilt2(Im))- Scl(1))/Scl(2);
                Im(Im>1) = 1.0;
                Im(Im<0) = 0;
                
                writeVideo(outputVideo,Im)
                Frame = Frame+1;
                waitbar(Frame/framenm, hw);
            end
            close(hw)
            
            close(outputVideo)
            disp('Video saved')
        end
    end

    function helpShowsbx(~,~)
        strfile = which('Showsbx');
        [strpath, strname] = fileparts(strfile);
        html = fullfile(strpath, [strname '.html']);
        
        web(html)
    end
end



