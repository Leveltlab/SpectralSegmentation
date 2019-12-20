function sbxaverage(~,~)

global Pos clim Chan info

if isempty(Pos)
    return;
end

bVers = info.bVers;
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

Img = mean(Img,3);
im = imagesc(Img); colormap gray, caxis(clim);

cm = uicontextmenu;
im.UIContextMenu = cm;
uimenu(cm,'Label','Save image (tiff)','Callback', {@savetif, Img, clim, info.strfp});
uimenu(cm,'Label','Create BImg in workspace','Callback', {@bimg, Img, clim});

function savetif(~,~,Img, Clim, strfp)


    Clim = double(Clim);
    Img = Img - Clim(1);
    Img = Img./(Clim(2)-Clim(1)).*65535;

    [fn, pn] = uiputfile([strfp '_av.tiff']);
    if ischar(fn)
        imwrite(uint16(Img),[pn fn])
    end
    

    
function bimg(~,~, Img, Clim)
    Clim = double(Clim);
    Img = Img - Clim(1);
    Img = Img./(Clim(2)-Clim(1));
    
    assignin('base', 'BImg', Img)
    
    