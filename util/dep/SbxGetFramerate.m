function framerate = SbxGetFramerate(info)
% From Showsbx
% 
% Leander de Kraker, Chris vd Togt
% Netherlands Institute for Neuroscience 
% 
if isfield(info, 'Freq')
    framerate = info.Freq;
elseif isfield(info, 'Tframe') %frame time
    framerate = round(1/info.Tframe);
else
    if info.scanmode
        framerate = round(info.resfreq/512); %unidirectional
    else
        framerate = round(info.resfreq/256); %bidirectional
    end
end
