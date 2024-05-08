function SbxCropTime(strPathIn, strSbxIn, strPathOut, frameStart, frameEnd)
% SbxCropTime(strPathIn, strSbxIn, strPathOut, frameStart, frameEnd)
% Cropping an sbx file in time, creating a copy of it.
% Converting it from a long file, to a shorter file.
% 
% Input:
%   - strPathIn (string) (ending with '\')
%   - strSbxIn (string) filename of sbx
%   - strPathOut (string) (folder to put output in, ending with '\')
%   - frameStart (scalar double): at which frame to start reading data
%   - frameEnd (scalar double): final frame to read data from
%  
% Leander de Kraker
% 2024-5-7
% 
%%

strSbxIn = strsplit(strSbxIn, {'.mat', '.sbx'});
strSbxIn = strSbxIn{1};

strIn  = [strPathIn,  strSbxIn];
strOut = [strPathOut, strSbxIn, '_cropT'];

nframes = frameEnd-frameStart+1;

fileID = fopen([strOut, '.sbx'], 'w'); % open sbx file
ttoreport = max(100, round(nframes/10));

clearvars -global info
clearvars info
load(strIn, 'info')

% Writing the data
tic
for i = 1:nframes
    frame = sbxread(strIn, frameStart+i-2, 1);
    fwrite(fileID, intmax('uint16')-frame', 'uint16');
    if mod(i,ttoreport)==0
        fprintf('done with frame %d/%d  (%.0f%%)  ETA %.1f minutes\n',...
                         i, nframes, (i/nframes*100), ETA(i, nframes, toc)/60)
    end
end
fprintf('Done\n')
fclose(fileID);

% Adjusting the info approperiately
if isfield(info, 'frame')
    info.frame = info.frame - frameStart + 1;
    toDel = info.frame<1 | info.frame>nframes;
    info.frame(toDel) = [];
    if isfield(info, 'line')
        info.line(toDel) = [];
    end
    if isfield(info, 'event_id')
        info.event_id(toDel) = [];
    end
end
info.max_idx = nframes;
rmfield(info, 'bytesPerBuffer');
save([strOut, '.mat'], 'info')



