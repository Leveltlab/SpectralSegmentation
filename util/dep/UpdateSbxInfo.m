function [info, nSlices] = UpdateSbxInfo(info, pn, fn)
% sbx info can be unreliable with missing fields.
% This function will:
%       set or check strfp (string filepath + name)
%       set Slices (different imaged depths in the frames dimension)
%       set bSplit (boolean saying if multiple slices are present)
%       set bVers (a boolean which says if the data is permuted because of
%                   GPU recording)
% 
% Leander de Kraker
% 2026-3-26
% 

% set the name of the recording
if ~isfield(info, 'strfp')
    info.strfp = fullfile(pn, fn);
elseif ~strcmp(info.strfp, [pn fn]) % inconsitent names, maybe it'll all go wrong
    info.strfp
    [pn fn]
    fprintf('expected name and name in info are inconsistent\nPROBLEM\n\n')
else
    fprintf('everything is going well\n')
end

%number of splits
if ~isfield(info, 'Slices')
    if isfield(info, 'otparam') && length(info.otparam)>=3
        info.Slices = info.otparam(3);
    else
        %  % % code commented because if it is bad, probably 1 depth anyway % % % 
        % !Maybe the info file is super bad, ask how many slices there are!
        % if ~isfield(info, 'Slices')
        %     answ = inputdlg('How many Slices?', 'Slice info!!...', [1 40],  {'1'});
        %     info.Slices = str2double(answ{1});
        % end
        
        % The info file could have no slices because only 1 depth is imaged
        info.Slices = 1;
    end
end

nSlices = info.Slices;
if ~isfield(info, 'bsplit')
    if nSlices > 1 % bsplit true means there are more than 1 split
        info.bsplit = true;
    else
        info.bsplit = false;
    end
end

if size(Img,1) < 3
    info.bVers = true; % bvers means data is permuted because of GPU recording..
else
    info.bVers = false;
end