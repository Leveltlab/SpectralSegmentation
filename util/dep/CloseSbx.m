function out = CloseSbx()
% Detects an active sbx file via its global variable and closes it
% 
% Leander de Kraker
% 2026-8-21
% 
global info


if ~isempty(info) & isfield(info, 'fid')
    fclose(info.fid);
    
    clearvars info
    clearvars -global info
    if nargout>0
        out = 'global info detected with file identifier. global info variable deleted and file closed.';
    end
elseif ~isempty(info)
    clearvars info
    clearvars -global info
    if nargout>0
        out = 'global info detected without file identifier. global info deleted.';
    end
elseif nargout>0
    out = 'no global info detected';
end
