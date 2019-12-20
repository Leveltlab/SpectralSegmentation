function x = sbxread(fname,k,N,varargin)

% img = sbxread(fname,k,N,varargin)
%
% Reads from frame k to k+N-1 in file fname
% 
% fname - the file name (e.g., 'xx0_000_001')
% k     - the index of the first frame to be read.  The first index is 0.
% N     - the number of consecutive frames to read starting with k.
%
% If N>1 it returns a 4D array of size = [#pmt rows cols N] 
% If N=1 it returns a 3D array of size = [#pmt rows cols]
%
%   gpu version [rows cols #pmt N]
%
% #pmts is the number of pmt channels being sampled (1 or 2)
% rows is the number of lines in the image
% cols is the number of pixels in each line
%
%
% The function also creates a global 'info' variable with additional
% informationi about the file

global info

persistent Perm Shape


% check if already loaded...
if isempty(info) || (~isempty(info) && ~strcmp(fname, info.strfp))
    if ~isempty(info)
        try
            fclose(info.fid);       
        catch
        end
        info = [];
    end

    load(fname);
    info.strfp = fname;
    
%     if(exist([fname ,'.align'])) % aligned?
%         info.aligned = load([fname ,'.align'],'-mat');
%     else
%         info.aligned = [];
%     end   

    if(~isfield(info,'sz'))
        info.sz = [512 796];    % it was only sz = .... 
    end
    
    switch info.channels
        case 1
            info.nchan = 2;      % both PMT0 & 1
        case 2
            info.nchan = 1;      % PMT 0
        case 3
            info.nchan = 1;      % PMT 1
    end
    

    if ~isfield(info,'Perm')
        Perm = [1 3 2 4]; %1st dim is channel
        Shape = [info.nchan info.sz(2) info.sz(1)];
        if info.scanbox_version == 2.5 || isfield(info, 'simon') %gpu or processed
            Perm = [2 1 3 4]; %3rd dim is channel
            Shape = [info.sz(2) info.sz(1) info.nchan];
        end        
        info.Perm = Perm;
        info.Shape = Shape;
    else
        if length(info.Perm) == 3 %there is a bug somewhere!
            disp('info.Perm was incorrect!!')
            info.Perm = [info.Perm 4];
%            save(fname, 'info');
        end
        Perm = info.Perm;
        Shape = info.Shape;
    end
    
    info.fid = fopen([fname '.sbx']);
    d = dir([fname '.sbx']); 
    info.max_idx =  round(d.bytes/info.sz(2)/info.sz(1)/info.nchan/2);
    info.nsamples = (info.sz(2) * info.sz(1) * 2 * info.nchan);   % bytes per frame 
    
elseif isempty(Perm) || isempty(Shape)
    Perm = info.Perm;
    Shape = info.Shape;
end

if(isfield(info,'fid') && info.fid ~= -1 && k + N > 0)
    
    % nsamples = info.postTriggerSamples * info.recordsPerBuffer;
    try
        fseek(info.fid,k*info.nsamples,'bof');
        x = fread(info.fid,info.nsamples/2 * N,'uint16=>uint16');        
        %changed for gpu version
        x = reshape(x,[Shape N]);
            
    catch
        error('Cannot read frame.  Index range likely outside of bounds.');
    end
    %changed for gpu version
        if isfield(info, 'simon')
            x = permute(x,Perm);
        else
            x = intmax('uint16')-permute(x,Perm);  
        end
     
else
    x = [];
end