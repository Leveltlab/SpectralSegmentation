function sbxsplit(~,~,info, fn)
%
    
%[fn, pn] = uigetfile('*.mat');
% sbxread(fn, 0, 1); %gets info and intializes other parameters

def = {'2', num2str(info.max_idx)};
if isfield(info, 'otparam')
    try
        def{1} = num2str(info.otparam(3));
        disp(['This stack contains ' num2str(info.otparam(3)) ' depth slices'])
    catch
        def{1} = '1';
    end
end
prompt = {'Slices', 'Stacklength'};
dlg_title = 'Split Sbx file in slices:';
num_lines = [1 50];

answer = inputdlg(prompt,dlg_title,num_lines,def);
if isempty(answer)
    return; %cancel button was pressed
end

nmslices = str2double(answer{1});
mxidx = str2double(answer{2});


stackln = floor(mxidx / nmslices);
N = info.nchan;
info.Slices = nmslices;

savefl = zeros(nmslices,1);
for i = 1:nmslices
    savefl(i) = fopen([fn '_split' num2str(i)  '.sbx'], 'a');
end


if(isfield(info,'fid') && info.fid ~= -1)
    strfn = strrep(fn, '\', '\\');
    hl = waitbar(1/stackln, ['Slicing ' strfn]);
    for i = 1:stackln   %the first image of the stack is always bad
        for j = 1:nmslices
            try  %nsamples = number of bytes
                fseek(info.fid,(i*nmslices+j-1)*info.nsamples,'bof');
                x = fread(info.fid,info.nsamples/2,'uint16=>uint16');

            catch
                error('Cannot read frame.  Index range likely outside of bounds.');
            end
            
            fwrite(savefl(j), x, 'uint16');
        end
        waitbar(i/stackln, hl);
    end
end
close(hl)

info.max_idx = stackln-1;
for i = 1:nmslices
    fclose(savefl(i));
    info.Section = i; %which section was this
    save([fn '_split' num2str(i) ], 'info' )
end

end

