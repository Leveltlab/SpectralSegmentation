function simonalign3

 global info
 
 disp(info.strfp)
  
%  fids = fopen('all');
%  if any(fids == info.fid)
%     fclose(info.fid);
%  end
 
 if isfield(info, 'bsplit') && info.bsplit == 1 
     disp(['Dimensions : ' num2str([info.sz info.max_idx/info.Slices] ) ', File will be split in ' num2str(info.Slices) ' slices'])
 else
     disp(['Dimensions : ' num2str([info.sz info.max_idx]) ', File will not be split.' ])
 end
 
 if isfield(info, 'crop')
     disp(['Crop x :' num2str(min(info.crop.x)) ':' num2str(max(info.crop.x)) ] ) 
     disp(['Crop y :' num2str(min(info.crop.y)) ':' num2str(max(info.crop.y)) ] ) 
 else
     info.crop.x = (1:info.sz(1));
     info.crop.y = (1:info.sz(2));
     disp(['Crop x :' num2str(min(info.crop.x)) ':' num2str(max(info.crop.x)) ] ) 
     disp(['Crop y :' num2str(min(info.crop.y)) ':' num2str(max(info.crop.y)) ] ) 
 end
 info.d1 = length(info.crop.x);
 info.d2 = length(info.crop.y);
 
if isfield(info, 'SelChan')
    nchan = info.SelChan; %either only green channel or both 1,2
else
    nchan = 1;
end

%% Set parameters 

%set nchan = 1 only green channel; nchan = 2 ; both green and red channel
%nchan = 1;
basefilepath = info.strfp; 

% set parameters (first try out rigid motion correction)
options_rigid = NoRMCorreSetParms('d1',info.d1,'d2',info.d2*nchan, 'd3', 1, 'output_type', 'dat', ...
                                   'bin_width',50,'max_shift',15,'us_fac',50);

options_rigid.scanbox_version = info.scanbox_version;

options_rigid.nsamples = info.nsamples;
options_rigid.Shape = info.Shape;
options_rigid.crop = info.crop;
options_rigid.Perm = info.Perm;
options_rigid.nchan = nchan;

if isfield(info, 'Skipframe') && (info.Skipframe > 0)
    options_rigid.Skipframe = info.Skipframe; %frame bug in scanbox
    info.max_idx = info.max_idx - floor(info.max_idx/info.Skipframe);
end

if isfield(info, 'bsplit') && info.bsplit == 1   
    max_idx = floor(info.max_idx/info.Slices);   
else
    max_idx = info.max_idx;
end
options_rigid.max_idx = max_idx;  %start reading from zero

if isfield(info, 'simon')
    options_rigid.simon = info.simon;  %flip pixel value over intmax not neccessary
end

% save info to mat file
infonew = info;
if nchan == 2
    infonew.channels = 1; %both 1, pmt0 2  pmt1 3
else
    infonew.channels = 2; %only pmt0: green channel
end
infonew.sz = [info.d1, info.d2];
infonew.Shape = [info.d1, info.d2 nchan];
infonew.Perm = [2 1 3 4];
infonew.nsamples = (info.d1 * info.d2 * 2 * nchan);
infonew.simon = 1;
infonew.nchan = nchan;
infonew.Slices = 1; %only one slice in a normcorre file.
if isfield(info, 'Freq')
    infonew.Freq = info.Freq;
else
    infonew.Freq = info.resfreq/info.recordsPerBuffer/info.Slices;
end
infonew.max_idx = max_idx;

infocopy = info;
strfp = [info.strfp '.sbx']; %copied file

%% perform motion correction
 %button = questdlg('Align Rigid or non rigid', 'Align', 'Rigid', 'Nonrigid', 'Rigid');
 if ~isfield(info, 'AlignMethod') || (isfield(info, 'AlignMethod') && strcmp(info.AlignMethod, 'Rigid')) 

    if isfield(infocopy, 'bsplit') && infocopy.bsplit == 1 %to do alignment on unsplit files

        options_rigid.bSplit = infocopy.bsplit;
        options_rigid.Slices = infocopy.Slices;
        for i = 0:infocopy.Slices-1
            options_rigid.Slice = i;
            options_rigid.h5_filename = [basefilepath '_Split' num2str(i+1) '_normcorr.sbx'];
            
            info = infonew;
            info.Slice = i;
            save([basefilepath '_Split' num2str(i+1) '_normcorr.mat'], 'info')
            info = infocopy;
            
            disp(['Aligning (Rigid) ' basefilepath '_Split' num2str(i+1) '_normcorr.sbx'])
            tic; [M1,shifts, template, options] = normcorreSpecSeg(strfp, options_rigid); toc
            
            save([basefilepath '_Split' num2str(i+1) '_Metadata.mat'], 'shifts', 'template', 'options')
        end
    else
        disp(['Aligning (Rigid) ' basefilepath '_normcorr.sbx'])
            options_rigid.h5_filename = [basefilepath '_normcorr.sbx'];
            info = infonew;
            save([basefilepath '_normcorr.mat'], 'info')
            info = infocopy;
            
        tic; [M1,shifts, template, options] = normcorreSpecSeg(strfp, options_rigid); toc
        save([basefilepath '_Metadata.mat'], 'shifts', 'template', 'options')
    end

 elseif (isfield(info, 'AlignMethod') && strcmp(info.AlignMethod, 'Nonrigid'))
%% nonrigid

    grid.x = infocopy.d1/4;
    grid.y = infocopy.d2/4;

    options = options_rigid;
    options.grid_size = [grid.x, grid.y, 1];
    options.mot_uf = [4 4 1];
    
    if isfield(infocopy, 'bsplit') %to do alignment on unsplit files
        options.bSplit = infocopy.bsplit;
        options.Slices = infocopy.Slices;
        for i = 0:infocopy.Slices-1
            options.Slice = infocopy.Slice;
            options.h5_filename = [basefilepath '_Split' num2str(i+1) '_normcorr.sbx'];
            
            info = infonew;
            info.Slice = i;
            save([basefilepath '_Split' num2str(i+1) '_normcorr.mat'], 'info')
            info = infocopy;
            
            disp(['Aligning (Nonrigid) ' basefilepath '_Split' num2str(i+1) '_normcorr.sbx'])
            tic; [M1,shifts, template, options] = normcorreSpecSeg(strfp ,options); toc
            
            save([basefilepath '_Split' num2str(i+1) '_Metadata.mat'], 'shifts', 'template', 'options')
        end
    else          
        disp(['Aligning (Nonrigid) ' basefilepath '_normcorr.sbx'])
        options.h5_filename = [basefilepath '_normcorr.sbx'];
        info = infonew;
        save([basefilepath '_normcorr.mat'], 'info')
        info = infocopy;
            
        tic; [M1,shifts,template, options] = normcorreSpecSeg(strfp ,options); toc
        save([basefilepath '_Metadata.mat'], 'shifts', 'template', 'options')
    end
 end
 
 



