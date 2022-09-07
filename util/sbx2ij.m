% Chris van der Togt, 2017, 
% Netherlands Institute for Neuroscience 

%Miji
%IJ = ij.IJ;

[fn , pn] = uigetfile('*.sbx');
filename = strsplit(fn, '.');
strfp = [pn filename{1}];
sbxread(strfp, 0,1);
global info

d = dir([strfp '.sbx']); 
lngth =  d.bytes/info.sz(2)/info.sz(1)/info.nchan/2;
if isfield(info, 'Shape') && info.scanbox_version == 2 && info.nchan == 1 && ~isfield(info, 'simon')
   dim = info.Shape([3 2]); 
elseif isfield(info, 'Shape')
    dim = info.Shape([1 2]);
else
    dim = info.sz;
end

if isfield(info, 'simon') %values have been subracted from maxint6
    IJ.run('Raw...', ['open=' strfp '.sbx image=[16-bit Unsigned] width=' num2str(dim(1)) ' height=' num2str(dim(2)) ' number=' num2str(lngth*info.nchan) ' little-endian use']);
elseif info.scanbox_version== 2.5
    IJ.run('Raw...', ['open=' strfp '.sbx image=[16-bit Unsigned] width=' num2str(dim(1)) ' height=' num2str(dim(2)) ' number=' num2str(lngth*info.nchan) ' white little-endian use']);
else
    IJ.run('Raw...', ['open=' strfp '.sbx image=[16-bit Unsigned] width=' num2str(dim(2)) ' height=' num2str(dim(1)) ' number=' num2str(lngth*info.nchan) ' white little-endian use']);
end


%has to be transposed for ImageJ
% IStack = ij.ImageStack(info.sz(2), info.sz(1));
% 
% for i = 0:lngth
%     pixels = sbxread(strfp,i,1);
%     if info.nchan == 2
%         pixels = pixels(:,:,1);   
%     end
%     pixels = pixels'; %transpose
%     ImageProc = ij.process.ShortProcessor( info.sz(2), info.sz(1));
%     ImageProc.setPixels(pixels(:));
%     
%     IStack.addSlice(num2str(i),ImageProc);
% end
% 
% IP = ij.ImagePlus('new', IStack);
% IP.show()
