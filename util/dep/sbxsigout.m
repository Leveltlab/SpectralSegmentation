function sbxsigout(~,~)

global info
persistent Sbx filenm

h = imrect();
pos = round(wait(h));
x = pos(1):pos(1)+pos(3);
y = pos(2):pos(2)+pos(4);

dim = [info.Shape, info.max_idx];

disp('Memory mapping file.....')
if isempty(Sbx) || ~strcmp(filenm, info.strfp)
    Sbx = memmapfile([info.strfp '.sbx'], 'Format',{'uint16', dim, 'Y'}, 'Writable', true);
    filenm = info.strfp;
end

disp('Retrieving data....')
if dim(1) < 3
    Sig = squeeze(mean(mean(Sbx.Data.Y(:,x,y,:),3),2));
else
    Sig = squeeze(mean(mean(Sbx.Data.Y(x,y,:,:),2)));
end

if ~isfield(info, 'simon')
    Sig = double(intmax('uint16'))- Sig;
end
if size(Sig,1) < size(Sig,2)
    Sig = Sig'; %put signal on first dim
end

if info.bsplit 
    disp('Splitting signals....')
    S = cell(info.Slices,1);
    for i = 1:info.Slices
        S{i} = Sig(i:info.Slices:end,:);
    end
    Sig = S;
end
disp('Done')

% tic
% Sig1 = zeros(info.max_idx,1);
% for i = 0:info.max_idx-1   
%     Img = sbxread(fn, i,1);
%     if dim(1) < 3
%         Sig1(i+1) = mean(mean(Img(:,y,x),3),2);
%     else
%         Sig1(i+1) = mean(mean(Img(y,x,:),2));
%     end
% end
% toc

assignin('base', 'Sig', Sig)

delete(h)

