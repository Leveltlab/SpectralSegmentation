function [cent, varargout]  = mylocalmax(d, edg)

%edge ignores border pixels 
sd=size(d);
[x, y]=find(d(edg:sd(1)-edg,edg:sd(2)-edg));

% initialize outputs
cent_map=zeros(sd);

x=x+edg-1;
y=y+edg-1;
for j=1:length(y)
    %area around this position may contain other maxima
    %we need to select the highest.
    wx = x(j)-5:x(j)+5;
    wy = y(j)-5:y(j)+5;
    M = cent_map(wx,wy);
    Pval = max(M(:));
    
    val = d(x(j),y(j));
    if      (val >d(x(j)-1,y(j)-1 ))  &&...
            (val >d(x(j)-1,y(j)))     &&...
            (val >d(x(j)-1,y(j)+1))   &&...
            (val >d(x(j),y(j)-1))     && ...
            (val >d(x(j),y(j)+1))     && ...
            (val >d(x(j)+1,y(j)-1))   && ...
            (val >d(x(j)+1,y(j)))     && ...
            (val >d(x(j)+1,y(j)+1))   && ...
            (val > Pval)
        
       % cent = [cent ;  y(j) ; x(j) ; val];
        cent_map(wx,wy) = 0;
        cent_map(x(j),y(j))=  val; 
        
    end
end

 [xv, yv] = find(cent_map);
 v = arrayfun(@(x,y) cent_map(x,y), xv, yv);
 cent = [xv,yv, v]; %x posistion, y position, max value
 [~, ix] = sort(cent(:,3),'descend');
 cent = cent(ix,:);
     
%cent = (reshape(cent, 3, length(cent)/3))';

if nargout>1 ;  varargout{1}=cent_map; end
