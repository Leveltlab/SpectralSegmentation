function cent  = mylocalmax(d, edge)

%edge ignores border pixels 
dsz=size(d);
[x, y]=find(d(edge:dsz(1)-edge,edge:dsz(2)-edge));

% initialize outputs
cent_map=zeros(dsz);

x=x+edge-1;
y=y+edge-1;
dx = -4:4;
for j=1:length(y)
    %area around this position may contain other maxima
    %we need to select the highest.
    wx = x(j)+ dx;
    wy = y(j)+ dx;
    M = cent_map(wx,wy);
    Pval = max(M(:));
    
    val = d(x(j),y(j));
    if  (val >d(x(j)-1,y(j)-1 ))  &&...
        (val >d(x(j)-1,y(j)))     &&...
        (val >d(x(j)-1,y(j)+1))   &&...
        (val >d(x(j),y(j)-1))     &&...
        (val >d(x(j),y(j)+1))     &&...
        (val >d(x(j)+1,y(j)-1))   &&...
        (val >d(x(j)+1,y(j)))     && ...
        (val >d(x(j)+1,y(j)+1))   && ...
        (val >= Pval)           
            
        
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
