function [Con, A, F, Pin, Ro] = getCon(Imgin, th, area, Rof, py, px, Iy, Ix)

dim = size(Imgin);
F = zeros(dim);
%    M = -Inf;
%    P = 0;

c = contourc(Imgin,[th th]);
s = getcontourlines(c);
v = arrayfun(@(x) eq(x.c,1),s); %closed
iv = find(v);
A = 0;
Ro = 1;
Pin = false; %is there a contour with peak inside
Con = [];
if ~isempty(iv)
    A = area(1); 
end
for j = 1:length(iv)
    vx = s(iv(j)).x;
    vy = s(iv(j)).y;
  % plot(vx, vy, 'r')
    In = find(inpolygon(px,py,vx,vy),1);

    if ~isempty(In) 
        Pin = true;
        At = polyarea(vx,vy);
        if At > area(1) 
            A = At;
            Ro = perimarea(vx, vy);
%            indices = poly2mask(vx,vy, dim(1), dim(2));
%             MxA = max(Imgin(indices));
%             if M < MxA
%                M = MxA;
%             end
            if At < area(2) && Ro > Rof
                Con.x = vx;
                Con.y = vy;            
                F(inpolygon(Ix,Iy,vx,vy)) = 1;
                %pixv = Imgin(F>0);                                     
            end
        end
        return;
    end
end
