function [Con, A, F, Pin, roundedness] = getCon(Imgin, th, area, Rof, py, px, Iy, Ix)
% Chris van der Togt, 2017
% Netherlands Institute for Neuroscience 


dim = size(Imgin);
F = zeros(dim);
%    M = -Inf;

c = contourc(Imgin,[th th]);
s = getcontourlines(c);
v = arrayfun(@(x) eq(x.c,1),s); %closed
iv = find(v);
A = 0;
roundedness = 1;
Pin = false; %is there a contour with peak inside
Con = [];
if ~isempty(iv)
    A = area(1); 
end
for j = 1:length(iv)
    vx = s(iv(j)).x;
    vy = s(iv(j)).y;
    In = find(inpolygon(px,py,vx,vy),1);

    if ~isempty(In) 
        Pin = true;
        At = polyarea(vx,vy);
        if At > area(1) 
            A = At;
            roundedness = perimarea(vx, vy);
            if At < area(2) && roundedness > Rof
                Con.x = vx;
                Con.y = vy;            
                F(inpolygon(Ix,Iy,vx,vy)) = 1;
            end
        end
        return;
    end
end
