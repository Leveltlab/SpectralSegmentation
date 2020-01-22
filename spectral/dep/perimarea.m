function F = perimarea(x, y)
%calculate area to perimeter of a contour
%F = 1 for perfect circle

fac = round(length(x)/50);
if fac > 0
    p = smoothG([x(:) y(:)], fac);
else
    p = [x(:) y(:)];
end
dx = diff(p(:,1));
dy = diff(p(:,2));
A = polyarea(x, y);
P = sum(sqrt(dx.^2 + dy.^2));
F = 4*pi*A/(P*P);
%F = sqrt(A)/P;