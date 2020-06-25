function roundedness = perimarea(x, y)
%calculate area to perimeter of a contour
% roundedness = 1 for perfect circle

fac = round(length(x)/50);
if fac > 0
    p = smoothG([x(:) y(:)], fac);
    p = [p(end,:); p];
else
    p = [x(:) y(:)];
end
dx = diff(p(:,1));
dy = diff(p(:,2));
A = polyarea(p(:,1), p(:,2));
P = sum(sqrt(dx.^2 + dy.^2));
roundedness = 4*pi*A/(P*P);
%roundedness = sqrt(A)/P;