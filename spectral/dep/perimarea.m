function F = perimarea(x, y)
%calculate area to perimeter of a contour
%F = 1 for perfect circle

dx = diff(x);
dy = diff(y);
A = polyarea(x, y);
P = sum(sqrt(dx.^2 + dy.^2));
F = 4*pi*A/(P*P);
%F = sqrt(A)/P;