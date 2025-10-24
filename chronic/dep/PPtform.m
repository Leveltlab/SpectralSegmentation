function PP = PPtform(PP, tform, offset)
% Update coordinates using transform
% 
% 
% 
% 
% Leander de Kraker
% 2025-10-22
% 
arguments
    PP struct
    tform
    offset double = [0 0];
end

[PP.P(1,:), PP.P(2,:)] = transformPointsForward(tform, PP.P(1,:), PP.P(2,:));
PP.P(1,:) = PP.P(1,:) + offset(1);
PP.P(2,:) = PP.P(2,:) + offset(2);
for i = 1:PP.Cnt
    [PP.Con(i).x, PP.Con(i).y] = transformPointsForward(tform, PP.Con(i).x, PP.Con(i).y);
    PP.Con(i).x = PP.Con(i).x + offset(1);
    PP.Con(i).y = PP.Con(i).y + offset(2);
end