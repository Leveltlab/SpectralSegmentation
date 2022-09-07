function IM = setminlevel(IM)
% Chris van der Togt, 2017, 
% Netherlands Institute for Neuroscience 

     H = sort(IM(:));
     Mn = find(H > -inf, 1, 'first');
     IM(IM < H(Mn)) = H(Mn);     
     IM = IM - H(Mn);
     
