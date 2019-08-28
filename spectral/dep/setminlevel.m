function IM = setminlevel(IM)

     H = sort(IM(:));
     Mn = find(H > -inf, 1, 'first');
     IM(IM < H(Mn)) = H(Mn);     
     IM = IM - H(Mn);
     
