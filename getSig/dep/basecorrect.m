function sig = basecorrect(sig, window)
% Convert fluorescence signal to Df/f
% Chris vd Togt.
% 2023-2-22. Leander changed calculation from (F-F0)/linearfitF0 + 1 
%                                          to (F-F0)/F0

% 10th quantile interpolated baseline correction
f0 = prctfilt(sig',10,window,window,0)';  
sig = (sig - f0) ./ f0;