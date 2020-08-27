function sigb = basecorrect(sig, window)


bsl = prctfilt(sig',10,window,window,0)';  % 10% = baseline percentile correction
x = (1:size(bsl,1))';

%scale roi traces to DF/F
parfor i = 1:size(bsl,2)
    p = polyfit(x,bsl(:,i),1); %linear baseline estimator: 1st deg polynomal fit = linear least squares fit
    sigb(:,i) = (sig(:,i) - bsl(:,i))./polyval(p,x) + 1.0;
end

