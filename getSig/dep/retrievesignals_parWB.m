function retrievesignals_parWB(~)
% Waitbar (WB) for use in a parfor loop in DeconvolveSignals.m
global p N h
waitbar(p/N, h, sprintf('Running MLspike (Runtime: %1.2f min)',toc/60));
p = p + 1;
end