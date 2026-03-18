function [s, imS, maxCor] = ShiftLinesCheck(im, trans, hax, doPlot, doPlotAll, shifts)
% Figure out how much a recording's alternating lines should be shifted to
% correct for bidirectional scanning misalignment
% 
% input: 
%    - im (2D, 3D or 4D double): The image where to figure out the shift.
%                                if 3D, please no more than 2 elements.
%    - trans (string): transpose data? (put true if vertical artifact)
%    - hax (axis graphic object): where to plot results which is 
%                                       always generated.
%    - doPlot (boolean): plotting the best found shift
%    - doPlotAll (boolean): plotting the different options of lineshift 
%                                       (default=false)
%    - shiftVal ([1 x n] double): which shifts to check (default -5:5)
% output: s (digid integer): shift in px that the lines should have
% 
% Leander de Kraker
% 2022-12-2
% 2026-3-18: inputs changed and parsed as arguments
% 
%%

arguments
    im
    trans = false;
    hax = gobjects(1);
    doPlot = true;
    doPlotAll = false;
    shifts = -5:5;
end
if isempty(shifts)
    shifts = -5:5;
end

if doPlot
    figure('WindowStyle','docked')
    fprintf('\nplotting the different line shifts.. check for plot\n')
end

% im(im>65500) = median(im(:));

imdim = size(im);
if length(imdim)>3 % if 4D or more, average the 4th dimension (should be frames)
    im = mean(im, 4);
end

if trans
    im = permute(im, [2 1 3]);
    strTrans = 'data is transposed';
else
    strTrans = '';
end

imdim = size(im);
xlims = [1 imdim(2)];
% xlims = [1 300];
ylims = [100 300]; % zoom in
% clims = [2000 20000]; % color limits
clims = [min(im(:)), max(im(:))];


buf = 20;
nshifts = length(shifts);
correl = zeros(nshifts, 2);
correls = zeros(nshifts, 1);

d1 = imdim(1);
if mod(d1, 2) % even: select data to 1 line less
    d1 = d1 - 1;
end


for i = 1:nshifts
    if shifts(i) < 1
        y = 1;
    else
        y = 2;
    end
    x = abs(shifts(i));
    
    imS = im; % image with shifted lines
    imS(y:2:d1, 1:end-x) = imS(y:2:d1, 1+x:end);
    
    imSH = imS(1:2:d1, buf:end-buf);
    imH1 = imS(2:2:d1,  buf:end-buf);
    
    correl(i, :) = [corr2(imH1, imSH), corr2(imH1(2:end), imSH(1:end-1))];
    correls(i) = mean(correl(i, :));
    % PLOTTING
    if doPlotAll
        if size(imS, 3)==2; imS = CreateRGB2_mat(imS, [1 0 0; 0 1 0]);end 
        imagesc(MidtoneBalance(imS, 0.2))
        colormap(cmapL('greenFancy', 256)); colorbar
        title(sprintf('shift %dpx. correl %.5f %s. Press enter to continue checking!',...
                shifts(i), correls(i), strTrans))
        if length(size(imS))==2; clim(clims); end
        xlim(xlims)
        ylim(ylims)
    %     pause(0.05)
        pause()
    end
end
[maxCor, idx] = max(correls);
s = shifts(idx);

% Plotting
if s < 1; y = 1;
else;     y = 2;
end
x = abs(s);
imS = im; % image with shifted lines
imS(y:2:d1, 1:end-x) = imS(y:2:d1, 1+x:end);
if size(imS, 3)==2; imS = CreateRGB2_mat(imS, [1 0 0; 0 1 0]); end


if trans
    imS = permute(imS, [2 1 3]);
end

if doPlot
    if isempty(hax)
        hax = subplot(1,1,1);
    end
    imagesc(MidtoneBalance(imS, 0.2), 'Parent', hax)
    % colormap(hax, cmapL('greenFancy', 256)); colorbar
    title(hax, sprintf('shift %dpx. correl %.5f', s, maxCor))
    if length(size(imS))==2; clim(hax, clims); end
    xlim(hax, xlims)
    ylim(hax, ylims)

    % printing
    fprintf('best shift = %d pixels\n', s)
    if s>0; stri = 'un';
    else;   stri = '';
    end
    if s~=0
        if trans
            strOri = 'columns';
            strDirection = 'up';
        else
            strOri = 'rows';
            strDirection = 'the left';
        end
        fprintf('going to shift %seven %s %d pixels to %s\n', stri, strOri, abs(s), strDirection)
        fprintf('which is %.5f better than no shift\n', maxCor-correls(shifts==0))
    else
        fprintf('no shift recommended I guess')
    end
end
