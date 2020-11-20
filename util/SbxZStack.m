% Visualize an sbx z stack recording
% Run only one of the two reading options: OPTOTUNE or KNOBBY
%
% Leander de Kraker
% 2018-10-15

[fn, pn] = uigetfile('%.mat', 'get info file of sbx z stack file');
filename = [pn fn];
load(filename)

%% OPTOTUNE  COMPATIBLE
% Optotune stacks are created by going through the different depths of the
% stack multiple times. 
% 
% 
% Number of splits
if ~isfield(info, 'Slices')
    try
        ndepths = info.otparam(3);
    catch
        answ = inputdlg('How many Slices?', 'Slice info!!...', [1 40],  {'1'});
        ndepths = str2double(answ{1});
    end
end
% Number of frames
if isfield(info, 'max_idx')
    nframes = info.max_idx;
elseif isfield(info, 'config')
    if isfield(info.config, 'frames')
        nframes = info.config.frames;
    end
end

% Reading the file and seperating red and green
data = sbxread(filename(1:end-4), 0, nframes);
dataRed = squeeze(data(2,:,:,:));
dataGre = squeeze(data(1,:,:,:));
dims = size(dataGre);

% Averaging the different frames with each other
nFA = floor(dims(3)/ndepths); % number of Frames Actually possible
if nFA*ndepths ~= nframes
    fprintf('Ask Prem wtf\n')
end
% Averaging the depths by reshaping the data so each frame gets put in a
% 4th dimension
dataRed = mean(reshape(dataRed(:,:,1:nFA*ndepths), dims(1), dims(2),ndepths ,nFA),4);
dataGre = mean(reshape(dataGre(:,:,1:nFA*ndepths), dims(1), dims(2),ndepths, nFA),4);

dataRed2 = dataRed;
lim1 = prctile(dataRed2(:),10);
lim2 = 20534; % Set the maxed out lines to lowest value for visibility
dataRed2(dataRed2<lim1 | dataRed2>lim2) = lim1;
dataRed2 = (dataRed2-min(dataRed2(:)))/range(dataRed2(:));
dataGre2 = dataGre;
lim1 = prctile(dataGre2(:),10);
lim2 = 20534;
dataGre2(dataGre2<lim1 | dataGre2>lim2) = lim1;
dataGre2 = (dataGre2-min(dataGre2(:)))/range(dataGre2(:));

zoom = 1.6; % FROM LOGBOOK
scalexy = (1000/zoom) / dims(2); % µm images / number of pixels DEPENDS ON ZOOM

x = [1:dims(2)]*scalexy; % The position of pixels in um
y = [1:dims(1)]*scalexy;
z = info.otwave_um; % [1:dims(3)]*scalexz % scale 


%% KNOBBY
% Knobby stacks are created by doing imaging the same depth x amount of
% times and then going to the next depth.
% 
% As far as I know almost all the info in the info files is wrong

nreps = 1000; % number of frames per depth (repetitions)
nframes = info.config.frames;
ndepths = nframes/nreps;
depthStart = 0:nreps:nframes-1;
fprintf('number of depths: %.1f\n',ndepths)

exampleFrame = sbxread(filename(1:end-4), 0, 3);
dims = size(exampleFrame);
nchannels = dims(1);
heig = dims(1);
wid = dims(2);
if nchannels ~= 2
    warning('There are %d color channels. This code would probably crash.', nchannels)
    error('There need to be 2 color channels')
end

dataBuffer = zeros(nchannels, heig, wid, nreps);% this will store nreps frames every depth
dataGre = zeros(heig, wid, ndepths); % 3D matrix for the zstack of channel 1
dataRed = zeros(heig, wid, ndepths); % 3D matrix for the zstack of channel 2

% Read through the sbx file and average the frames for every depth
tic
for i = 1:ndepths
    dataBuffer = sbxread(filename(1:end-4), depthStart(i), nreps);
    dataGre(:,:,i) = squeeze(mean(dataBuffer(1, :, :, :),4));
    dataRed(:,:,i) = squeeze(mean(dataBuffer(2, :, :, :),4));
    if mod(i,5)==0
        fprintf('processing depth %d/%d, elapsed time %.1f minutes\n', i, ndepths, toc/60)
    end
end

dataRed2 = dataRed;
lim1 = prctile(dataRed2(:),50);
lim2 = 65534;
dataRed2(dataRed2<lim1 | dataRed2>lim2) = lim1;
dataRed2 = dataRed2-min(dataRed2(:))/range(dataRed2(:));
dataGre2 = dataGre;
lim1 = prctile(dataGre2(:),3);
lim2 = 65534;
dataGre2(dataGre2<lim1 | dataGre2>lim2) = lim1;
dataGre2 = dataGre2-min(dataGre2(:))/range(dataGre2(:));

zoom = 1.6;
standardSize = 1000; % According to the website 1mm->1000µm is the standard FOV
scalexy = (1000/zoom) / wid; % µm images / number of pixels | DEPENDS ON ZOOM

x = [1:wid]*scalexy; % The position of pixels in um
y = [1:heig]*scalexy;
z = 1:5:5*92; % z scale | INFO FROM LOGBOOK



%% Plot the Z-stack

m = 0.05; % margin to the sides of the figure (relative)
ys = double(max(y)); % dimensions of the stack in um
xs = double(max(x));
zs = double(max(z));
h = (ys + zs);
w = (xs + zs); % Total width of the image when they are plotted next to each other

% position = [horstart, vertstart, width, height] (start at bottomleft)
position1 = [m,    ys/h, xs/w-m, zs/h-m]; % The first vertical slice on top of the x direction
position2 = [m,    m,    xs/w-m, ys/h-m]; % The image as we usually see it.
position3 = [xs/w, m,    zs/w-m, ys/h-m]; % The vertical slice 90 degrees rot of 1st one

% Plot the figure with both red and green colorchannel____________________
clearvars subplot % No tight subplotting please. But manual subplot coordinates
handvat = figure();
handvat.Units = 'pixels';
m2 = 100; % A margin for where to plot the figure (pixels)
screenSize = get(0,'screensize'); % Set the aspect ratio of the figure correctly!
if screenSize(3)<(w+m2) || screenSize(4)<(h+m2) % reduce size of screen if fig too large
    scaleImg = min((screenSize(3)/(w+m2*2)), (screenSize(4)/(h+m2*2)));
    handvat.InnerPosition = [m2, m2, w*scaleImg h*scaleImg];
else
    handvat.InnerPosition = [m2 m2 w h];
end
subplot('position', position1)
greXZ = log1p(squeeze(max(dataGre2,[],1)))';
redXZ = log1p(squeeze(max(dataRed2,[],1)))';
RGBXZ = CreateRGB({redXZ, greXZ}, 1:2, 'r g');
imagesc(x,z,RGBXZ)
title('max projections of zstack, axis in µm')
subplot('position', position2)
greXY = log1p(squeeze(max(dataGre2,[],3)));
redXY = log1p(squeeze(max(dataRed2,[],3)));
RGBXY = CreateRGB({redXY, greXY}, 1:2, 'r g');
imagesc(x,y,RGBXY)
subplot('position', position3)
greZY = log1p(squeeze(max(dataGre2,[],2)));
redZY = log1p(squeeze(max(dataRed2,[],2)));
RGBZY = CreateRGB({redZY, greZY}, 1:2, 'r g');
imagesc(z,y,RGBZY)

% Plot only the red channel__________________________________________
handvat = figure();
handvat.Units = 'pixels';
m2 = 100; % A margin for where to plot the figure (pixels)
screenSize = get(0,'screensize'); % Set the aspect ratio of the figure correctly!
if screenSize(3)<(w+m2) || screenSize(4)<(h+m2) % reduce size of screen if fig too large
    scaleImg = min((screenSize(3)/(w+m2*2)), (screenSize(4)/(h+m2*2)));
    handvat.InnerPosition = [m2, m2, w*scaleImg h*scaleImg];
else
    handvat.InnerPosition = [m2 m2 w h];
end
subplot('position', position1)
imagesc(x,z,redXZ)
title('max projections of zstack, axis in µm')
subplot('position', position2)
imagesc(x,y,redXY)
subplot('position', position3)
imagesc(z,y,redZY)
colormap(cmapL([1 1 1; 1 1 0; 1 0 0; 0.3 0 0; 0 0 0],256))

% Plot only the green channel__________________________________________
handvatSub = gobjects(3,1);
handvat = figure();
handvat.Units = 'pixels';
m2 = 100; % A margin for where to plot the figure (pixels)
screenSize = get(0,'screensize'); % Set the aspect ratio of the figure correctly!
if screenSize(3)<(w+m2) || screenSize(4)<(h+m2) % reduce size of screen if fig too large
    scaleImg = min((screenSize(3)/(w+m2*2)), (screenSize(4)/(h+m2*2)));
    handvat.InnerPosition = [m2, m2, w*scaleImg h*scaleImg];
else
    handvat.InnerPosition = [m2 m2 w h];
end
handvatSub(1) = subplot('position', position1);
imagesc(x,z,greXZ)
title('max projections of zstack, axis in µm')
handvatSub(2) = subplot('position', position2);
imagesc(x,y,greXY)
handvatSub(3) = subplot('position', position3);
imagesc(z,y,greZY)
colormap(cmapL([1 1 1; 0 0.7 0; 0 0 0],256))
linkaxes([handvatSub(1) handvatSub(2)], 'x')
linkaxes([handvatSub(2) handvatSub(3)], 'y')

%% Save figure

[savename, savepath] = uiputfile();
oldpath = cd;
cd(savepath)
print(gcf, '-dpng',[savename '.png'], '-r300')
cd(oldpath) % navigate back to original folder

%% Colorful depth image like in Glas, Goltstein 2019: miniaturized 2P imaging in mouse V

dataGre2Cell = squeeze(squeeze(num2cell(dataGre2,[1 2])));
dataRed2Cell = squeeze(squeeze(num2cell(dataRed2,[1 2])));
colorsDep = cmapL([1 0 0; 1 1 0; 0 1 0; 0 1 1; 0.3 0.3 1; 0 0 1], ndepths);

figure
imagesc(CreateRGB2(dataGre2Cell, 1:ndepths, colorsDep))
colormap(colorsDep)
title('green Z stack colored by depth')
hand = colorbar;
title(hand, 'depth')

figure
imagesc(CreateRGB2(dataRed2Cell, 1:ndepths, colorsDep))
title('red Z stack colored by depth')
colormap(colorsDep)
hand = colorbar;
title(hand, 'depth')


%% plot 3 slices of the data
depthsQuest = [80 160 250]; % Requested depths to plot
index = zeros(1,3);
for i = 1:3 % Find closest imaged plane
    [~, index(i)] = min(abs(double(z)-depthsQuest(i)));
end
depthmu = double(z(index));

subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0], [0.05 0.05], [0.05 0.05]);

handvat = gobjects(1,5);
figure
handvat(1) = subplot(3,2,[1 3]);
imagesc(x,z,RGBXZ); hold on
for i = 1:3
    plot([x(1) x(end)], [depthmu(i) depthmu(i)], 'y', 'linewidth', 2)
end
imgDepth1 = dataGre2(:,:,index(1));
imgDepth2 = dataGre2(:,:,index(2));
imgDepth3 = dataGre2(:,:,index(3));
handvat(2) = subplot(3,2,5);
imagesc(x,y,CreateRGB({imgDepth1 imgDepth2 imgDepth3}, 1:3, 'r g b'))
title(sprintf('red = %dµm, green=%dµm, blue=%dµm', depthmu(1),depthmu(2),depthmu(3)))
handvat(3) = subplot(3, 2, 2);
imagesc(x,y,imgDepth1)
title(sprintf('depth %d um', depthmu(1)))
handvat(4) = subplot(3, 2, 4);
imagesc(x,y,imgDepth2)
title(sprintf('depth %d um', depthmu(2)))
handvat(5) = subplot(3, 2, 6);
imagesc(x,y,imgDepth3)
title(sprintf('depth %d um', depthmu(3)))
colormap(cmapL([1 1 1; 0 0.8 0; 0 0 0],256))
linkaxes(handvat,'x')


%% Walk through the Zstack
% With the tool made by Tobias v.d. Bijl

viewVolume(log1p(dataGre2))

%% Walk through the Z stack
figure
for i = 1:size(dataGre,3)
    imagesc(dataGre(:,111:end,i));
    title(sprintf('%d',i))
    pause()
end
%% Register the two different stacks
% Lots of registering
%
% Backup the dataGre etc
dataGre2Backup = dataGre2;
dataRed2Backup = dataRed2;

%%
% First register the knobby stack because the recording lasted very long
% Register by registering adjecent depths
dataGre2Regi = dataGre2;
dataRed2Regi = dataRed2;
for i = 2:ndepths
    dataGre2Regi(:,:,i) = Register2Imgs(dataGre2Regi(:,:,i-1), dataGre2Regi(:,:,i));
    dataRed2Regi(:,:,i) = Register2Imgs(dataRed2Regi(:,:,i-1), dataRed2Regi(:,:,i));
    if mod(i,5)==0
        fprintf('registering depth %d\n',i)
    end
end
fprintf('done!\n')

dataGre2 = dataGre2Regi;
dataRed2 = dataRed2Regi;

% Luminocity normalization per depth
meansGre = zeros(1,ndepths);
meansRed = zeros(1,ndepths);
for i = 1:ndepths
    depthi = dataGre2(:,:,i);
    depthi = (depthi - min(depthi(:)))/ std(depthi(:));
    dataGre2(:,:,i) = depthi;
    meansGre(i) = mean(depthi(:));
    depthi = dataRed2(:,:,i);
    depthi = (depthi - min(depthi(:)))/ std(depthi(:));
    dataRed2(:,:,i) = depthi;
    meansRed(i) = mean(depthi(:));
end

figure;
plot(z, meansRed,'r'); hold on
plot(z, meansGre,'g')
xlabel('depth µm'), ylabel('mean luminocity')

% Register the different stacks
% dataGre2 is from optotune stack, dataGre2Regi from the knobby

% Subsample the data
dataGre2sub = dataGre2;
dataRed2sub = dataRed2;
dataGre2Knobsub = dataGre2Regi;
dataRed2Knobsub = dataRed2Regi;
% Before subsampling apply a median filter
for i = 1:size(dataGre2, 3)
    dataGre2sub(:,:,i) = medfilt2(dataGre2sub(:,:,i));
    dataRed2sub(:,:,i) = medfilt2(dataRed2sub(:,:,i));
end
for i = 1:size(dataGre2Knobsub, 3)
    dataGre2Knobsub(:,:,i) = medfilt2(dataGre2Knobsub(:,:,i));
    dataRed2Knobsub(:,:,i) = medfilt2(dataRed2Knobsub(:,:,i));
end
% subsample data
dataGre2sub = dataGre2sub(1:3:end,1:3:end,:);
dataRed2sub = dataRed2sub(1:3:end,1:3:end,:);
dataGre2Knobsub = dataGre2Knobsub(1:3:end,1:3:end,:);
dataRed2Knobsub = dataRed2Knobsub(1:3:end,1:3:end,:);


%% Get the brightness distribution of cells in the z-direction

subplot = @(m,n,p) subtightplot (m, n, p, [0.06 0.03], [0.05 0.05], [0.05 0.05]);

% Coordinates of some neurons (row, col)
% pos = [137, 143;...
%        196, 257;...
%        271, 279;...
%        294, 494;...
%        340, 343;...
%        371, 257;...
%        387, 500;...
%        390, 357;...
%        394, 100;...
%        410, 184;...
%        444, 269;...
%        461, 341;...
%        537, 365;...
%        540, 100;...
%        572, 146;...
%        576, 359;...
%        639, 265;...
%        645, 308;...
%        682, 469;...
%        732, 283];
% pos  = [174 54; 271 152; 227 185; 166 309; 207 438; 325 82; 558 203; 660 80; 650 276;
%         515 33; 544 214; 586 263; 227 135; 601 213; 622 71; 174 276];
% pos = [131 118; 254 125; 257 270; 251 291; 317 188; 372 129; 387 189; 426 324;...
%        459 93;  563 221; 578 169; 561 268];
pos = [353 168; 231 193; 242 227; 250 265; 106 316; 197 388; 328 291; 500 162;...
    603 105; 697 131; 540 232; 736 405];
extra = 1; % Get extra to smooth the signals out
extraR = extra*2+1; % The width of the area that will be extracted per pos

colmap = cmapL([0.8 0 0; 0 0.8 0; 0 0 0.8], extraR^2);
colmap(:,end+1) = 0.6; % the alpha

npos = size(pos,1);
zBright = zeros(npos,length(z));

% Calculate how many subplots need to be made for this many positions
subx = ceil(sqrt(npos));
if (subx-1) * subx >= npos
    suby = subx-1;
else
    suby = subx;
end

figure % figure of zstack
imagesc(x,y,RGBXY); hold on
plot(x(pos(:,1)),y(pos(:,2)),'xy')
text(x(pos(:,1)),y(pos(:,2)),string(1:npos),'Color','y')

% Get the data of the Z stack 
figure; 
for i = 1:npos
    row = pos(i,2);
    col = pos(i,1);
    dataColumn = dataGre(row-extra:row+extra, col-extra:col+extra, :); % column of data
    dataColumn = (dataColumn - min(dataColumn(:)))/ range(dataColumn(:));
    
    zBright(i,:) = squeeze(mean(squeeze(mean(dataColumn))));
    zBright(i,:) = (zBright(i,:)-min(zBright(i,:)))/range(zBright(i,:));
    
    subplot(subx, suby, i) % plotting is important ok
    hold on
    for j = 1:extraR
        for k = 1:extraR
            lineCounter = (j-1)*extra*2+k;
            plot(z, squeeze(dataColumn(j,k,:))', 'color', colmap(lineCounter,:))
        end
    end
    plot(z, zBright(i,:),'linewidth',5)
    xlim([z(1) z(end)])
    xlabel('µm')
end

% Calculate the height of the neurons
[~, maxPos] = max(zBright,[],2);
zBrightCum = sum(zBright,2);

decayPos = zeros(npos,2);
fraction = 2; % 1/fraction van de area wordt gebruikt als cellhoogte threshold
for i = 1:npos
    for j = 1:length(z)
        lim1 = maxPos(i)-j;
        lim2 = maxPos(i)+j-1;
        if lim1<1
            lim1 = 1;
        end
        if lim2>length(z)
            lim2 = length(z);
        end
        areaj = sum(zBright(i,lim1:lim2));
        if areaj > zBrightCum(i)/fraction % als 1/2e van area tussen de limieten
            decayPos(i, :) = [lim1 lim2];
            break
        end
    end
end
cellHeight = [z(decayPos(:,2))-z(decayPos(:,1))]';


% Get brightness based on when the brightness passes a certain value
% threshold = 0.5;
% zBright = zBright - threshold;
% zBright = abs(zBright);
% [~, decayPos] = min(zBright, [], 2);
% 
% cellHeight = abs(double(z(maxPos))-double(z(decayPos)))';

% Finish by plotting the found limits and printing and plotting
for i = 1:npos
    subplot(subx, suby, i)
    line([z(decayPos(i,1)) z(decayPos(i,1))], [0 1], 'color', 'b')
    line([z(decayPos(i,2)) z(decayPos(i,2))], [0 1], 'color', 'b')
    line([z(maxPos(i))   z(maxPos(i))], [0 1], 'color', 'k')
    title(sprintf('%d) Brightness at position [%d, %d, :] h=%dµm',...
        i, row, col, cellHeight(i)))
end
fprintf('minimal cellheight based on 1/%d of signal: %dµm\n',...
    fraction, min(cellHeight))
fprintf('maximal cellheight based on 1/%d of signal: %dµm\n',...
    fraction, max(cellHeight))
fprintf('average cellheight based on 1/%d of signal: %dµm\n',...
    fraction, round(mean(cellHeight)))
figtitle(sprintf('average cellheight: %dµm', round(mean(cellHeight))))

