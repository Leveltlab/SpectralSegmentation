function EyeData2Avi(varargin)
% Eye mat file to avi file
% eye video is stored in variable called Data and is expected to be a 4D
% uint8. [height x width x color x frames]
% 
% Can be run either as a script or a function
%
% Input:
%   filepath: [string] filepath where the file is and will be located
%   filenameMat: [string] mat file name which holds the eye data
%   deFlicker: [1x1 boolean (true/false)] Try to de-flicker the video
%   invert: [1x1 boolean (true/ false)] Invert the brightness of the video
% 
% 
% Leander de Kraker
% 2021-3-18
% 2021-7-29 Implemented deflickering
% 

if exist('varargin', 'var') && nargin >= 2
    filepath    = varargin{1};
    filenameMat = varargin{2};
else
    [filenameMat, filepath] = uigetfile('*eye.mat', 'select file with eye tracking data');
end

fprintf('Loading eye data, this can take a while...\n')
load([filepath, filenameMat])
filename = filenameMat(1:end-4);
nframes = size(data, 4);
frameDim = [size(data, 1), size(data, 2)];
nColChannel = size(data, 3);
fprintf('done\n')

% data2 = data;

%% Process the file. step 1: increase brightness if maximum is way too low

if max(data(:))<180
    % Increase max of video to nearly value 255
    fprintf('increasing data values/brightness\n')
    data = data * uint8(floor(255/double(max(data(:)))));
end

%% Process the file. step 2: De-flicker the video based on mean brightness
% Some of the data suffers from brightness flickering, correct that

framebrightness = mean(squeeze(mean(squeeze(data))));
if exist('varargin','var') && nargin >= 3
    deFlicker = varargin{3};
else
    % Plot average brightness of frames to give indication of flickering
    figure
    n = 400;
    plot(time(1:n), framebrightness(1:n),'.-')
    xlabel('time (sec)'); ylabel('frame average brightness')
    
    deFlicker = questdlg('De-flicker?', 'Brightness should not change too much..',...
                            'yes', 'no', 'no');
    if strcmp(deFlicker,'yes')
        deFlicker = true;
    else
        deFlicker = false;
    end
end


tic
if deFlicker
    baseLine = movmean(framebrightness, 6);
    baseLine = baseLine + prctile(framebrightness, 99)-prctile(baseLine, 99);
    correction = baseLine ./ framebrightness;
    correction = 1*0.1 + correction*0.9; % Don't overcorrect
    for i = 1:nframes
        data(:,:,:,i) = uint8(single(data(:,:,:,i)) .* correction(i));
        if mod(i, 10000)==1
            fprintf('Deflickering frame %d/%d, elapsed time: %.1f sec\n', i, nframes, toc)
        end
    end
    framebrightnessNew = mean(squeeze(mean(squeeze(data))));
    hold off; plot(time, framebrightness, 'color', [0 0 0 0.9])
    hold on;  plot(time, baseLine, 'color', [0.1 0.4 0.7], 'LineWidth', 2)
    plot(time, framebrightnessNew, 'color', [0.2, 0.7, 0.4], 'LineWidth', 1.3);
    legend({'original brightnesses', 'aimed brightnesses', 'new brightnesses'})
end

% % Activate much tighter subplots
% % [subplot margin top&side],[figure bottomspace,topspace],[leftspace,rightspace]
% subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.01], [0.1 0.04], [0.1 0.1]);
% for i = 1:nframes
%     img1 = data(:,:,:,i);
%     img2 = data2(:,:,:,i);
%     subplot(2,2,1)
%     imagesc(img1); caxis([0 255]); title(sprintf('correction: %.2f', correction(i)))
%     
%     subplot(2,2,2)
%     imagesc(img2); caxis([0 255]); title(sprintf('original frame %d', i))
%     
%     subplot(2,2,3)
%     imagesc(img1 - img2); colorbar; caxis([-15 15])
%     colormap(cmapL('italian roast', 256))
%     pause
% end

%%  Process the file. Step 3: Invert brightness: Pupil needs to be dark
if exist('varargin','var')&& nargin == 4
    invert = varargin{4};
else
    figure
    imagesc(data(:,:,:,3)); colormap(hot)
    title(sprintf('frame %d',i))
    colorbar; caxis([0 255]);
    
    invert = questdlg('Invert brightness?', 'Pupil needs to be black, invert?',...
                            'yes', 'no', 'no');
    if strcmp(invert,'yes')
        invert = true;
    else
        invert = false;
    end
end

% Invert brightness
if invert
    data = (max(data(:))+1)-data;
end

fps = 1 / ((time(end)-time(1)) / length(time));

figure
imagesc(data(:,:,:,3))
colormap(gray); caxis([0 255])
title(sprintf('example frame. video framerate = %.5f fps', fps))

%% Save small mat file without the video if the pupil tracking has been done
if exist('eye', 'var')
    save([filepath, filename, '-noVid.mat'], 'time', 'abstime', 'eye')
    fprintf('saved new mat file without eye video\n')
else
    fprintf('no eye tracking present yet.\n')
end


%% Save the avi file
savename = [filepath, filename, '.mp4'];

V = VideoWriter(savename, 'MPEG-4');
V.Quality = 90; % best/ largest file = 100.
V.FrameRate = fps;
open(V)

tic
fprintf('Writing video...\n')
for i = 1:nframes % Per frame (takes slightly longer)
    writeVideo(V, data(:,:,:,i))
end
close(V)
fprintf('done writing video, in %d sec\n', round(toc))


