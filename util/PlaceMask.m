% placing a mask on a maskless recording requires some registration
% 
% Select a recording that has a Mask.
% Then select a recording that requires that Mask.
% Execute the script step by step.
% 
% See also: PlaceMaskAuto, which attempts to do the registration automated
% 
% Leander de Kraker
% 2019-4-1: 17:09-18:09, update 2021-7-13
%


% The recording that has a nice mask
[fnTemplate, pnTemplate] = uigetfile('*SPSIG.mat','select file that has the mask');
fileTemplate = [pnTemplate , fnTemplate];
load(fileTemplate,'Mask','PP','spar');

% The maskless recording that needs the same mask as the template recording
[fnNew, pnNew] = uigetfile('*SPSIG.mat','select file that needs the mask');
fileNew = [pnNew, fnNew];

% Check if the mask that doesn't have a mask yet secretly has a mask
varInfo = who('-file', fileNew);
if any(ismember(varInfo,'Mask'))
    warning('Mask already present in the recording that required a mask! Will overwrite existing Mask')
end

% Load in an image of the maskless recording
if ismember('BImg', varInfo)
    load(fileNew, 'BImg');
    img = BImg;
    fprintf('loaded variable as img: BImg (spectral image)\n')
elseif ismember('BImgMax', varInfo)
    load(fileNew, 'BImgMax');
    img = BImgMax;
    fprintf('loaded variable as img: BImg  (maximum fluorescence)\n')
elseif ismember('BImgAverage', varInfo)
    load(fileNew, 'BImgAverage');
    img = BImgAverage;
    fprintf('loaded variable as img: BImgAverage (average fluorescence)\n')
else
    fprintf('No image found in maskless recording to register Mask to.\n')
    fprintf('Load in an image of the maskless recording manually and call it img\n')
    fprintf('An image is necessary to check correct placement.\n')
end



clearvars pnTemplate fnTemplate pnNew fnNew i Slow BImgMax BImgAverage BImg


%% start the viewimageAlign of Tobias v.d. Bijl
h = viewImageAlign(img, Mask);

%% Apply the selected alignment (DO NOT CLOSE viewImageALign!)

% The transformation applied according to the user interface
transformation.x = h.g.rowOffset;
transformation.y = h.g.columnOffset;
transformation.rotAng = h.g.rotation;

% figure before registration
figure;
imagesc(CreateRGB({img, Mask>0},'bg r'))
hold on
plot(PP.P(1,:),PP.P(2,:),'xk')
for i = 1:PP.Cnt
    plot(PP.Con(i).x, PP.Con(i).y, 'r')
end
title('before registration')


% Apply the transformation to the Mask
MaskNew = RegistrationApply(Mask, transformation);

% Apply the transformation to the contours
PPNew = RegistrationApplyPP(PP, transformation, size(Mask));


% figure after registration
figure;
imagesc(CreateRGB({img, MaskNew>0},'bg r'))
hold on
plot(PPNew.P(1,:),PPNew.P(2,:),'xk')
for i = 1:PPNew.Cnt
    plot(PPNew.Con(i).x, PPNew.Con(i).y, 'r')
end
title('after registration')


%% Check if all ROIs are still present. DELETING ROIs IF NEEDED

presentROIs = unique(MaskNew(:));
presentROIs(presentROIs==0) = [];
missingROIs = 1:PP.Cnt;
missingROIs(presentROIs) = [];
nMissing = length(missingROIs);
if nMissing > 0
    nMissingbefore = length(find(diff(unique(Mask(:)))>1));
    warning('%d ROIs are missing after registration!!!!!!! this WILL result in errors later!', nMissing)
    warning('before registration %d ROIs were missing', nMissingbefore)
    
    % Visualize missing ROIs
    MaskMissingROIs = zeros(size(Mask));
    for i = 1:length(missingROIs)
        MaskMissingROIs(Mask==missingROIs(i)) = 1;
    end
    figure
    imagesc(CreateRGB({MaskMissingROIs, (Mask>1)-MaskMissingROIs}, 'r rgb'))
    
    % Remove missing ROIs from ORIGINAL recording
    answer = questdlg('Delete ROIs from ORIGINAL recording?? Overwriting that file (including signals)',...
                      'Missing ROIs found',...
                      'Delete & save', 'ignore problems', 'Delete & save');
    if strcmp(answer, 'Delete & save')
        fprintf('Removing ROIs that were missing after registration, from original recording\n')
        fprintf('NO PROBLEMS expected in later analysis')
        
        % Remove signals from original template file
        RemoveSigs(fileTemplate, missingROIs)
        % Replace Mask & PP from original template file
        [Mask, PP] = RemoveROIs(Mask, PP, missingROIs);
        save(fileTemplate, 'Mask', 'PP', '-append')
        
        
        % Apply the transformation to the Mask
        MaskNew = RegistrationApply(Mask, transformation);
        % Apply the transformation to the contours
        PPNew = RegistrationApplyPP(PP, transformation, size(Mask));

    else
        warning('Ignoring problems now.') 
        warning('This WILL result in problems in later analysis for registered file which is missing ROIs')
    end
else
    fprintf('No ROIs missing after moving recording. Everything seems okay.\n')
end

%% Save the moved mask and contours and some other variables to the maskless new file

% finalize the variable names
Mask = MaskNew;
PP = PPNew;

% Calculate SpatialCor and PP.Rvar (mean R^2)
filenameTrans = [fileNew(1:end-9), 'DecTrans.dat'];
if ~isfile(filenameTrans)
    filenameTrans = [fileNew(1:end-9), 'Trans.dat'];
end
[sbxt, dim, freq] = transmemap(filenameTrans);
originPos = PP.P([1 2],:);
calcRvar = true;
SpatialCorr = zeros(size(Mask));
fprintf('Calculating mean R2 and SpatialCorr for all ROIs')
for i = 1:PP.Cnt
    [cor, idx, PP.Rvar(i)] = SpatialCorrCalcFun(sbxt, freq, Mask, i, originPos, calcRvar);
    SpatialCorr(idx) = cor(idx); % 1D indexing in 2D matrix
    if mod(i, 25) == 1
        fprintf('Done with %d/%d\n', i, PP.Cnt)
    end
end

save(fileNew, 'PP', 'Mask', 'SpatialCorr', 'spar', '-append')
fprintf('saved.\n')

if ~ismember('BImg', varInfo)
    BImg = img;
    fprintf('saving the reference image as BImg to make use of file easier in the future')
    save(fileNew, 'BImg','-append')
end


