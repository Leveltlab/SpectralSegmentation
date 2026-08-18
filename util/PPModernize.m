function PP = PPModernize(pn, fn, doSave)
% PP = PPModernize(pathname, filename) saves an updated PP (contour & ROI
%   information) to the [pathname filename] SPSIG.mat file
% 
% Script requires the transposed .dat file in order to retrieve essential data  
% 
% Code can be run as function or as script
% As function: 
%   input:
%       pathname [string] e.g. 'E:\myData\mousy\'
%       filename [string] *_SPSIG.mat'
%       doSave [boolean] (default = true) save the new PP to file directly
%       
%   output: PP [struct].
% 
% Modernize the PP variable.
% The PP variable has changed over the years. This script will make sure it
% corresponds to the latest version. It has to calculate the properties
% using the DecTrans data
% 
% Essential correct fields in PP are: 
%       PP.Con.x, PP.Con.y, PP.P([1,2],:), PP.Cnt
% 
% If the contour coordinates are bad/ inconsistent with the Mask, run:
%   PPfromMask.m
% 
% Leander de Kraker
% 2020-4-3
% 2026-8-5: optional save
% 

%% Load file
arguments
    pn = [];
    fn = [];
    doSave = true;
end
if isempty(fn)
    [fn, pn] = uigetfile('*SPSIG.mat');
end
if isempty(fn)
    fprintf('\nno valid file given, quitting PPModernize\n')
    return
end
fprintf('loading %s %s\n', pn, fn)
load([pn fn], 'PP', 'BImg', 'Mask', 'SPic', 'Sax')

% obtain average spectral density for each image
imgStack = log(SPic(:,:,2:end));
Sax(1) = []; %first spectral component is the average power over al components

imgStack = setminlevel(imgStack);
imgStackT = permute(imgStack,[2 1 3]); % transpose the SPic variable so it's same as BImg


%% Check SpecProfile (& peakFreq)

calcSpecProfile = false;
if isfield(PP,'SpecProfile')
    fprintf('SpecProfile present\n')
    if size(PP.SpecProfile, 2) ~= PP.Cnt | length(PP.peakFreq) ~= PP.Cnt
        fprintf('But SpecProfile does not have correct number of ROIs\n')
        calcSpecProfile = true;
    end
else
    calcSpecProfile = true;
    fprintf('going to calculate specprofiles...\n')
end

if calcSpecProfile
    % Retrieve the average spectral profile of the ROI
    [PP.SpecProfile, PP.peakFreq, PP.peakVal] = SpecProfileCalcFun(imgStackT, Mask, 1:PP.Cnt, Sax); 
    fprintf('done\n')
end

%% Calcualte A (Area)

Anew = zeros(1, PP.Cnt);
for i = 1:PP.Cnt
    Anew(i) = sum(Mask(:)==i);
end

if isfield(PP, 'A')
    if ~(all(Anew)==PP.A)
        fprintf('When calculating PP.A, results were different then old PP.A\n')
    end 
else
    fprintf('PP.A calculated\n')
end

PP.A = Anew;

%% Check RVar & SPatialCorr

calcRvar = false;
if isfield(PP, 'Rvar')
    fprintf('Rvar is present\n')
    if length(PP.Rvar) ~= PP.Cnt
        calcRvar = true;
        fprintf('But Rvar has different amount of ROIs then PP.Cnt, so calculating new\n')
    end
    if size(PP.Rvar, 1) > size(PP.Rvar, 2)
        % Transpose PP.Rvar
        PP.Rvar = PP.Rvar';
    end
else
    calcRvar = true;
    fprintf('Going to calculate Rvar and spatialcorr\n')
end

if calcRvar
    % Load transposed sbx data
    pnSbx = pn;
    fnSbx = [fn(1:end-10), '_DecTrans.dat']; % Search for decimated trans file first
    if ~exist([pnSbx fnSbx], 'file') % decimated not found -> Search for regular transposed file
        fnSbx = [fn(1:end-10), '_Trans.dat'];
        if ~exist([pnSbx fnSbx], 'file')   % regular transposed not found -> ask user
            [fnSbx, pnSbx] = uigetfile('*TRANS.dat', 'load raw data (_DecTrans.dat file preferred)');
            if isempty(fnSbx)
                warning('cannot calculate ROI signal correlations without transposed data, quiting')
                return
            end
        else
            fprintf('using transposed data (non decimated) to calculate signal ROI correlations\n')
        end
    else
        fprintf('using decimated transposed data to calculate signal ROI correlations\n')
    end
    [sbxt, dim, freq] = transmemap([pnSbx, fnSbx]);
    
    % Calculate SpatialCorr
    SpatialCorrNew = zeros(size(Mask));
    rvar = zeros(PP.Cnt, 1);
    tic
    for i = 1:PP.Cnt
        [Corri, idx, rvari] = SpatialCorrCalcFun(sbxt, freq, Mask, i, PP.P([1 2], i), calcRvar);
        SpatialCorrNew(idx) = Corri(idx); % 1D indexing in 2D matrix
        rvar(i) = rvari;
        if mod(i, 10)==1
            fprintf('done with ROI %d/%d: elapsed time: %.2f minutes\n', i, PP.Cnt, toc/60)
            pause(0.01)
        end
    end
    toc
    PP.Rvar = rvar;
    
    % If SpatialCorr is not present in SPSIG file, save the SpatialCorr
    % that was just calculated
    clear SpatialCorr
    load([pn fn], 'SpatialCorr')
    if ~exist('SpatialCorr', 'var')
        SpatialCorr = SpatialCorrNew;
        if doSave
            fprintf('Saving new spatialCorr\n')
            save([pn fn], 'SpatialCorr', '-append')
        end
    end
end


%% Calculate Roundedness

% check if roundedness was already present
roundednessPresent = false;
if isfield(PP, 'Roundedness')
    if length(PP.Roundedness) == PP.Cnt
        roundednessPresent = true;
        roundednessOld = PP.Roundedness;
    else
        fprintf('Roundedness was present but did not have correct nr ROIs\n')
    end
end

% Calculate the roundness of all ROIs
PP.Roundedness = zeros(1,PP.Cnt);
for i = 1:PP.Cnt
    PP.Roundedness(i) = perimarea(PP.Con(i).x, PP.Con(i).y);
end

% report changed roundedness
if roundednessPresent
    differentRound = sum(PP.Roundedness~=roundednessOld);
    fprintf('%d/%d values of roundedness changed after recalculation\n', differentRound, PP.Cnt)
else
    fprintf('Calculated the roundednesses\n')
end


%% Check PP.P, should only contain seedpoint coordinates

% Check for error in number of ROIs
if size(PP.P,2) ~= PP.Cnt
    fprintf('catastrophic looking errors in PP.P,')
    fprintf(' incorrect number of coordinates for amount of ROIs according to PP.Cnt\n')
    if max(Mask(:)) ~= PP.Cnt
        fprintf('PP.Cnt has different number than the amount of ROIs present in the mask\n')
    end
end

% Check if the 5row exists, which should have corresponded to PP.Rvar
if size(PP.P,1)>=5 && size(PP.P,2) == PP.Cnt
    differentRvar = sum(PP.P(5,:)~=PP.Rvar');
    fprintf('PP.P 5th row had %d\\%d changes compared to PP.Rvar\n', differentRvar, PP.Cnt)
end
if size(PP.P,1)>2
    fprintf('removing rows from PP.P. Should only contain ROI center coordinates (2 rows)\n')
    PP.P(3:end, :) = [];
end


%% Save the updated PP to the file
if doSave
    save([pn fn], 'PP', '-append')
    fprintf('saved PP\n')
end



