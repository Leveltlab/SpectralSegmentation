function SummaryGetRois(rejlog, spar)
% Report a summary of the automatic ROI creation from getSpectrois
% 
% Chris v.d. Togt 
% 2020
% 

total = size(rejlog,1);

% valid rois

vidx = (rejlog(:,1)== 2);
valnm = sum(vidx);
if valnm > 10
    valid = rejlog(vidx, :);

    %valid stats:
    vMd = round(valnm/2);
    v5 = ceil(valnm*0.05);
    v95 = floor(valnm*0.95);

    vSort = sort(valid(:,2)); % Areas
    vAm = vSort(vMd);
    Arange = [vSort(v5) vSort(v95)];

    vSort = sort(valid(:,3)); % Roundedness
    vRm = vSort(vMd);
    Rrange = [vSort(v5) vSort(v95)];

    vSort = sort(valid(:,4)); % Rvariance
    vRvarm = vSort(vMd);
    Rvarange = [vSort(v5) vSort(v95)];
end

% invalid maxima
rejected = sum(rejlog(:,1) ~= 2);  % All
notsignificant = sum(rejlog(:,1)== -1); % Not significant height

remidx = ( rejlog(:,1)== 0 | rejlog(:,1)== 1 ); % rejected
remaining = sum(remidx);
notvalid = rejlog(remidx,:);

tosmall = sum(notvalid(:,2) <= spar.areasz(1));
tolarge = sum(notvalid(:,2) > spar.areasz(2));

belowround = sum(notvalid(:,3) < spar.roundedness);
% bad = sum(notvalid(:,4) > 0.3);

someidx = isnan(rejlog(:,1));
Inother = sum(someidx);

disp("Of " + total + " peaks(maxima), " + valnm + " were valid and " + rejected + " were rejected.")

if valnm > 10
    disp("valid stats:")
    disp("Area(pixels) median and 95% range: " + vAm + " : " + Arange(1) + " - " + Arange(2))
    disp("Roundedness median and 95% range: " + vRm + " : " + Rrange(1) + " - " + Rrange(2))
    disp("Rvar median and 95% range: " + vRvarm + " : " + Rvarange(1) + " - " + Rvarange(2))
end

disp("Of " + rejected + " rejected peaks:")
disp(notsignificant + " were not significant")
disp(Inother + " were within previously made rois" )
disp("Of the remaining " + remaining + " peaks the contours were not valid because;")
disp(tosmall + " were too small, " + tolarge + " too large, " + belowround + " not round enough,")
% disp("and " + bad + " had more than 40% below threshold pixel values") 
disp(" ")