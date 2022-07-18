function filedate = FindDate(filepath, filename)
% Find a filedate that is hidden somewhere in filepath or filename
% Format: yyyyMMdd  = 20220708
% 
% Leander de Kraker
% 2022-7-8
% 

filestr = [filepath, filename];

% Find exactly 8 consecutive digids, that's probably the date
desiredLength = 8;

strdate = FindNNumInStr(filestr, desiredLength);
if isempty(strdate) % Hope the date was not found because it has dashes
    filestr = strrep(filestr,'-','');
    strdate = FindNNumInStr(filestr, desiredLength);
    if isempty(strdate) % Last hope, maybe different dashes?
        filestr = strrep(filestr,'_','');
        strdate = FindNNumInStr(filestr, desiredLength);
    end
end

try
    filedate = datetime(strdate, 'inputformat','yyyyMMdd');
catch ME % Date extraction failed.  Ask user for input
    if strcmp(ME.identifier, 'MATLAB:datetime:ParseErr') %
        filedate = AskDate(filepath, filename);
    else % extraction failed because of unkown reason
        throw(ME)
    end
end

if isnat(filedate) 
    filedate = AskDate(filepath, filename);
end