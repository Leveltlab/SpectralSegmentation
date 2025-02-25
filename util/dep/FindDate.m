function filedate = FindDate(filepath, filename)
% Find a filedate that is hidden somewhere in filepath or filename
% Format: yyyyMMdd  = 20220708
% 
% Leander de Kraker
% 2022-7-8
% 
% Update from the original code
% Augustijn Vrolijk
% 2025-02-13

desiredLength = 8; %hardcoded as this is the expected format
desiredFormat = 'yyyyMMdd';
filestr = [filepath, filename];


%find date automatically
strdate = findDateAuto(filestr);

% can fail to find it automatically so check the found date
isValidDate = checkDate(strdate);

% if it failed ask user 
if ~ isValidDate
    strdate = AskForDate(filepath, filename);
end

try % try to decipher user input as a datetime
    filedate = datetime(strdate, 'inputformat','yyyyMMdd');
catch METOO % another error: No date for this recording.
    warning('User input not recognized as date')
    filedate = strdate;
end


    function isValid = checkDate(inputDate)
        try
            tempVal = datetime(inputDate, 'inputformat',desiredFormat);
            if isnat(tempVal) 
                isValid = false;
            else    
                isValid = true;
            end
        catch METOO % Date extraction failed.  Ask user for input
            isValid = false;
        end
    end
    
    function outputDate = AskForDate(filepath, filename)
        questionStr = sprintf('No date detected for file, \nfile:%s\nfolder: %s\n',...
                              filename, filepath); 
        questionStr = [questionStr, 'Please type its date here as yyyyMMdd:'];
        outputDate = inputdlg(questionStr);
    end

    function outputDate = findDateAuto(filestr)
        %findNNumInStr returns [], an empty array if it can't find a
        %sequence of consecutive numbers of desiredLength

        outputDate = FindNNumInStr(filestr, desiredLength);
        if isempty(outputDate) % Hope the date was not found because it has dashes
            filestr = strrep(filestr,'-','');
            outputDate = FindNNumInStr(filestr, desiredLength);
            if isempty(outputDate) % Last hope, maybe different dashes?
                filestr = strrep(filestr,'_','');
                outputDate = FindNNumInStr(filestr, desiredLength);
            end
        end
    end
end