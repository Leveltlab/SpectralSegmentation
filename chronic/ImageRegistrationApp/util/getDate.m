function filedates = getDate(filenames)
    % Find filedates that is hidden somewhere in filename
    % Format: yyyyMMdd  = 20220708
    % 
    % Leander de Kraker
    % 2022-7-8
    % 
    % Update from the original code
    % Augustijn Vrolijk
    % 2025-02-13
    arguments (Input)
        filenames (:, 1) string 
    end
    arguments (Output)
        filedates (:, 1) cell
    end
    
    desiredLength = 8; %hardcoded as this is the expected format
    desiredFormat = 'yyyyMMdd';
    nFiles = length(filenames);
    
    filedates = cell(nFiles, 1);
    %NaT(nFiles, 1, 'Format',desiredFormat);
    for i=1:nFiles
        filename = filenames{i};
        
        %find date automatically
        strdate = findDateAuto(filename);
        
        % can fail to find it automatically so check the found date
        isValidDate = checkDate(strdate);
        
        % if it failed ask user 
        if ~ isValidDate
            strdate = AskForDate(filename);
        end
        
        try % try to decipher user input as a datetime
            filedates{i} = datetime(strdate, 'inputformat','yyyyMMdd');
        catch METOO % another error: No date for this recording.
            warning('User input not recognized as date')
        end
    
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
    
    function outputDate = AskForDate(filename)
        questionStr = sprintf('No date detected for file, \nfile:%s\n',...
                              filename); 
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

    function strOut = FindNNumInStr(str, n)
        % Find a number that is exactly n digits long in a string
        % 
        % 
        % Leander de Kraker
        % 2022-7-15
        % 
        
        % Get info about the digits in the string
        pos = regexp(str, '\d');
        a=diff(pos);
        b=find([a inf]>1);
        c=diff([0 b]); % length of the sequences
        good = pos(b(c==n));
        
        if isscalar(good)
            strOut = str(good-n+1:good);
        else
            strOut = [];
        end
    end

end