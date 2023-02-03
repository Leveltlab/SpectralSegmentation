function filedate = AskDate(filepath, filename)
% Prompt for getting a date
% If only it was this easy
% 
% Leander de Kraker
% 2022-7-18
% 

questionStr = sprintf('No date detected for file, \nfile:%s\nfolder: %s\n',...
               filename, filepath); 
questionStr = [questionStr, 'Please type its date here as yyyyMMdd:'];

answerStr = inputdlg(questionStr);

try % try to decipher user input as a datetime
    filedate = datetime(answerStr, 'inputformat','yyyyMMdd');
catch METOO % another error: No date for this recording.
    warning('User input not recognized as date')
    filedate = answerStr;
end