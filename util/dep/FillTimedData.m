function FillTimedData(timedDatai, timedFile, pn, fn, nRois, nSlices, alignMethod)
fileInfo = dir([pn fn '.sbx']);
timedDatai(1).fileSize = fileInfo.bytes * 10^-9;
timedDatai(1).nRois      = nRois;
timedDatai(1).nSplits    = nSlices;
timedDatai(1).computerName = getenv('COMPUTERNAME');
timedDatai(1).fileName   = fn;
timedDatai(1).filePath   = pn;
timedDatai(1).date       = datetime;
timedDatai(1).alignMethod= alignMethod;
timedDatai(1).memory = memory;
load(timedFile)
n = length(timedData);
timedData(n+1) = timedDatai;
save(timedFile, 'timedData')
