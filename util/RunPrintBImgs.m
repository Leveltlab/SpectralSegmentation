% Run printBImgs on many files
% 
% 



picFolder = 'C:\Users\Leander\Documents\Baan\PICLe\Analysis\BImgs\';

% Data on helero2p
miceChronicPath = {'C:\Users\Leander\Documents\Baan\PICLe\17.20.12\Line\';...
                   'C:\Users\Leander\Documents\Baan\PICLe\17.20.12\Monster\';...
                   'C:\Users\Leander\Documents\Baan\PICLe\17.20.12\Venus\';...
                   'C:\Users\Leander\Documents\Baan\PICLe\14.90\Mickey\';...
                   'C:\Users\Leander\Documents\Baan\PICLe\14.90\Timor\';...
                   'C:\Users\Leander\Documents\Baan\PICLe\14.90\Zebedeus\'};
miceChronicName = {'Line_The12_chronic';...
                   'Monster_The12_Chronic';...
                   'Venus_The12(9)_chronic';...
                   'Mickey_The12_Chronic';...
                   'Timor_The12_Chronic';...
                   'Zebedeus_The12_Chronic'};
% miceTodo = [1 3 4];
miceTodo = [1 2 3 4 5 6];
labelsThe12 = {'ori 1', 'passive 1', 'passive 2',...
               'ori 2', 'active 1', 'active 2',...
               'ori 3', 'delay 1', 'delay 2',...
               'ori 4', '8 ori 1', '8 ori 2'};

miceChronicPath = miceChronicPath(miceTodo);
miceChronicName = miceChronicName(miceTodo);
nMice = length(miceChronicPath);

doScalebarPlot = false;
doPlotForOurEyes = false;
zoom = 1.3;



for m = 6
    load([miceChronicPath{m}, miceChronicName{m}], 'filenames', 'filepaths', 'nfiles')
    shiftAmount = zeros(nfiles, 2); shiftMethod = zeros(nfiles, 2);
    for i = 1:nfiles
        [shiftAmount(i,:), shiftMethod(i,:)] = PrintBImgs([], filenames{i}, picFolder, zoom, doScalebarPlot, doPlotForOurEyes);
    end
end