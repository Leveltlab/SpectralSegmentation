% Run printBImgs on many files
% 
% Leander de Kraker
% updated 2024-12-20 to adjust file selection
% 


picFolders = 'D:\2Pdata\Leander\Koen\';

%% Option 1: Select many chronic files
% % Data on helero2p
% miceChronicPath = {'C:\Users\Leander\Documents\Baan\PICLe\17.20.12\Line\';...
%                    'C:\Users\Leander\Documents\Baan\PICLe\17.20.12\Monster\';...
%                    'C:\Users\Leander\Documents\Baan\PICLe\17.20.12\Venus\';...
%                    'C:\Users\Leander\Documents\Baan\PICLe\14.90\Mickey\';...
%                    'C:\Users\Leander\Documents\Baan\PICLe\14.90\Timor\';...
%                    'C:\Users\Leander\Documents\Baan\PICLe\14.90\Zebedeus\'};
% miceChronicName = {'Line_The12_chronic';...
%                    'Monster_The12_Chronic';...
%                    'Venus_The12(9)_chronic';...
%                    'Mickey_The12_Chronic';...
%                    'Timor_The12_Chronic';...
%                    'Zebedeus_The12_Chronic'};

% Koen data
miceChronicPath = {'D:\2Pdata\Koen\Muckli\Kopernik\'};
miceChronicName = {'Kopernik_gray_v3_chronic'};



%% Option 2: Select all .sbx files
nMice = length(miceChronicPath);

doScalebarPlot = false;
doPlotForOurEyes = false;
zoom = 2; % ZOOM FACTOR from logbook, info.config.magnification_list(info.config.magnification)

if ~icell(picFolders)
    picFolders = repmat({picFolders}, [nfiles, 1]);
else
    picFolders = filepaths;
end

for m = 1
    load([miceChronicPath{m}, miceChronicName{m}], 'filenames', 'filepaths', 'nfiles')
    shiftAmount = zeros(nfiles, 2); shiftMethod = zeros(nfiles, 2);
    for i = 1:nfiles
        [shiftAmount(i,:), shiftMethod(i,:)] = PrintBImgs(filepaths{i}, filenames{i}, picFolders, zoom, doScalebarPlot, doPlotForOurEyes);
    end
end