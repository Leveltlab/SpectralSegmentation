function ChooseConfig(choice)
% Choose which config imgReg2P app will use
% 
% Leander de Kraker
% 2025-11-21
% 
arguments (Input)
    choice = []
end

if ispref('imageReg2P') % Only try if the imageReg2P preference is present
    configs = getpref('imageReg2P');
    if ispref('imageReg2P', 'config_name')
        currentConfig = configs.config_name;
    else
        currentConfig = 'NONE';
    end
    configNames = fieldnames(configs);
    configNames(strcmp(configNames, 'config_name')) = [];
    nconfigs = length(configNames);
    
    if nconfigs>1 % Only bother chosing when multiple configs are present
        if isempty(choice) % if no name was given in the input choose now
            choice = listdlg('Name', 'Which config to use?',...
                             'PromptString', sprintf('(currently: %s)', currentConfig),...
                             'ListSize',[300,100],...
                             'SelectionMode','single', ...
                             'ListString',configNames);
        end
        
        % Set the config_name to make imageReg2P load that one on startup
        if ~isempty(choice)
            if ispref('imageReg2P', 'config_name')
                rmpref('imageReg2P', 'config_name')
            end
            addpref('imageReg2P', 'config_name', configNames{choice})
            fprintf('Has set config to %s\n', configNames{choice})
        else
            fprintf('Cancelled choosing of the config. stopping\n')
        end
    else
        fprintf('Only one config present: %s\n', configNames{1})
    end
else
    fprintf('No config preference is set at all for imageReg2P.\n')
    fprintf('Start imageReg2P app or run createConfig("myNew_config")\n')
end
end