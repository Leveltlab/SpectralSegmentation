function config = resetConfig()
    % Default Setting config for ImageReg2P app for if none exists
    % Augustijn Vrolijk
    
    config = struct;
    [optimiser, metric]  = imregconfig('monomodal');
    
    config.register = struct;
    config.register.metric = metric;
    config.register.optimiser = optimiser;
    config.register.levels = 4;
    config.register.transform = "affine";
    config.register.crossCTransform = "translation";
    
    % WHICH VARIABLES TO LOAD
    config.data = struct;
    config.data.ExtraColumns = {"IDName = getIDName(Name, Date)";...
                                "Date = getDate(Name)"};
    config.data.FetchDataCols = {"Image = getImg(Name, Path)";...
                                 "Mask = getMask(Name, Path)"};
    
    
    % SETTINGS FOR FILE DETECTION WHEN LOADING
    config.fileType = '_SPSIG.mat'; 
    config.fileRemoveDuplicates = false;
    config.folderAskDepth = true;
    
    % Save the config preferences
    name = "default_config";
    if ispref('imageReg2P', name)
        rmpref('imageReg2P', name)
    end
    addpref('imageReg2P', name, config)
    
    % Set this config to be the active one
    if ispref('imageReg2P', 'config_name')
        rmpref('imageReg2P', 'config_name')
    end
    addpref('imageReg2P', 'config_name', name)
    
    fprintf('\imageReg2P Preference config has been set for your MATLAB via resetConfig()\n')
end