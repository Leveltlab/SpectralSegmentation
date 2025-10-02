function config = resetConfig()
    % Setting config for ImageReg2P app
    % Augustijn Vrolijk

    config = struct;
    [optimiser, metric]  = imregconfig('monomodal');
    
    config.register = struct;
    config.register.metric = metric;
    config.register.optimiser = optimiser;
    config.register.levels = 4;
    config.register.transform = "affine";
    config.register.crossCTransform = "translation";
    
    
    config.data = struct;
    config.data.ExtraColumns = {"IDName = getIDName(Name, Date)",...
                                "Date = getDate(Name)"};
    config.data.FetchDataCols = {"Image = getImg(Name, Path)",...
                                 "Mask = getMask(Name, Path)"};
    
    if ispref('imageReg2P', 'default_config')
        rmpref('imageReg2P', 'default_config')
    end
    addpref('imageReg2P', 'default_config', config)
    
    name = "default_config";
    if ispref('imageReg2P', 'config_name')
        rmpref('imageReg2P', 'config_name')
    end
    addpref('imageReg2P', 'config_name', name)
    
    clc
    enda