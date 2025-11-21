function config = createConfig(name)
% This function creates a config file which will get used by the
% ImgRegMain app.
% 
% 
% 
    arguments (Input)
        name string = "NOTVALIDNAME"
    end
    if strcmp(name, "NOTVALIDNAME")
        error("no name given for the configuration, either add one manually or call createConfig with a name")
    end
    
    config = struct;
    config.register = struct;
    config.data = struct;
    
    %ALL IMAGE REGISTRATION INFO CAN BE FOUND AT:
    %https://nl.mathworks.com/help/images/image-registration.html
    
    [optimiser, metric]  = imregconfig('monomodal');
    %for default configurations check out: https://www.mathworks.com/help/images/ref/imregconfig.html
    
    config.register.metric = metric;
    %options include:
    %https://nl.mathworks.com/help/images/ref/registration.metric.mattesmutualinformation.html
    %https://nl.mathworks.com/help/images/ref/registration.metric.meansquares.html
    
    config.register.optimiser = optimiser;
    %options include:
    %https://nl.mathworks.com/help/images/ref/registration.optimizer.oneplusoneevolutionary.html
    %https://nl.mathworks.com/help/images/ref/registration.optimizer.regularstepgradientdescent.html
    
    config.register.levels = 5;
    config.register.transform = "affine"; 
    %options include: "translation", "rigid", "similarity", "affine"
    config.register.crossCTransform = "translation"; 
    %options include: "translation", "rigid", "similarity"
    
    
    % WHICH VARIABLES TO LOAD
    % data which result in images which will be registerable
    config.data.FetchDataCols = {"Image = getImg(Name, Path)";...
                                 "Mask  = getMask(Name, Path)"};
    % additional file information
    config.data.ExtraColumns = {"IDName = getIDName(Name, Date)";...
                                "Date   = getDate(Name)"};  
    
    
    % SETTINGS FOR FILE DETECTION WHEN LOADING
    config.fileType = '_SPSIG.mat'; % ONLY files with this extension 
                                    % get loaded when using 'Load Folder',
                                    % and this file extension is proposed
                                    % when using 'Load File'
    config.fileRemoveDuplicates = true; % detect duplicates and delete 
                                  % from folder selection? See
                                  % the function CollectFiles to change
                                  % how duplicates are detected
    config.folderAskDepth = false; % Ask how many folders deep to search 
                            % (assumes 0 (selected folder only) if false)
    
    
    % APPLY CONFIG
    if ispref('imageReg2P', name)
        rmpref('imageReg2P', name)
    end
    addpref('imageReg2P', name, config)
    
    % Set this config to be the active one
    if ispref('imageReg2P', 'config_name')
        rmpref('imageReg2P', 'config_name')
    end
    addpref('imageReg2P', 'config_name', name)
    fprintf('Active config set to %s\n', name)
end