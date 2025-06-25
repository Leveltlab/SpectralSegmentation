function config = createConfig(name)
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

    config.register.levels = 4;
    config.register.transform = "affine"; 
    %options include: "translation", "rigid", "similarity", "affine"
    config.register.crossCTransform = "translation"; 
    %options include: "translation", "rigid", "similarity"
    
    config.data.ExtraColumns = {"IDName = getIDName(Name, Date)", "Date = getDate(Name)"};  
    config.data.FetchDataCols = {"Image = getImg(Name, Path)", "Mask = getMask(Name, Path)"};
   
    if ispref('imageReg2P', name)
        rmpref('imageReg2P', name)
    end
    addpref('imageReg2P', name, config)

    if ispref('imageReg2P', 'config_name')
        rmpref('imageReg2P', 'config_name')
    end
    addpref('imageReg2P', 'config_name', name)
    
    clc
end