function options = mergeStructs(defaults, userOptions)
    % Get field names of both structs
    userFields = fieldnames(userOptions);
    
    % Copy default values
    options = defaults;
    
    % Overwrite defaults with user-specified values
    for i = 1:numel(userFields)
        options.(userFields{i}) = userOptions.(userFields{i});
    end
end