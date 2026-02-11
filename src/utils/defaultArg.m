function v = defaultArg(s, field, vDefault)
    % DEFAULTARG  Return s.(field) if it exists and 
    % is nonempty; otherwise vDefault.
    if isstruct(s) ...
        && isfield(s, field) ...
        && ~isempty(s.(field))
        v = s.(field);
    else
        v = vDefault;
    end
end