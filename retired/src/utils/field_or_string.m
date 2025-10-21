function v = field_or_string(S, k)
    if isfield(S,k) && ~isempty(S.(k))
        v = string(S.(k));
    else
        v = "";
    end
end