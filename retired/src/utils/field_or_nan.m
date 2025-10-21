function v = field_or_nan(S, k)
    if isfield(S,k) && ~isempty(S.(k)) && isnumeric(S.(k))
        v = double(S.(k));
    else
        v = NaN;
    end
end