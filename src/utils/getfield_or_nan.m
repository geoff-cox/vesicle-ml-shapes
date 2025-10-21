function v = getfield_or_nan(e, a, b)
    if isstruct(e) && isfield(e,a) && isstruct(e.(a)) && isfield(e.(a), b)
        v = double(e.(a).(b));
    else
        v = NaN;
    end
end