% --- change signature to return cache too
function [lbl, cache] = getLabel(S, p, cache)
    S.H0 = p;
    [~, meta] = solveAtParams(S, cache);
    lbl   = meta.label;
    cache = meta.cache;   % propagate
end