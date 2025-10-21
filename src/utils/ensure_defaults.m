function cache = ensure_defaults(cache)
    if ~isfield(cache,'config'), cache.config = struct(); end
    if ~isfield(cache.config,'bounds')
        % Default box if not provided: [-1, +1] x [-1, +1]
        cache.config.bounds = [-1 1; -1 1];
    end
    cache.config.maxDepth = defaultArg(cache.config,'maxDepth', 6);
    cache.config.maxCells = defaultArg(cache.config,'maxCells', 2000);
    cache.config.eTol     = defaultArg(cache.config,'eTol',     5e-3);
    cache.config.pTol     = defaultArg(cache.config,'pTol',     5e-3);
    cache.config.shapeTau = defaultArg(cache.config,'shapeTau', 0.08);
end