function cache = loadCacheIfExists(simDir)
% LOADCACHEIFEXISTS  MAT-only cache loader for the unified layout.
% Reads sim-results/cache.mat if present; otherwise returns a fresh struct.

    if isstruct(simDir) && isfield(simDir,'paths') && isfield(simDir.paths,'root')
        % Backward-compat: someone passed the old "S" struct
        simDir = fullfile(simDir.paths.root, 'sim-results');
    end

    f = fullfile(simDir, 'cache.mat');
    if exist(f,'file') == 2
        S = load(f, 'cache');
        if isfield(S, 'cache'), cache = S.cache; else, cache = struct(); end
    else
        cache = struct();
    end

    % Ensure required fields exist
    if ~isfield(cache,'frontier'), cache.frontier = []; end
    if ~isfield(cache,'failures'), cache.failures = {}; end
end
