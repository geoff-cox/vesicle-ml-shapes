function saveCache(simDir, cache)
% SAVECACHE  Atomic save of the unified cache at sim-results/cache.mat.

    if isstruct(simDir) && isfield(simDir,'paths') && isfield(simDir.paths,'root')
        % Backward-compat: old "S" struct
        simDir = fullfile(simDir.paths.root, 'sim-results');
    end

    f   = fullfile(simDir, 'cache.mat');
    tmp = [f '.tmp'];
    save(tmp, 'cache', '-v7');
    movefile(tmp, f, 'f');
end