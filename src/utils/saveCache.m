function saveCache(simDir, cache)
% SAVECACHE  Atomic save of the unified cache at SimResults/cache.mat.

    if isstruct(simDir) && isfield(simDir,'paths') && isfield(simDir.paths,'root')
        % Backward-compat: old "S" struct
        simDir = fullfile(simDir.paths.root, 'SimResults');
    end

    f   = fullfile(simDir, 'cache.mat');
    tmp = [f '.tmp'];
    save(tmp, 'cache', '-v7');
    movefile(tmp, f, 'f');
end