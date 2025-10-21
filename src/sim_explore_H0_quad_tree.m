function sim_explore_H0_quad_tree(sim)
% SIM_EXPLORE_H0_QUAD_TREE  Driver for exploring (H0_1,H0_2) with fixed physics in sim.MP
%
% Uses utils/ for everything:
%   - catalog_load.m, catalog_append.m, catalog_save.m
%   - simpleDataHash.m
%   - loadCacheIfExists.m, saveCache.m
%   - processQuadtree.m
%   - solveAtParams.m
%   - logmsg.m  (optional but recommended)
%
% Conventions:
%   SimResults/
%     ├─ hashed_results/      % immutable blobs <hash>.mat
%     ├─ catalog.mat          % 3-col or legacy wide; handled by utils
%     └─ cache.mat            % small, resumable quadtree state

    % -- paths & logging
    projRoot = fileparts(fileparts(mfilename('fullpath')));
    simDir   = fullfile(projRoot,'SimResults');
    solDir   = fullfile(simDir,'hashed_results');
    if ~exist(simDir,'dir'), mkdir(simDir); end
    if ~exist(solDir,'dir'), mkdir(solDir); end %#ok<*NOPTS>
    opfile   = fullfile(simDir,'OPfile.txt');

    say = @(fmt,varargin) fprintf([fmt '\n'], varargin{:});
    if sim.SP.LogToFile, try logmsg(opfile, '[start]'); catch, end, end

    % -- load state
    T     = catalog_load(simDir);
    cache = loadCacheIfExists(simDir);
    if ~isfield(cache,'frontier'), cache.frontier = []; end

    say('[driver] start | version=%s | maxIters=%g', sim.SP.ModelVersion, sim.SP.MaxIters);

    % -- main loop
    iters = 0;
    while iters < sim.SP.MaxIters
        iters = iters + 1;

        % scheduler now takes physics for "seen?" checks if needed
        [task, cache] = processQuadtree(cache, T, sim.MP);
        if isempty(task)
            say('[driver] no pending tasks; exiting.'); break;
        end

        % Build physics-aware key for hashing (H0 + MP + ABI)
        H0  = [task.params.H0_1, task.params.H0_2];
        key = struct('model_version', string(sim.SP.ModelVersion), ...
                     'H0_1', H0(1), 'H0_2', H0(2), ...
                     'A', sim.MP.A, 'V', sim.MP.V, ...
                     'KA', sim.MP.KA, 'KB', sim.MP.KB, 'KG', sim.MP.KG);
        hash = simpleDataHash(key, 'SHA-256');

        fmat = fullfile(solDir, hash + ".mat");

        % fast path
        if exist(fmat,'file')
            say('[skip] %s (exists)', hash);
            try S = load(fmat,'meta'); catch, S = struct('meta',struct()); end
            entry = struct('params', merge_params(task.params, sim.MP), 'meta', S.meta);
            T = catalog_append(simDir, hash, entry);
            saveCache(simDir, cache);
            continue
        end

        % warm start (from catalog solves or seeds)
        warm = pickWarmStart(task.params, sim, simDir, T);

        % solve with explicit sim context (physics + thresholds/opts)
        try
            [result, meta] = solveAtParams(task.params, sim, warm);
            meta.hash    = hash;
            if ~isfield(meta,'version') || strlength(string(meta.version))==0
                meta.version = string(sim.SP.ModelVersion);
            end
        catch ME
            say('[error] %s :: %s', hash, ME.message);
            if ~isfield(cache,'failures'), cache.failures = {}; end
            cache.failures{end+1} = struct('hash',hash,'msg',ME.message);
            saveCache(simDir, cache);
            continue
        end

        % atomic write
        tmp = [fmat '.tmp'];
        save(tmp, 'result','meta','-v7.3');
        movefile(tmp, fmat, 'f');

        % append catalog with H0 + physics in params
        entry = struct('params', merge_params(task.params, sim.MP), 'meta', meta);
        T = catalog_append(simDir, hash, entry);

        % persist cache
        saveCache(simDir, cache);
    end

    say('[driver] done | iterations=%d | catalogRows=%d', iters, height(catalog_load(simDir)));
end

function P = merge_params(h0params, MP)
    P = h0params;
    fields = fieldnames(MP);
    for i=1:numel(fields), P.(fields{i}) = MP.(fields{i}); end
    % derived scales for auditability
    if ~isfield(P,'aS') || ~isfield(P,'bS')
        [P.aS,P.bS] = computePhaseScales(P.A);
    end
end

