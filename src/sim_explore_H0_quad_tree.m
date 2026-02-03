function sim_explore_H0_quad_tree(sim)

    % -- paths & logging
    projRoot = fileparts(fileparts(mfilename('fullpath')));
    simDir   = fullfile(projRoot,'SimResults');
    solDir   = fullfile(simDir,'hashed_results');
    if ~exist(simDir,'dir'), mkdir(simDir); end
    if ~exist(solDir,'dir'), mkdir(solDir); end %#ok<*NOPTS>
    opfile   = fullfile(simDir,'OPfile.txt');

    say = @(fmt,varargin) fprintf([fmt '\n'], varargin{:});
    if isfield(sim.SP,'LogToFile') && sim.SP.LogToFile
        try logmsg(opfile, '[start]'); catch, end
    end

    % -- load state
    T     = catalog_load(simDir);
    cache = loadCacheIfExists(simDir);
    if ~isfield(cache,'frontier'), cache.frontier = []; end

    % ---------------- failure registry (driver-level backoff) ----------------
    % Canonical per-hash record:
    %   hash, count, lastIter, blockUntil, lastMsg
    cache = normalize_failure_registry(cache);

    % Policy knobs (allow overrides from sim.SP if present)
    MAX_RETRIES   = getfield_default(sim.SP, 'MaxRetries', 3);     % permanent skip after this many failures
    BACKOFF_STEPS = getfield_default(sim.SP, 'BackoffSteps', [0 1 2 4 8 16]); % cooldown steps by failure count
    MIN_BACKOFF   = getfield_default(sim.SP, 'MinBackoff', 1);

    say('[driver] start | version=%s | maxIters=%g', sim.SP.ModelVersion, sim.SP.MaxIters);

    % -- main loop
    iters = 0;
    while iters < sim.SP.MaxIters
        iters = iters + 1;

        % scheduler selects next (H0_1,H0_2) task
        [task, cache] = processQuadtree(cache, T, sim.MP);
        if isempty(task)
            say('[driver] no pending tasks; exiting.');
            break;
        end

        % Build physics-aware key for hashing (H0 + MP + ABI)
        H0  = [task.params.H0_1, task.params.H0_2];
        key = struct('model_version', string(sim.SP.ModelVersion), ...
                     'H0_1', H0(1), 'H0_2', H0(2), ...
                     'A', sim.MP.A, 'V', sim.MP.V, ...
                     'KA', sim.MP.KA, 'KB', sim.MP.KB, 'KG', sim.MP.KG);
        hash = simpleDataHash(key, 'SHA-256');
        fmat = fullfile(solDir, hash + ".mat");

        % --------- 0) Skip if result already exists ----------
        if exist(fmat,'file')
            say('[skip] %s (exists)', hash);
            try S = load(fmat,'meta'); catch, S = struct('meta',struct()); end
            entry = struct('params', merge_params(task.params, sim.MP), 'meta', S.meta);
            T = catalog_append(simDir, hash, entry);
            saveCache(simDir, cache);
            continue
        end

        % --------- 1) Consult failure registry (cooldown / cap) ----------
        [hasRec, recIdx] = failure_lookup(cache.failures, hash);
        if hasRec
            rec = cache.failures(recIdx);

            % permanent skip
            if rec.count >= MAX_RETRIES
                say('[skip] %s (permanent after %d failures: "%s")', ...
                    hash, rec.count, scrub(rec.lastMsg));
                % IMPORTANT: rotate queue so we don't immediately re-hit the same unsolved corner
                cache = rotate_quadtree_queue(cache);
                saveCache(simDir, cache);
                continue
            end

            % cooldown skip
            if iters < rec.blockUntil
                say('[skip] %s (cooldown: wait %d iters, last="%s")', ...
                    hash, rec.blockUntil - iters, scrub(rec.lastMsg));
                cache = rotate_quadtree_queue(cache);
                saveCache(simDir, cache);
                continue
            end
        end

        % --------- 2) Warm start (from catalog solves or seeds) ----------
        warm = pickWarmStart(task.params, sim, simDir, T);

        % --------- 3) Solve ----------
        try
            [result, meta] = solveAtParams(task.params, sim, warm);
            meta.hash = hash;
            if ~isfield(meta,'version') || strlength(string(meta.version))==0
                meta.version = string(sim.SP.ModelVersion);
            end

            % On success, clear any failure record for this hash
            if hasRec
                cache.failures(recIdx) = [];   % remove entry
            end

        catch ME
            say('[error] %s :: %s', hash, ME.message);

            % update failure registry with backoff
            if hasRec
                rec.count      = rec.count + 1;
                rec.lastIter   = iters;
                cool           = BACKOFF_STEPS(min(rec.count, numel(BACKOFF_STEPS)));
                rec.blockUntil = iters + max(cool, MIN_BACKOFF);
                rec.lastMsg    = ME.message;
                cache.failures(recIdx) = rec;
            else
                rec = struct('hash',hash, 'count',1, 'lastIter',iters, ...
                             'blockUntil', iters + max(BACKOFF_STEPS(1), MIN_BACKOFF), ...
                             'lastMsg', ME.message);
                cache.failures(end+1) = rec; %#ok<AGROW>
            end

            cache = rotate_quadtree_queue(cache);
            saveCache(simDir, cache);
            continue
        end

        % --------- 4) Persist success ----------
        tmp = fmat + ".tmp";
        save(tmp, 'result','meta','-v7.3');
        movefile(tmp, fmat, 'f');

        % catalog append with full params
        entry = struct('params', merge_params(task.params, sim.MP), 'meta', meta);
        T = catalog_append(simDir, hash, entry);

        saveCache(simDir, cache);
    end

    say('[driver] done | iterations=%d | catalogRows=%d', iters, height(catalog_load(simDir)));
end

% ---------------- local helpers ----------------

function P = merge_params(h0params, MP)
    P = h0params;
    fields = fieldnames(MP);
    for i=1:numel(fields), P.(fields{i}) = MP.(fields{i}); end
    if ~isfield(P,'aS') || ~isfield(P,'bS')
        [P.aS,P.bS] = computePhaseScales(P.A);
    end
end

function cache = rotate_quadtree_queue(cache)
    % Prevent immediate reselection of the same unsolved corner after skip/error.
    if isfield(cache,'QT') && isfield(cache.QT,'queue') && numel(cache.QT.queue) > 1
        cache.QT.queue = [cache.QT.queue(2:end), cache.QT.queue(1)];
    end
end

function cache = normalize_failure_registry(cache)
    % Canonical struct array: hash,count,lastIter,blockUntil,lastMsg
    if ~isfield(cache,'failures') || isempty(cache.failures)
        cache.failures = struct('hash',{},'count',{},'lastIter',{},'blockUntil',{},'lastMsg',{});
        return;
    end

    % If prior runs used a cell array (your current active code appends to a cell) :contentReference[oaicite:6]{index=6}
    if iscell(cache.failures)
        F = struct('hash',{},'count',{},'lastIter',{},'blockUntil',{},'lastMsg',{});
        for k=1:numel(cache.failures)
            x = cache.failures{k};
            if isstruct(x) && isfield(x,'hash')
                msg = "";
                if isfield(x,'msg'), msg = string(x.msg); end
                F(end+1) = struct('hash',string(x.hash), 'count',1, 'lastIter',0, 'blockUntil',0, 'lastMsg',msg); %#ok<AGROW>
            end
        end
        cache.failures = F;
        return;
    end

    % If already struct array, enforce fields
    if isstruct(cache.failures)
        need = {'hash','count','lastIter','blockUntil','lastMsg'};
        for i=1:numel(need)
            if ~isfield(cache.failures, need{i})
                [cache.failures.(need{i})] = deal([]);
            end
        end
        % normalize types
        for k=1:numel(cache.failures)
            cache.failures(k).hash = string(cache.failures(k).hash);
            if isempty(cache.failures(k).count),      cache.failures(k).count      = 1;  end
            if isempty(cache.failures(k).lastIter),   cache.failures(k).lastIter   = 0;  end
            if isempty(cache.failures(k).blockUntil), cache.failures(k).blockUntil = 0;  end
            if isempty(cache.failures(k).lastMsg),    cache.failures(k).lastMsg    = ""; end
        end
    end
end

function [found, idx] = failure_lookup(F, hash)
    if isempty(F), found = false; idx = 0; return; end
    idx = find(strcmp({F.hash}, string(hash)), 1, 'first');
    found = ~isempty(idx);
    if ~found, idx = 0; end
end

function s = scrub(msg)
    s = string(msg);
    if strlength(s)==0, return; end
    s = regexprep(s, '\s+', ' ');
    if strlength(s) > 120, s = extractBefore(s, 121) + "..."; end
end

function v = getfield_default(S, name, def)
    v = def;
    if isstruct(S) && isfield(S, name)
        v = S.(name);
    end
end

function cache = loadCacheIfExists(simDir)
    % LOADCACHEIFEXISTS  MAT-only cache loader for the unified layout.
    % Reads SimResults/cache.mat if present; otherwise returns a fresh struct.

    if isstruct(simDir) && isfield(simDir,'paths') && isfield(simDir.paths,'root')
        % Backward-compat: someone passed the old "S" struct
        simDir = fullfile(simDir.paths.root, 'SimResults');
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

function warm = pickWarmStart(params, sim, simDir, T)

    warm = struct();
    if nargin < 4 || isempty(T), T = catalog_load(simDir); end
    MP = sim.MP; Htgt = [params.H0_1 params.H0_2];

    tol = 1e-12;

    % 1) nearest solved neighbor with matching physics
    if ~isempty(T) && any(T.Properties.VariableNames=="entry")
        H1 = cellfun(@(e) fget(e,'params','H0_1'), T.entry);
        H2 = cellfun(@(e) fget(e,'params','H0_2'), T.entry);
        Av = cellfun(@(e) fget(e,'params','A'),   T.entry);
        Vv = cellfun(@(e) fget(e,'params','V'),   T.entry);
        KAv= cellfun(@(e) fget(e,'params','KA'),  T.entry);
        KBv= cellfun(@(e) fget(e,'params','KB'),  T.entry);
        KGv= cellfun(@(e) fget(e,'params','KG'),  T.entry);

        isSolve = cellfun(@(e) ~isstruct(e.meta) || ~isfield(e.meta,'type') ...
                                   || ~strcmpi(string(e.meta.type),'seed'), T.entry);

        mask = isSolve ...
            & abs(Av-MP.A)<tol & abs(Vv-MP.V)<tol ...
            & abs(KAv-MP.KA)<tol & abs(KBv-MP.KB)<tol & abs(KGv-MP.KG)<tol ...
            & isfinite(H1) & isfinite(H2);

        if any(mask)
            d2 = (H1(mask)-Htgt(1)).^2 + (H2(mask)-Htgt(2)).^2;
            [~,krel] = min(d2); idx = find(mask); ix = idx(krel);
            if isfield(T.entry{ix}.meta,'hash')
                f = fullfile(simDir,'hashed_results', T.entry{ix}.meta.hash + ".mat");
                if exist(f,'file')==2
                    S = load(f,'result');
                    if isfield(S,'result') && isfield(S.result,'sol')
                        warm.result = S.result;
                        warm.sol = S.result.sol; % compatibility
                        warm.fromParams = struct('H0_1',H1(ix),'H0_2',H2(ix));
                        return
                    end
                end
            end
        end
    end

    % 2) seed from initial-shapes
    if ~isempty(T) && any(T.Properties.VariableNames=="entry")
        isSeed = cellfun(@(e) isstruct(e.meta) && isfield(e.meta,'type') ...
                                 && strcmpi(string(e.meta.type),'seed'), T.entry);
        Av = cellfun(@(e) fget(e,'params','A'),   T.entry);
        Vv = cellfun(@(e) fget(e,'params','V'),   T.entry);
        KAv= cellfun(@(e) fget(e,'params','KA'),  T.entry);
        KBv= cellfun(@(e) fget(e,'params','KB'),  T.entry);
        KGv= cellfun(@(e) fget(e,'params','KG'),  T.entry);

        mask = isSeed & abs(Av-MP.A)<tol & abs(Vv-MP.V)<tol & abs(KAv-MP.KA)<tol & abs(KBv-MP.KB)<tol & abs(KGv-MP.KG)<tol;
        if any(mask)
            ix = find(mask,1,'first');
            src = ""; if isfield(T.entry{ix}.meta,'source'), src = string(T.entry{ix}.meta.source); end
            ishapesDir = fullfile(fileparts(simDir),'initial-shapes');
            f = src; if ~isempty(src) && exist(f,'file')~=2, f = fullfile(ishapesDir, src); end
            if exist(f,'file')==2
                tmp = load(f);
                if isfield(tmp,'Version') && numel(tmp.Version)>=1 && isfield(tmp.Version(1),'Solution')
                    warm.sol = tmp.Version(1).Solution;
                    warm.result = struct('sol', warm.sol); % critical fix
                    warm.fromParams = struct('H0_1',0,'H0_2',0);
                    return
                end
            end
        end
    end

    function v = fget(e, a, b)
        if isstruct(e) && isfield(e,a) && isstruct(e.(a)) && isfield(e.(a),b)
            v = double(e.(a).(b));
        else
            v = NaN;
        end
    end

end

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