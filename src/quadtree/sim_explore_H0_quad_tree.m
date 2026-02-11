function sim_explore_H0_quad_tree(sim)

    % --- project root and output dirs (allow override for tests) ---
    here     = fileparts(mfilename('fullpath'));   % src/quadtree
    srcRoot  = fileparts(here);                    % src
    projRoot = fileparts(srcRoot);                 % repo root

    simDir = getfield_default(sim.SP, 'SimDir', fullfile(projRoot,'sim-results'));
    solDir = fullfile(simDir,'hashed_results');
    if ~exist(simDir,'dir'), mkdir(simDir); end
    if ~exist(solDir,'dir'), mkdir(solDir); end %#ok<*NOPTS>

    opfile = fullfile(simDir,'OPfile.txt');
    say = @(fmt,varargin) fprintf([fmt '\n'], varargin{:});

    if getfield_default(sim.SP,'LogToFile',true)
        try logmsg(opfile, '[start]'); catch, end %#ok<CTCH>
    end

    % --- load catalog and cache ---
    T     = catalog_load(simDir);
    cache = loadCacheIfExists(simDir);
    if ~isfield(cache,'frontier'), cache.frontier = []; end

    % --- apply config overrides from sim.SP (bounds/depth/etc) ---
    cache = apply_cache_overrides(cache, sim.SP);

    % --- failure registry normalization ---
    cache = normalize_failure_registry(cache);

    MAX_RETRIES   = getfield_default(sim.SP, 'MaxRetries', 3);
    BACKOFF_STEPS = getfield_default(sim.SP, 'BackoffSteps', [0 1 2 4 8 16]);
    MIN_BACKOFF   = getfield_default(sim.SP, 'MinBackoff', 1);

    SAVE_EVERY    = getfield_default(sim.SP, 'CacheSaveEvery', 50);  % reduce I/O
    USE_WARMSTART = getfield_default(sim.SP, 'UseWarmStart', true);

    % ---- hysteresis: branch tag for multi-solution support ----
    BRANCH_TAG    = string(getfield_default(sim.SP, 'BranchTag', ""));

    solverFcn = getfield_default(sim.SP, 'SolverFcn', @solveAtParams);

    say('[driver] start | version=%s | maxIters=%g', string(sim.SP.ModelVersion), sim.SP.MaxIters);

    % Store physics context for hash-aware failure matching in processQuadtree
    cache.config.modelVersion = string(sim.SP.ModelVersion);
    cache.config.MP = sim.MP;
    cache.config.branchTag = BRANCH_TAG;

    iters = 0;
    while iters < sim.SP.MaxIters
        iters = iters + 1;
        cache.iter = iters;  % <-- critical for cooldown-aware scheduling

        [task, cache] = processQuadtree(cache, T, sim.MP);
        if isempty(task)
            say('[driver] no pending tasks; exiting.');
            break;
        end

        H0  = [task.params.H0_1, task.params.H0_2];

        % --- physics-aware hash key (model_version must be included) ---
        key = struct( ...
            'model_version', string(sim.SP.ModelVersion), ...
            'H0_1', H0(1), 'H0_2', H0(2), ...
            'A', sim.MP.A, 'V', sim.MP.V, ...
            'KA', sim.MP.KA, 'KB', sim.MP.KB, 'KG', sim.MP.KG);

        % hysteresis: include branch_tag when non-empty
        if strlength(BRANCH_TAG) > 0
            key.branch_tag = BRANCH_TAG;
        end

        hash = simpleDataHash(key, 'SHA-256');  % now safe with patched simpleDataHash
        fmat = fullfile(solDir, hash + ".mat");

        % --- skip if already computed ---
        if exist(fmat,'file')
            say('[skip] %s (exists)', hash);
            try S = load(fmat,'meta'); catch, S = struct('meta',struct()); end %#ok<CTCH>
            entry = struct('params', merge_params(task.params, sim.MP), 'meta', S.meta);
            T = catalog_append(T, hash, entry);

            if mod(iters, SAVE_EVERY)==0
                catalog_save(simDir, T);
                saveCache(simDir, cache);
            end
            continue
        end

        % --- consult failure registry ---
        [hasRec, recIdx] = failure_lookup(cache.failures, hash);
        if hasRec
            rec = cache.failures(recIdx);

            if rec.count >= MAX_RETRIES
                say('[skip] %s (permanent after %d failures: "%s")', hash, rec.count, scrub(rec.lastMsg));
                % Mark permanently blocked so processQuadtree's isBlocked()
                % will skip this corner instead of returning it again.
                cache.failures(recIdx).blockUntil = Inf;
                cache = rotate_quadtree_queue(cache);
                if mod(iters, SAVE_EVERY)==0
                    catalog_save(simDir, T);
                    saveCache(simDir, cache);
                end
                continue
            end

            if iters < rec.blockUntil
                say('[skip] %s (cooldown: wait %d iters, last="%s")', hash, rec.blockUntil - iters, scrub(rec.lastMsg));
                cache = rotate_quadtree_queue(cache);
                if mod(iters, SAVE_EVERY)==0
                    catalog_save(simDir, T);
                    saveCache(simDir, cache);
                end
                continue
            end
        end

        % --- warm start (optional) ---
        warm = struct();
        if USE_WARMSTART
            warm = pickWarmStart(task.params, sim, simDir, T);
        end

        % --- solve ---
        try
            [result, meta] = solverFcn(task.params, sim, warm);
            meta.hash = hash;
            if ~isfield(meta,'version') || strlength(string(meta.version))==0
                meta.version = string(sim.SP.ModelVersion);
            end
            if strlength(BRANCH_TAG) > 0
                meta.branch_tag = BRANCH_TAG;
            end

            if hasRec
                cache.failures(recIdx) = []; % clear failures on success
            end

        catch ME
            say('[error] %s :: %s', hash, ME.message);

            % update failure registry with backoff (store H0 for scheduler)
            if hasRec
                rec.count      = rec.count + 1;
                rec.lastIter   = iters;
                cool           = BACKOFF_STEPS(min(rec.count, numel(BACKOFF_STEPS)));
                rec.blockUntil = iters + max(cool, MIN_BACKOFF);
                rec.lastMsg    = ME.message;
                rec.H0_1       = H0(1);
                rec.H0_2       = H0(2);
                cache.failures(recIdx) = rec;
            else
                rec = struct('hash',hash, 'count',1, 'lastIter',iters, ...
                             'blockUntil', iters + max(BACKOFF_STEPS(1), MIN_BACKOFF), ...
                             'lastMsg', ME.message, 'H0_1', H0(1), 'H0_2', H0(2));
                cache.failures(end+1) = rec; %#ok<AGROW>
            end

            cache = rotate_quadtree_queue(cache);
            if mod(iters, SAVE_EVERY)==0
                catalog_save(simDir, T);
                saveCache(simDir, cache);
            end
            continue
        end

        % --- persist success ---
        tmp = fmat + ".tmp";
        save(tmp, 'result','meta','-v7.3');
        movefile(tmp, fmat, 'f');

        entry = struct('params', merge_params(task.params, sim.MP), 'meta', meta);
        T = catalog_append(T, hash, entry);

        if mod(iters, SAVE_EVERY)==0
            catalog_save(simDir, T);
            saveCache(simDir, cache);
        end
    end

    catalog_save(simDir, T);
    saveCache(simDir, cache);
    say('[driver] done | iterations=%d | catalogRows=%d', iters, height(T));
end

% ---------------- helpers (unchanged where possible) ----------------
function P = merge_params(h0params, MP)
    P = h0params;
    fields = fieldnames(MP);
    for i=1:numel(fields), P.(fields{i}) = MP.(fields{i}); end
    if ~isfield(P,'aS') || ~isfield(P,'bS')
        [P.aS,P.bS] = computePhaseScales(P.A);
    end
end

function cache = rotate_quadtree_queue(cache)
    if isfield(cache,'QT') && isfield(cache.QT,'queue') && numel(cache.QT.queue) > 1
        cache.QT.queue = [cache.QT.queue(2:end), cache.QT.queue(1)];
    end
end

function cache = apply_cache_overrides(cache, SP)
    if ~isfield(cache,'config'), cache.config = struct(); end
    if isfield(SP,'H0Bounds') && isnumeric(SP.H0Bounds) && all(size(SP.H0Bounds)==[2 2])
        cache.config.bounds = SP.H0Bounds;
    end
    if isfield(SP,'QTmaxDepth'), cache.config.maxDepth = SP.QTmaxDepth; end
    if isfield(SP,'QTmaxCells'), cache.config.maxCells = SP.QTmaxCells; end
    if isfield(SP,'eTol'), cache.config.eTol = SP.eTol; end
    if isfield(SP,'pTol'), cache.config.pTol = SP.pTol; end
    if isfield(SP,'shapeTau'), cache.config.shapeTau = SP.shapeTau; end
end

function cache = normalize_failure_registry(cache)
    if ~isfield(cache,'failures') || isempty(cache.failures)
        cache.failures = struct('hash',{},'count',{},'lastIter',{},'blockUntil',{},'lastMsg',{},'H0_1',{},'H0_2',{});
        return;
    end

    if iscell(cache.failures)
        F = struct('hash',{},'count',{},'lastIter',{},'blockUntil',{},'lastMsg',{},'H0_1',{},'H0_2',{});
        for k=1:numel(cache.failures)
            x = cache.failures{k};
            if isstruct(x) && isfield(x,'hash')
                msg = ""; if isfield(x,'msg'), msg = string(x.msg); end
                F(end+1) = struct('hash',string(x.hash), 'count',1, 'lastIter',0, 'blockUntil',0, ...
                                  'lastMsg',msg,'H0_1',NaN,'H0_2',NaN); %#ok<AGROW>
            end
        end
        cache.failures = F;
        return;
    end

    need = {'hash','count','lastIter','blockUntil','lastMsg','H0_1','H0_2'};
    for i=1:numel(need)
        if ~isfield(cache.failures, need{i})
            [cache.failures.(need{i})] = deal([]);
        end
    end
    for k=1:numel(cache.failures)
        cache.failures(k).hash = string(cache.failures(k).hash);
        if isempty(cache.failures(k).count),      cache.failures(k).count      = 1;  end
        if isempty(cache.failures(k).lastIter),   cache.failures(k).lastIter   = 0;  end
        if isempty(cache.failures(k).blockUntil), cache.failures(k).blockUntil = 0;  end
        if isempty(cache.failures(k).lastMsg),    cache.failures(k).lastMsg    = ""; end
        if isempty(cache.failures(k).H0_1),       cache.failures(k).H0_1       = NaN; end
        if isempty(cache.failures(k).H0_2),       cache.failures(k).H0_2       = NaN; end
    end
end

function [found, idx] = failure_lookup(F, hash)
    found = false; idx = 0;
    if isempty(F), return; end
    hashes = string({F.hash});
    idx = find(hashes == string(hash), 1, 'first');
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
    f = fullfile(simDir, 'cache.mat');
    if exist(f,'file') == 2
        S = load(f, 'cache');
        if isfield(S, 'cache'), cache = S.cache; else, cache = struct(); end
    else
        cache = struct();
    end
    if ~isfield(cache,'frontier'), cache.frontier = []; end
    if ~isfield(cache,'failures'), cache.failures = []; end
end

function saveCache(simDir, cache)
    f   = fullfile(simDir, 'cache.mat');
    tmp = f + ".tmp";
    save(tmp, 'cache', '-v7');
    movefile(tmp, f, 'f');
end

function warm = pickWarmStart(params, sim, simDir, T)

    warm = struct();
    if nargin < 4 || isempty(T), T = catalog_load(simDir); end
    MP = sim.MP; Htgt = [params.H0_1 params.H0_2];

    branchTag = string(getfield_default(sim.SP, 'BranchTag', ""));

    tol = 1e-12;

    % 1) nearest solved neighbor with matching physics (and branch)
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

        % branch-aware filtering: prefer same branch tag
        if strlength(branchTag) > 0 && any(mask)
            btags = cellfun(@(e) string(fget_str(e,'meta','branch_tag')), T.entry);
            branchMask = mask & (btags == branchTag);
            if any(branchMask)
                mask = branchMask;
            end
            % if no same-branch neighbor found, fall through to any-branch
        end

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
            ishapesDir = fullfile(fileparts(simDir),'src','initial-shapes');
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

    function v = fget_str(e, a, b)
        if isstruct(e) && isfield(e,a) && isstruct(e.(a)) && isfield(e.(a),b)
            v = string(e.(a).(b));
        else
            v = "";
        end
    end

end