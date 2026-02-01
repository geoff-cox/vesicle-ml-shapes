% function sim_explore_H0_quad_tree(sim)
% % SIM_EXPLORE_H0_QUAD_TREE  Driver for exploring (H0_1,H0_2) with fixed physics in sim.MP
% %
% % New: failure registry with backoff to avoid repeated failed solves.
% 
%     % ---------------- paths & logging ----------------
%     projRoot = fileparts(fileparts(mfilename('fullpath')));
%     simDir   = fullfile(projRoot,'sim-results');
%     solDir   = fullfile(simDir,'hashed_results');
%     if ~exist(simDir,'dir'), mkdir(simDir); end
%     if ~exist(solDir,'dir'), mkdir(solDir); end
%     opfile   = fullfile(simDir,'OPfile.txt');
% 
%     say = @(fmt,varargin) fprintf([fmt '\n'], varargin{:});
%     if sim.SP.LogToFile, try logmsg(opfile, '[start]'); catch, end, end
% 
%     % ---------------- load state ----------------
%     T     = catalog_load(simDir);
%     cache = loadCacheIfExists(simDir);
%     if ~isfield(cache,'frontier'), cache.frontier = []; end
% 
%     % Failure registry shape (per-hash):
%     %   count: # failures so far
%     %   lastIter: last iteration index we tried
%     %   blockUntil: iteration index before which we will skip this hash
%     %   lastMsg: last error message
%     if ~isfield(cache,'failures') || ~isstruct(cache.failures)
%         cache.failures = struct('hash',{},'count',{},'lastIter',{},'blockUntil',{},'lastMsg',{});
%     end
% 
%     % Policy knobs (you can move these to sim.TH if you prefer):
%     MAX_RETRIES   = 1;             % after this, skip permanently (until manual clear)
%     BACKOFF_STEPS = [0 1 2 4 8];   % cool-down (in iterations) after kth failure
%     MIN_BACKOFF   = 1;             % enforce at least this many iterations of cooldown
% 
%     say('[driver] start | version=%s | maxIters=%g', sim.SP.ModelVersion, sim.SP.MaxIters);
% 
%     % ---------------- main loop ----------------
%     iters = 0;
%     while iters < sim.SP.MaxIters
%         iters = iters + 1;
% 
%         % scheduler selects next (H0_1,H0_2) task
%         [task, cache] = processQuadtree(cache, T, sim.MP);
%         if isempty(task)
%             say('[driver] no pending tasks; exiting.'); 
%             break;
%         end
% 
%         % Build physics-aware hash key
%         H0  = [task.params.H0_1, task.params.H0_2];
%         key = struct('model_version', string(sim.SP.ModelVersion), ...
%                      'H0_1', H0(1), 'H0_2', H0(2), ...
%                      'A', sim.MP.A, 'V', sim.MP.V, ...
%                      'KA', sim.MP.KA, 'KB', sim.MP.KB, 'KG', sim.MP.KG);
%         hash = simpleDataHash(key, 'SHA-256');
%         fmat = fullfile(solDir, hash + ".mat");
% 
%         % --------- 0) Skip if result already exists ----------
%         if exist(fmat,'file')
%             say('[skip] %s (exists)', hash);
%             try S = load(fmat,'meta'); catch, S = struct('meta',struct()); end
%             entry = struct('params', merge_params(task.params, sim.MP), 'meta', S.meta);
%             T = catalog_append(simDir, hash, entry);
%             saveCache(simDir, cache);
%             continue
%         end
% 
%         % --------- 1) Consult failure registry (cooldown / cap) ----------
%         [hasRec, recIdx] = failure_lookup(cache.failures, hash);
%         if hasRec
%             rec = cache.failures(recIdx);
%             % permanent skip
%             if rec.count >= MAX_RETRIES
%                 say('[skip] %s (permanent after %d failures: "%s")', hash, rec.count, scrub(rec.lastMsg));
%                 % We still record in catalog that this point is "skipped" if you want:
%                 % (comment out if you don't want catalog noise)
%                 % metaSkip = struct('status',"failed-permanent",'retries',rec.count,'msg',rec.lastMsg);
%                 % T = catalog_append(simDir, hash, struct('params',merge_params(task.params,sim.MP),'meta',metaSkip));
%                 saveCache(simDir, cache);
%                 continue
%             end
%             % cooldown skip
%             if iters < rec.blockUntil
%                 say('[skip] %s (cooldown: wait %d iters, last="%s")', hash, rec.blockUntil - iters, scrub(rec.lastMsg));
%                 saveCache(simDir, cache);
%                 continue
%             end
%         end
% 
%         % --------- 2) Warm start (from catalog solves or seeds) ----------
%         warm = pickWarmStart(task.params, sim, simDir, T);
% 
%         % --------- 3) Solve ----------
%         try
%             [result, meta] = solveAtParams(task.params, sim, warm);
%             meta.hash = hash;
%             if ~isfield(meta,'version') || strlength(string(meta.version))==0
%                 meta.version = string(sim.SP.ModelVersion);
%             end
% 
%             % on success, clear any failure record for this hash
%             if hasRec
%                 cache.failures(recIdx) = [];   % remove entry
%             end
% 
%         catch ME
%             say('[error] %s :: %s', hash, ME.message);
% 
%             % update failure registry with backoff
%             if hasRec
%                 rec.count     = rec.count + 1;
%                 rec.lastIter  = iters;
%                 cool          = BACKOFF_STEPS(min(rec.count, numel(BACKOFF_STEPS)));
%                 rec.blockUntil= iters + max(cool, MIN_BACKOFF);
%                 rec.lastMsg   = ME.message;
%                 cache.failures(recIdx) = rec;
%             else
%                 rec = struct('hash',hash, 'count',1, 'lastIter',iters, ...
%                              'blockUntil', iters + max(BACKOFF_STEPS(1),MIN_BACKOFF), ...
%                              'lastMsg', ME.message);
%                 cache.failures(end+1) = rec; %#ok<AGROW>
%             end
% 
%             saveCache(simDir, cache);
%             continue
%         end
% 
%         % --------- 4) Persist success ----------
%         tmp = fmat + ".tmp";
%         save(tmp, 'result','meta','-v7.3');
%         movefile(tmp, fmat, 'f');
% 
%         % catalog append with full params
%         entry = struct('params', merge_params(task.params, sim.MP), 'meta', meta);
%         T = catalog_append(simDir, hash, entry);
% 
%         saveCache(simDir, cache);
%     end
% 
%     say('[driver] done | iterations=%d | catalogRows=%d', iters, height(catalog_load(simDir)));
% end
% 
% % ---------- helpers (local to driver) ----------
% 
% function [found, idx] = failure_lookup(F, hash)
%     if isempty(F), found = false; idx = 0; return; end
%     idx = find(strcmp({F.hash}, hash), 1, 'first');
%     found = ~isempty(idx);
%     if ~found, idx = 0; end
% end
% 
% function s = scrub(msg)
%     % short, single-line reason for logs
%     s = msg;
%     if isempty(s), return; end
%     s = regexprep(s, '\s+', ' ');
%     if strlength(s) > 120, s = extractBefore(s, 121) + "..."; end
% end
% 
% function P = merge_params(h0params, MP)
%     P = h0params;
%     fields = fieldnames(MP);
%     for i=1:numel(fields), P.(fields{i}) = MP.(fields{i}); end
%     if ~isfield(P,'aS') || ~isfield(P,'bS')
%         [P.aS,P.bS] = computePhaseScales(P.A);
%     end
% end

function sim_explore_H0_quad_tree(sim)

    % -- paths & logging
    projRoot = fileparts(fileparts(mfilename('fullpath')));
    simDir   = fullfile(projRoot,'sim-results');
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
        tmp = fmat + ".tmp";
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