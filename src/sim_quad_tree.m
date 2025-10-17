% File: src/sim_driver_quad_tree.m
% Summary of changes:
%   - Removes per-run dated folders and modes; always-continue is default
%   - Uses content-addressed results (hash) to skip re-solves
%   - Appends to catalog.csv and keeps a MATLAB mirror catalog.mat
%   - Maintains a small cache.mat (quadtree frontier + warm starts)

function sim_quad_tree()
    % --- setup ---
    simroot = project_root();                 % or hardcode if you prefer
    simdir  = fullfile(simroot,'SimResults');
    ensure_dir(fullfile(simdir,'solutions'));
    catMat  = fullfile(simdir,'catalog.mat');
    cacheF  = fullfile(simdir,'cache.mat');

    % Load catalog mirror
    T = load_or_build_catalog(catMat);

    % Load quadtree/warm cache
    cache = load_cache(cacheF);

    % Model ABI/version — bump this when ODE/BC/state layout changes
    sim.model_version = "BVP-v3.1";
    sim.solver        = @solveAtParams;  %#ok<NASGU> placeholder ref

    % --- main exploration loop (pseudo-quadtree) ---
    while true
        task = next_task_from_quadtree(cache, T);
        if isempty(task)
            fprintf('[driver] No pending tasks. Exiting.\n');
            break;
        end

        params = task.params;

        % Fast skip if already solved
        hash = simpleDataHash(params, sim.model_version);
        fmat = result_file(simdir, hash);
        if exist(fmat,'file')
            % Ensure catalog row exists
            if ~is_hash_in_catalog(T, hash)
                % load meta to fill summary columns
                S = load(fmat,'meta');
                T = catalog_append(catMat, T, hash, params, S.meta);
            end
            % Mark task as done and continue
            cache = mark_done(cache, task);
            save_cache(cacheF, cache);
            continue;
        end

        % Warm start (optional): pick nearest neighbor by (H0_1,H0_2,...) from T
        warm = pick_warm_start(params, T, simdir);

        % Solve
        try
            [result, meta] = solveAtParams(params, warm); % your existing bvp6c wrapper
            meta.version = sim.model_version;
            meta.hash    = hash;
        catch ME
            cache = mark_failed(cache, task, ME);
            save_cache(cacheF, cache);
            continue;
        end

        % Atomic write
        tmp = fmat + ".tmp";
        save(tmp, 'result','meta','-v7.3');
        movefile(tmp, fmat, 'f');

        % Append catalog
        T = catalog_append(catCsv, catMat, T, hash, params, meta);

        % Advance quadtree / bookkeeping
        cache = mark_done(cache, task);
        save_cache(cacheF, cache);
    end
end

% ---------------- local helpers (keep in this file or split into a util) ----------------

function root = project_root()
    root = fileparts(fileparts(mfilename('fullpath'))); % src/ -> project_root
end

function ensure_dir(p)
    if ~exist(p,'dir'); mkdir(p); end
end

function f = result_file(simdir, hash)
    f = fullfile(simdir,'solutions', hash + ".mat");
end

function T = load_or_build_catalog(catMat)
    if exist(catMat,'file')
        S = load(catMat,'T');
        T = S.T;
        return;
    end
    % empty
    T = cell2table(cell(0,15), 'VariableNames', ...
        ["hash","timestamp","H0_1","H0_2","x1","v","k1","k2","kG","energy","pressure","residualMax","meshN","shape","model_version"]);
end

function cache = load_cache(cacheF)
    if exist(cacheF,'file'); S = load(cacheF); cache = S.cache; else; cache = struct(); end
    if ~isfield(cache,'frontier'); cache.frontier = []; end
end

function save_cache(cacheF, cache)
    tmp = cacheF + ".tmp";
    save(tmp, 'cache','-v7');
    movefile(tmp, cacheF, 'f');
end

function tf = is_hash_in_catalog(T, hash)
    if isempty(T); tf = false; return; end
    tf = any(T.hash == string(hash));
end

function T = catalog_append(catCsv, catMat, T, hash, params, meta)
    row = table( ...
        string(hash), datetime('now','Timezone','UTC','Format','yyyy-MM-dd HH:mm:ss'), ...
        field_or_nan(params,'H0_1'), field_or_nan(params,'H0_2'), field_or_nan(params,'x1'), field_or_nan(params,'v'), ...
        field_or_nan(params,'k1'), field_or_nan(params,'k2'), field_or_nan(params,'kG'), ...
        field_or_nan(meta,'energy'), field_or_nan(meta,'pressure'), field_or_nan(meta,'residualMax'), field_or_nan(meta,'meshN'), ...
        field_or_string(meta,'shapeLabel'), field_or_string(meta,'version'), ...
        'VariableNames', ["hash","timestamp","H0_1","H0_2","x1","v","k1","k2","kG","energy","pressure","residualMax","meshN","shape","model_version"]);
    if ~is_hash_in_catalog(T, hash)
        % append to CSV on disk
        if ~exist(catCsv,'file'); writetable(row, catCsv);
        else; writetable(row, catCsv, 'WriteMode','Append');
        end
        % append in-memory
        T = [T; row];
        save(catMat, 'T');
    end
end

function v = field_or_nan(S,k)
    if isfield(S,k) && ~isempty(S.(k)) && isnumeric(S.(k)), v = double(S.(k)); else, v = NaN; end
end

function v = field_or_string(S,k)
    if isfield(S,k) && ~isempty(S.(k)), v = string(S.(k)); else, v = ""; end
end

function hash = simpleDataHash(params, model_version)
    key = struct('solverABI', string(model_version));
    keys = {'H0_1','H0_2','x1','v','k1','k2','kG'};
    for i=1:numel(keys)
        k = keys{i};
        if isfield(params,k), key.(k) = params.(k); end
    end
    bytes = getByteStreamFromArray(key);
    algorithm = 'MD5'; % Other options: SHA-256
    hash = hash_bytes(bytes, algorithm);

end

function h = hash_bytes(bytes, algo)
    if nargin < 2, algo = 'MD5'; end
    md = java.security.MessageDigest.getInstance(algo);
    md.update(uint8(bytes));
    raw = typecast(md.digest, 'uint8');
    h = lower(reshape(dec2hex(raw)',1,[]));
end

function warm = pick_warm_start(params, T, simdir)
    % Optional: nearest neighbor on selected parameter space; returns struct or []
    warm = [];
    if isempty(T); return; end
    need = ["H0_1","H0_2","x1","v"];
    hasP = all(isfield(params, cellstr(need)));
    if ~hasP; return; end
    R = rmmissing(T(:, ["hash","H0_1","H0_2","x1","v"]));
    if isempty(R); return; end
    q = [params.H0_1, params.H0_2, params.x1, params.v];
    M = [R.H0_1, R.H0_2, R.x1, R.v];
    [~,ix] = min(vecnorm(M - q,2,2));
    cand = R.hash(ix);
    fmat = fullfile(simdir,'solutions', cand + ".mat");
    if exist(fmat,'file')
        S = load(fmat,'result','meta');
        warm = S; % or a subset needed by your solver
    end
end

function cache = mark_done(cache, task)
    % remove from frontier; advance children, etc.
    cache.frontier = setdiff(cache.frontier, task.id);
end

function cache = mark_failed(cache, task, ME)
    if ~isfield(cache,'failures'); cache.failures = struct('id',{},'msg',{}); end
    cache.failures(end+1) = struct('id', task.id, 'msg', ME.message);
end

function task = next_task_from_quadtree(cache, T)
    % Placeholder: implement your frontier/queue logic.
    % Return [] to stop.
    task = [];
    if ~isempty(cache.frontier)
        task = cache.frontier(1);
        % For demo, attach params
        task.params = cache.frontier(1).params;
    end
end

% ---------------- your existing solver wrapper ----------------
function [result, meta] = solveAtParams(params, warm)
    %#ok<INUSD>
    % Call your bvp6c pipeline here. Must return:
    %   result: struct with fields {s, y, bcResidual, deResidual}
    %   meta:   struct with quick scalars {energy, pressure, residualMax, meshN, shapeLabel}
    error('solveAtParams not wired in this stub. Connect to your existing implementation.');
end
