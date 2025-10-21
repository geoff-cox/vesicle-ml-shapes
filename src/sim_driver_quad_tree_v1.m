function sim_driver_quad_tree(varargin)
% SIM_DRIVER_QUAD_TREE  Unified, MAT-only driver (single "continue" mode).
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
%
% Optional name-value args:
%   'MaxIters'     (int)   : safety cap per run (default 1e9)
%   'ModelVersion' (string): ABI tag to bump when ODE/BC/state changes (default "BVP-v3.1")
%   'LogToFile'    (logical): log to SimResults/OPfile.txt (default true)

    % --------- parse options ----------
    p = inputParser;
    addParameter(p,'MaxIters', 1e9, @(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(p,'ModelVersion', "BVP-v3.1", @(s)isstring(s)||ischar(s));
    addParameter(p,'LogToFile', true, @(b)islogical(b)&&isscalar(b));
    parse(p, varargin{:});
    opt = p.Results;
    model_version = string(opt.ModelVersion);
    
    % --------- paths ----------
    projRoot = fileparts(fileparts(mfilename('fullpath'))); % src/ -> project root
    simDir   = fullfile(projRoot,'SimResults');
    solDir   = fullfile(simDir,'hashed_results');
    if ~exist(simDir,'dir');   mkdir(simDir);   end
    if ~exist(solDir,'dir');   mkdir(solDir);   end
    opfile   = fullfile(simDir,'OPfile.txt');
    
    % --------- logging helper ----------
    function say(fmt, varargin)
        msg = sprintf(fmt, varargin{:});
        fprintf('%s\n', msg);
        if opt.LogToFile
            try logmsg(opfile, msg); catch, end %#ok<CTCH>
        end
    end
    
    say('[driver] start | version=%s | maxIters=%g', model_version, opt.MaxIters);
    
    % --------- load catalog + cache ----------
    T = catalog_load(simDir);              % utils: robust to legacy/wide forms
    cache = loadCacheIfExists(simDir);     % utils: returns struct or fresh default
    
    % Ensure cache has a frontier if empty (your initQuadtree may do this inside processQuadtree)
    if ~isfield(cache,'frontier'); cache.frontier = []; end
    
    % --------- main loop ----------
    iters = 0;
    while iters < opt.MaxIters
        iters = iters + 1;
    
        % Ask quadtree for the next task (contract: returns struct with .params or [])
        [task, cache] = processQuadtree(cache, T);
        if isempty(task)
            say('[driver] no pending tasks; exiting.');
            break
        end
        if ~isfield(task,'params') || isempty(task.params)
            say('[driver] task missing params; skipping.');
            cache = saveAndContinue(cache); %#ok<NASGU>
            continue
        end
        params = task.params;
    
        % Stable content hash for this parameter set (include model ABI)
        % simpleDataHash signature in your utils may be (params, model_version) or (params) — pass both safely
        try
            hash = simpleDataHash(params, model_version);
        catch
            hash = simpleDataHash(params); % fallback if utils expect one arg
        end
        fmat = fullfile(solDir, hash + ".mat");
    
        % Fast path: already solved => ensure catalog row exists/updated
        if exist(fmat, 'file')
            say('[skip] %s (exists)', hash);
            try
                S = load(fmat, 'meta'); % prefer rich metadata saved with the result
                entry = struct('params', params, 'meta', S.meta);
            catch
                entry = struct('params', params, 'meta', struct('hash',hash,'version',model_version));
            end
            T = catalog_append(simDir, hash, entry); % utils: MAT-only append/merge
            cache = saveAndContinue(cache);
            continue
        end
    
        % Optional: warm start suggestion from cache or catalog
        warm = struct(); % leave empty unless your solveAtParams uses it
    
        % Solve
        try
            say('[solve] %s', hash);
            [result, meta] = solveAtParams(params, warm);  % utils: your bvp6c pipeline
            % ensure required bits on meta
            meta.hash    = hash;
            if ~isfield(meta,'version') || strlength(string(meta.version))==0
                meta.version = model_version;
            end
        catch ME
            say('[error] %s :: %s', hash, ME.message);
            % Optionally mark failure on the cache so the frontier advances reasonably
            if ~isfield(cache,'failures'), cache.failures = {}; end
            cache.failures{end+1} = struct('hash',hash,'msg',ME.message); %#ok<AGROW>
            cache = saveAndContinue(cache);
            continue
        end
    
        % Atomic write
        tmp = fmat + ".tmp";
        save(tmp, 'result','meta','-v7.3');
        movefile(tmp, fmat, 'f');
    
        % Catalog append (MAT-only)
        entry = struct('params', params, 'meta', meta);
        T = catalog_append(simDir, hash, entry); % utils handle timestamp merge
    
        % Advance + persist cache
        cache = saveAndContinue(cache);
    end
    
    say('[driver] done | iterations=%d | catalogRows=%d', iters, height(catalog_load(simDir)));

    % --------- nested helper ---------
    function cache2 = saveAndContinue(cache1)
        % Let your quadtree/process logic update cache/frontier as needed upstream.
        saveCache(simDir, cache1);      % utils: atomic save to SimResults/cache.mat
        cache2 = cache1;
    end
end
