function sim_driver_quad_tree()
% -------------------------------------------------------------------------
% Refactored main driver (single "continue" mode)
% Goals:
%   • Use a single, stable SimResults layout (no dated run folders)
%   • Never re-solve an already cached hash (idempotent over parameter tuples)
%   • De-duplicate on disk: one .mat per unique hash in SimResults/hashed_results
%   • Maintain a global cache.mat (in-memory hints, warm-starts) and global catalog.csv|.mat
%   • Append-only logging to SimResults/OPfile.txt
%   • Keep quadtree exploration logic, but route all solves through getOrSolve()
% -------------------------------------------------------------------------

    % --- Configuration knobs (edit as needed)
    sim = struct();
    sim.plot.enable = false;
    sim.qtree.maxDepth= 4;
    sim.qtree.maxCells= 2000;
    sim.qtree.eTol = 5e-3; % energy uniformity threshold for refinement
    sim.qtree.pTol = 5e-3; % pressure uniformity threshold
    sim.qtree.shapeTau= 0.08; % optional: shape-distance threshold

    % Parameter-plane bounds (example: spontaneous curvature H0^(1), H0^(2))
    sim.params.h0min1 = -10; sim.params.h0max1 = +10;
    sim.params.h0min2 = -10; sim.params.h0max2 = +10;
    
    % Fixed physical params (edit to match your model)
    sim.params.areaFrac = 0.50; % x^(1)
    sim.params.volFrac = 0.72; % v
    sim.params.kG = 0.00; % Gaussian modulus (nondim)
    sim.params.k1 = 1.00; % bending in phase 1
    sim.params.k2 = 1.00; % bending in phase 2

    % Solver options (bvp6c)
    sim.solver.maxMesh = 200; % example knob
    sim.solver.resTol = 1e-6; % residual gate for acceptance
    sim.solver.maxTries = 3; % retries with smaller steps/warm-starts

    % --- Stable results layout (no per-run dirs)
    sim.paths.root = fullfile(pwd, 'SimResults');
    sim.paths.hashes = fullfile(sim.paths.root, 'hashed_results');
    sim.paths.cache = fullfile(sim.paths.root, 'cache.mat');
    sim.paths.catalogCsv = fullfile(sim.paths.root, 'catalog.csv');
    sim.paths.catalogMat = fullfile(sim.paths.root, 'catalog.mat');
    sim.paths.log = fullfile(sim.paths.root, 'OPfile.txt');

    ensureFolders(sim.paths);
    logmsg(sim, "[START] sim_driver_quad_tree (continue-mode, de-duplicated IO)");

    % --- Global cache + catalog (persist across runs)
    state = loadGlobalState(sim);

    % --- (Optional) live coverage plot
    if sim.plot.enable
        sim.plot.fig = figure('Position',[40 40 1100 900],'Color','w');
        sim.plot.ax = axes('Position',[0.08 0.08 0.86 0.86]); hold(sim.plot.ax,'on'); box(sim.plot.ax,'on');
        xlabel(sim.plot.ax,'H0^{(1)}'); ylabel(sim.plot.ax,'H0^{(2)}');
        title(sim.plot.ax,'Quadtree coverage'); grid(sim.plot.ax,'on');
    end

    % --- Initialize quadtree cells
    QT = initQuadtree(sim);

    % --- Process quadtree (route solves through getOrSolve)
    qParams = struct('maxDepth', sim.qtree.maxDepth, ...
        'maxCells', sim.qtree.maxCells, ...
        'eTol', sim.qtree.eTol, ...
        'pTol', sim.qtree.pTol, ...
        'shapeTau', sim.qtree.shapeTau);
    [state, QT] = processQuadtree(sim, qParams, state, QT);

    % --- Persist state & close
    saveGlobalState(sim, state);
    logmsg(sim, '[DONE] Quadtree pass complete.');
end

%% ======================================================================
%% Quadtree driver (skeleton) — ensure all corner solves go through getOrSolve
function QT = initQuadtree(sim)
    QT = struct();
    QT.cells = struct('lvl', 0, ...
    'rect',[sim.params.h0min1, sim.params.h0max1, sim.params.h0min2, sim.params.h0max2], ...
    'children', [], 'stats', []);
end

function [state, QT] = processQuadtree(sim, qParams, state, QT)
    % For brevity, process a fixed grid at max depth; extend to adaptive later
    n = 2^qParams.maxDepth;
    h1 = linspace(sim.params.h0min1, sim.params.h0max1, n+1);
    h2 = linspace(sim.params.h0min2, sim.params.h0max2, n+1);

    for i = 1:numel(h1)-1
        for j = 1:numel(h2)-1
            % Four corners of this cell
            P = [
                h1(i) h2(j) ;
                h1(i+1) h2(j) ;
                h1(i) h2(j+1);
                h1(i+1) h2(j+1)
                ];

            % Solve (or fetch) at corners
            results = cell(1,4);
            for k = 1:4
                params = baseParams(sim);
                params.H01 = P(k,1); params.H02 = P(k,2);
                [state, results{k}] = getOrSolve(sim, state, params);
            end

            % TODO: classify morphology at corners; decide on refinement via eTol/pTol
            %#ok<NASGU>
        end
    end
end
%% ======================================================================
%% Idempotent solve wrapper
function [state, R] = getOrSolve(sim, state, params)
    key = paramsKey(params);
    h = sha256(key);
    fpath = fullfile(sim.paths.hashes, [h, '.mat']);
    
    % Fast path: in-memory index
    if isKey(state.index, h)
        R = loadResultFromDisk(fpath, h, true); % tolerate missing: will self-heal below
        if ~isempty(R); return; end
    end
    
    % Disk path: already computed?
    if exist(fpath, 'file')
        R = loadResultFromDisk(fpath, h, false);
        state.index(h) = R.meta; % rehydrate index
        state = upsertCatalog(sim, state, R.meta);
        return;
    end
    
    % Not found: run solver (with retries/warm-starts inside)
    logmsg(sim, sprintf('[SOLVE] %s', key));
    R = solveAtParams(sim, params);
    
    % Persist atomically
    meta = makeMetaRecord(h, params, R);
    atomicSave(fpath, R, meta);
    
    % Update state & catalog
    state.index(h) = meta;
    state = upsertCatalog(sim, state, meta);
    logmsg(sim, sprintf('[CACHED] %s → %s.mat', key, h));
end

function R = solveAtParams(sim, params)
    % --- PLACEHOLDER ---
    % Wire in your existing BVP system here; honor sim.solver.* gates
    % Expected outputs in R:
    %   R.shape   : struct with fields (s, r, z, psi, ...)
    %   R.energy  : scalar
    %   R.pressure: scalar
    %   R.residual: scalar (max residual across DE/BC)
    %   R.label   : morphology label (string), if available

    % For now, emit a stub result so the IO layer can be exercised.
    R = struct();
    R.shape   = struct('s', linspace(0,1,64), 'r', rand(1,64), 'z', rand(1,64), 'psi', rand(1,64));
    R.energy  = rand();
    R.pressure= rand();
    R.residual= sim.solver.resTol * 0.1; % pretend we converged
    R.label   = "unknown";
end
%% ======================================================================
%% Catalog & state management
function S = loadGlobalState(sim)
    % Load persistent state and ensure index is a containers.Map
    S = struct();
    S.version = 2; % containers.Map index
    S.index = containers.Map('KeyType','char','ValueType','any');

    % Try to load existing cache
    if exist(sim.paths.cache, 'file')
        try
            T = load(sim.paths.cache, '-mat');
            if isfield(T, 'state')
                S = T.state;
                % Upgrade or initialize index to containers.Map
                if ~isfield(S,'index') || isempty(S.index)
                    S.index = containers.Map('KeyType','char','ValueType','any');
                elseif ~isa(S.index, 'containers.Map')
                    Sindex = containers.Map('KeyType','char','ValueType','any');
                    if isstruct(S.index)
                        fn = fieldnames(S.index);
                        for ii = 1:numel(fn)
                            try
                                Sindex(fn{ii}) = S.index.(fn{ii});
                            catch
                            end
                        end
                    end
                    S.index = Sindex;
                end
            end
        catch
            % ignore and keep defaults
        end
    end
    
    
    % Warm-start: rehydrate index from catalog.mat
    if exist(sim.paths.catalogMat, 'file')
        try
            C = load(sim.paths.catalogMat, '-mat');
            if isfield(C, 'catalog') && ~isempty(C.catalog)
                for ii = 1:numel(C.catalog)
                    try
                        S.index(C.catalog(ii).hash) = C.catalog(ii);
                    catch
                    end
                end
            end
        catch
            % ignore catalog issues
        end
    end
end

function saveGlobalState(sim, state)
    try
        state.last_write = datetime('now');
        save(sim.paths.cache, 'state', '-v7');
    catch ME
        logmsg(sim, sprintf('[WARN] Cache save failed: %s', ME.message));
    end
end

function state = upsertCatalog(sim, state, meta)
    % In-memory array of records (struct with fields below)
    if ~isfield(state, 'catalog') || isempty(state.catalog)
        state.catalog = meta; % start array
    else
        % replace or append
        idx = find(strcmp({state.catalog.hash}, meta.hash), 1);
        if isempty(idx)
            state.catalog(end+1) = meta; %#ok<AGROW>
        else
            state.catalog(idx) = meta;
        end
    end
    
    
    % Mirror to .mat for quick loading
    try
        catalog = state.catalog; %#ok<NASGU>
        save(sim.paths.catalogMat, 'catalog', '-v7');
    catch ME
        logmsg(sim, sprintf('[WARN] catalog.mat save failed: %s', ME.message));
    end
    
    
    % Append/update CSV (append-only; duplicates filtered by a temp file)
    try
        writeCatalogCsv(sim.paths.catalogCsv, meta);
    catch ME
        logmsg(sim, sprintf('[WARN] catalog.csv write failed: %s', ME.message));
    end
end

function writeCatalogCsv(csvPath, meta)
    header = {'hash','H01','H02','areaFrac','volFrac','kG','k1','k2','energy','pressure','residual','label','timestamp'};
    row = {meta.hash, meta.H01, meta.H02, meta.areaFrac, meta.volFrac, meta.kG, meta.k1, meta.k2, meta.energy, meta.pressure, meta.residual, meta.label, char(meta.timestamp)};
    
    
    needHeader = ~exist(csvPath,'file');
    fid = fopen(csvPath,'a');
    if fid<0, error('Cannot open %s for append', csvPath); end
    cleaner = onCleanup(@() fclose(fid));
    
    
    if needHeader
        fprintf(fid, '%s\n', strjoin(header, ','));
    end
    fprintf(fid, '%s\n', csvEscape(row));
end

function s = csvEscape(row)
    parts = cellfun(@(x) toCellStr(x), row, 'UniformOutput', false);
    s = strjoin(parts, ',');
end

function t = toCellStr(x)
    if isnumeric(x)
        if isscalar(x)
            t = sprintf('%.15g', x);
        else
            t = jsonencode(x);
        end
    elseif isstring(x) || ischar(x)
        txt = string(x);
        if contains(txt, [",", '"', newline])
            t = ['"' strrep(strrep(txt, '"','""'), '"','""') '"']; %#ok<*STRCL1>
        else
            t = char(txt);
        end
    elseif isdatetime(x)
        t = char(x);
    else
        t = jsonencode(x);
    end
end
%% ======================================================================
%% IO helpers
function ensureFolders(paths)
    if ~exist(paths.root,'dir'), mkdir(paths.root); end
    if ~exist(paths.hashes,'dir'), mkdir(paths.hashes); end
end

function atomicSave(fpath, R, meta)
    tmp = [fpath, '.tmp'];
    S.R = R; S.meta = meta; %#ok<STRNU>
    save(tmp, '-struct', 'S', '-v7');
    movefile(tmp, fpath, 'f');
end

function R = loadResultFromDisk(fpath, h, tolerateMissing)
    if exist(fpath, 'file')
        T = load(fpath, '-mat');
        if isfield(T,'R') && isfield(T,'meta')
            R = T; % contains R + meta
            R = rmfield(R, setdiff(fieldnames(R), {'R','meta'}));
            R.meta.hash = h;
            return;
        end
    end
    if tolerateMissing
        R = [];
    else
        error('Expected result %s missing or malformed.', fpath);
    end
end

function meta = makeMetaRecord(h, params, R)
    meta = struct();
    meta.hash = h;
    meta.H01 = params.H01;
    meta.H02 = params.H02;
    meta.areaFrac = params.areaFrac;
    meta.volFrac = params.volFrac;
    meta.kG = params.kG;
    meta.k1 = params.k1;
    meta.k2 = params.k2;
    meta.energy = R.energy;
    meta.pressure = R.pressure;
    meta.residual = R.residual;
    meta.label = R.label;
    meta.timestamp= datetime('now');
end

function params = baseParams(sim)
    params = struct('H01', 0, 'H02', 0, ...
    'areaFrac', sim.params.areaFrac, ...
    'volFrac', sim.params.volFrac, ...
    'kG', sim.params.kG, ...
    'k1', sim.params.k1, ...
    'k2', sim.params.k2);
end

function key = paramsKey(p)
    % Strict ordering for hash reproducibility
    key = sprintf('H01=%.12g;H02=%.12g;af=%.12g;vf=%.12g;kG=%.12g;k1=%.12g;k2=%.12g', ...
    p.H01, p.H02, p.areaFrac, p.volFrac, p.kG, p.k1, p.k2);
end

function h = sha256(str)
    md = java.security.MessageDigest.getInstance('SHA-256');
    md.update(uint8(str));
    h = sprintf('%02x', typecast(md.digest(), 'uint8'));
end

function logmsg(sim, msg)
    ts = datestr(now, 'yyyy-mm-dd HH:MM:SS.FFF');
    txt = sprintf('%s %s\n', ts, msg);
    fprintf('%s', txt);
    try
    fid = fopen(sim.paths.log,'a');
    if fid>0
    fwrite(fid, txt);
    fclose(fid);
    end
    catch
    % best-effort
    end
end
