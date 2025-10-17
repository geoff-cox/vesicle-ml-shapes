function sim_driver_quad_tree_full()
% PHASEDIAGRAMDRIVER
% New, minimal driver for sampling the (H0_1, H0_2) plane using:
%   1) Quadtree refinement of parameter cells, and
%   2) Lightweight boundary tracing between morphology classes.
%
% Numerics: bvp6c
% ODE/BC:   BendV_Lag_EIGp_DE_lam / BendV_Lag_EIGp_BC_lam   (your existing files)

    % --- Paths (adjust if needed)
    addpath('src\bvp6c-solver');
    addpath('InitialShapes');

    % --- Build config/state and output folders/log
    sim = makeSimConfig();           % all knobs in one place (mode, grid limits, solver opts)
    sim = initFoldersAndLogging(sim);  % sets S.paths.* and S.log.fid
    sim = initResultsIO(sim);

    % --- (Optional) initialize a figure for live coverage
    if sim.plot.enable
        sim.plot.fig = figure('Position',[40 40 1100 900],'Color','w');
        sim.plot.ax  = axes('Position',[0.08 0.08 0.86 0.86]); hold(sim.plot.ax,'on'); box(sim.plot.ax,'on');
        xlabel(sim.plot.ax,'H0^{(1)}'); ylabel(sim.plot.ax,'H0^{(2)}');
        title(sim.plot.ax,'Quadtree coverage'); grid(sim.plot.ax,'on');
    end

    cache = bootstrapCache(sim);
    cache = loadCache(sim, cache);

    % --- Initialize quadtree over (H0_1, H0_2)
    QT = initQuadtree(sim);

    % --- Quadtree processing (solves corners, refines mixed cells, traces boundaries)
    qParams = struct('maxDepth', 4, ...     % refine up to 2^6 resolution per axis
                     'maxCells', 2000, ...
                     'eTol',     5e-3, ...  % energy uniformity threshold
                     'pTol',     5e-3, ...  % pressure uniformity threshold
                     'shapeTau', 0.08);     % (optional) shape-distance threshold if used
    processQuadtree(sim, QT, qParams, cache);

    % --- Done
    saveCache(sim, cache);
    logmsg(sim,'[DONE] Quadtree pass complete.');
    if sim.log.fid > 0, fclose(sim.log.fid); end
end

function sim = makeSimConfig()
    % MAKESIMCONFIG  Build configuration/state struct.

    % ---- High-level SIM params ----
    A  = 0.75;
    V  = 0.72;
    KG = 0;
    KA = 1;
    KB = 1;
    S_star = acos(1-2*A);
    aS = S_star/pi;
    bS = (S_star - pi)/pi;

    sim.params = struct(...
        'A', A, 'V', V, 'KG', KG, 'KA', KA, 'KB', KB, ...
        'aS', aS, 'bS', bS...
        );

    % Starting initial-guess file (lives in InitialShapes/)
    sim.initialGuess = sprintf('SIM_Node_%s_%s_%s_%s_%s_+00_+00.mat', ...
        num2str(round(A*100)), num2str(round(V*100)), ...
        num2str(KG), num2str(KA), num2str(KB));

    ig_src = fullfile('InitialShapes', sim.initialGuess);
    assert(exist(ig_src,'file')==2, 'Initial guess not found: %s', ig_src);

    % Mode: 'new' | 'continue' | 'rescan'
    sim.mode = 'new';           % change to 'continue' or 'rescan' as needed
    sim.contDate   = '15-Oct-2024';
    sim.rescanDate = '18-Oct-2024';

    sim.debug.short = true;  % short run mode

    sim.print2cmdwin = false;

    % Composite base title for MAT files
    sim.baseTitle = sprintf('SIM_Node_%s_%s_%s_%s_%s', ...
        num2str(round(A*100)), num2str(round(V*100)), ...
        num2str(KG), num2str(KA), num2str(KB));

    % ---- Phase-plane grid ----
    A_lim   = [-5, 5];
    B_lim   = [-5, 5];
    A_nodes = [0:A_lim(2), -1:-1:A_lim(1)];  % same ordering as original
    B_nodes = [0:B_lim(2), -1:-1:B_lim(1)];

    sim.grid = struct(...
        'A_lim', A_lim, 'B_lim', B_lim, ...
        'A_nodes',A_nodes,'B_nodes',B_nodes...
        );

    % ---- Constants / Directions ----
    sim.const.dirNames = [' E';' W';' N';' S';'NE';'SW';'NW';'SE'];
    sim.const.director = [ 1,0; -1,0; 0,1; 0,-1; 1,1; -1,-1; -1,1; 1,-1];

    % ---- Solver options ----
    opts            = bvpset('RelTol',1e-6,'Stats','on', 'NMax', 2000);
    delta           = 0.01;     % Taylor expansion width near poles
    maxArr          = 1000;
    stepTol         = 1/500;    % stopping criterion on step size magnitude
    stepGrow        = 40;       % matches your harmonic update scheme
    saveSolutions   = true;
    goodCountGrowThreshold = 100;

    sim.solver = struct(...
        'opts',                     opts, ...
        'delta',                    delta, ...
        'maxArr',                   maxArr, ...
        'stepTol',                  stepTol, ...
        'stepGrow',                 stepGrow, ...
        'saveSolutions',            saveSolutions, ...
        'goodCountGrowThreshold',   goodCountGrowThreshold ...
        );

    % ---- Plotting ----
    sim.plot.enable = false;

    % ---- Runtime mutable fields (filled later) ----
    sim.paths = struct();
    sim.log   = struct('fid',-1,'path','');
    sim.H0    = [];
end

function sim = initFoldersAndLogging(sim)
    % INITFOLDERSANDLOGGING  Prepare output directories and logfile.

    switch lower(sim.mode)
        case 'new'
            base = fullfile('SimResults', [datestr(now,'dd-mmm-yyyy'), filesep]);
            sim.paths.root    = base;
            sim.paths.mats    = fullfile(base, 'MATS', filesep);
            sim.paths.matname = sim.paths.mats;   % original var name parity
            
            if ~exist(sim.paths.mats, 'dir'), mkdir(sim.paths.mats); end
            copyfile(fullfile('InitialShapes', sim.initialGuess), ...
                     fullfile(sim.paths.mats, sim.initialGuess));
            logName = fullfile(base, 'OPfile.txt');

            sim.paths.results = fullfile(base, 'results');      % hashed solutions
            sim.paths.catalog = fullfile(base, 'catalog.csv');  % index
            if ~exist(sim.paths.results,'dir'), mkdir(sim.paths.results); end

        case 'continue'
            baseOld = fullfile('SimResults', [sim.contDate, filesep]);
            sim.paths.root = baseOld;
            sim.paths.mats = fullfile(baseOld, 'Expand_MATS', filesep);
            if ~exist(sim.paths.mats,'dir'), mkdir(sim.paths.mats); end
            copyfile(fullfile(baseOld,'MATS'), sim.paths.mats);
            logName = fullfile(baseOld, 'Expand_OPfile.txt');

        case 'rescan'
            baseOld = fullfile('SimResults', [sim.rescanDate, filesep]);
            sim.paths.root = baseOld;
            % Find incrementing ReScan_MATS_N directory
            d = dir(baseOld);
            count = 1;
            for k = 1:numel(d)
                if numel(d(k).name) > 10 && strncmp(d(k).name,'ReScan_MATS',11)
                    count = count + 1;
                end
            end
            sim.paths.mats = fullfile(baseOld, sprintf('ReScan_MATS_%d',count), filesep);
            if ~exist(sim.paths.mats, 'dir'), mkdir(sim.paths.mats); end
            if count > 1
                copyfile(fullfile(baseOld, sprintf('ReScan_MATS_%d',count-1), filesep), sim.paths.mats);
            else
                copyfile(fullfile(baseOld,'MATS', filesep), sim.paths.mats);
            end
            logName = fullfile(baseOld, 'ReScan_OPfile.txt');

        otherwise
            error('Unknown mode: %s', sim.mode);
    end

    % Open log
    sim.log.path = logName;
    sim.log.fid  = fopen(logName, 'wt');
    logmsg(sim, '*** Node Simulation ***');
    logmsg(sim, 'Phase-plane: [%g,%g] x [%g,%g] | kappa_A=%g', ...
        sim.grid.A_lim, sim.grid.B_lim, sim.params.KA);
end

% ============================
% Quadtree data structures
% ============================

function QT = initQuadtree(S)
    root = makeCell(S.grid.A_lim(1), S.grid.A_lim(2), ...
                    S.grid.B_lim(1), S.grid.B_lim(2), 0);
    QT.root  = root;
    QT.queue = {root};
    % empty struct array with same fields as a cell
    QT.cells = repmat(root, 0, 1);
end

function C = makeCell(a1, a2, b1, b2, depth)
    % A quadtree cell covering [a1,a2] x [b1,b2].
    C = struct( ...
        'a1',a1,'a2',a2, 'b1',b1,'b2',b2, ...
        'depth',depth, ...
        'corners', [a1 b1; a2 b1; a2 b2; a1 b2], ... % [SW; SE; NE; NW]
        'cornerSolved', false(4,1), ...
        'cornerKey', strings(4,1), ...
        'cornerLabel', nan(4,1), ...
        'cornerEnergy', nan(4,1), ...
        'cornerPressure', nan(4,1), ...
        'isUniform', false, ...
        'mixedEdges', zeros(0,2) ... % pairs of corner indices with differing labels
    );
end

function processQuadtree(S, QT, params, cache)
    % Main loop: pop a cell, solve corners, test stopping, subdivide if needed.
    
    % Debug overrides
    if isfield(S,'debug') && isfield(S.debug,'short') && S.debug.short
        params.maxDepth = 2;           % very shallow refinement
        params.maxCells = 30;          % small number of cells
        S.solver.opts   = bvpset(S.solver.opts,'NMax', 500,'RelTol',5e-6);
    end

    % params: configuration for stopping criteria & max depth, e.g.:
    %   params.maxDepth, params.maxCells, params.eTol, params.pTol, params.shapeTau
    maxCells = defaultArg(params,'maxCells', 2000);
    eTol     = defaultArg(params,'eTol',     5e-3);
    pTol     = defaultArg(params,'pTol',     5e-3);
    tau      = defaultArg(params,'shapeTau', 0.08);
    maxDepth = defaultArg(params,'maxDepth', 6);

    cellCount = 0;

    progressbar = char('-'*ones(1,100));
    while ~isempty(QT.queue) && cellCount < maxCells
        dx = floor((cellCount+1)/maxCells*100);
        C = QT.queue{1};  QT.queue(1) = [];
        progressbar(1:dx) = '#';
        fprintf('Progress: %s', progressbar);
        logmsg(S, 'Progress:\nFinalized Cells: %i (max = %i)\n', cellCount+1, maxCells);

        % 1) Solve corners (if needed)
        failCount = 0;
        for i = 1:4
            if ~C.cornerSolved(i)
                try
                    S.H0 = C.corners(i,:);
                    [sol, meta] = solveAtParams(S, cache);
                    C.cornerSolved(i)   = true;
                    C.cornerLabel(i)    = meta.label;
                    C.cornerEnergy(i)   = meta.E_total;
                    C.cornerPressure(i) = meta.P_osm;
                    cache = meta.cache;

                    if S.plot.enable
                        plot(S.plot.ax, S.H0(1), S.H0(2), '.', 'MarkerSize',8);
                        drawnow limitrate;
                    end
                catch
                    failCount = failCount + 1;
                    C.cornerSolved(i) = false;   % mark unsolved
                    logmsg(S, 'Corner (%+.2f,%+.2f) failed; will try after subdivision.', S.H0(1), S.H0(2));
                end
            end
        end
        
        if failCount >= 2
            % too risky to judge; subdivide immediately
            [C1,C2,C3,C4] = subdivideCell(C);
            QT.queue(end+1:end+4) = {C1,C2,C3,C4};
            continue
        end

        % 2) Stopping test
        [uniform, mixedEdges] = uniformTest(C, eTol, pTol, tau);
        C.isUniform  = uniform;
        C.mixedEdges = mixedEdges;

        % 3) If uniform or max depth reached, finalize cell
        if uniform || C.depth >= maxDepth
            QT.cells = [QT.cells; C]; 
            cellCount = cellCount + 1;
            % Optional: attempt boundary tracing from each mixed edge discovered earlier
            continue
        end

        % 4) Otherwise subdivide into 4 children
        [C1,C2,C3,C4] = subdivideCell(C);
        QT.queue{end+1} = C1; 
        QT.queue{end+1} = C2; 
        QT.queue{end+1} = C3; 
        QT.queue{end+1} = C4; 

        cellCount = cellCount + 1;

        % 5) Kick off boundary tracing from each mixed edge (lightweight)
        for e = 1:size(mixedEdges,1)
            i = mixedEdges(e,1); j = mixedEdges(e,2);
            pa = C.corners(i,:);  pb = C.corners(j,:);
            try
                traceBoundary(S, pa, pb, cache); % see block 2
            catch ME
                fprintf('Boundary trace failed at [%g,%g]-[%g,%g]: %s\n', ...
                        pa(1),pa(2),pb(1),pb(2), ME.message);
            end
        end
    end
end

function cache = bootstrapCache(S)
    cache = initSolveCache();
    seeds = [ 0 0;
              1 0; -1 0; 0 1; 0 -1;
              1 1; -1 1; 1 -1; -1 -1 ];
    for i = 1:size(seeds,1)
        S.H0 = seeds(i,:);
        try
            [~, meta] = solveAtParams(S, cache);
            cache = meta.cache;
            if S.plot.enable
                plot(S.plot.ax, S.H0(1), S.H0(2), '.', 'MarkerSize',10); drawnow limitrate;
            end
        catch ME
            fprintf('Bootstrap skip at (%g,%g): %s\n', S.H0(1), S.H0(2), ME.message);
        end
    end
end

function sol = continuationSolve(S, p0, sol0, pT, stepMax)
    % Two-leg path: (p0_1,p0_2) -> (pT_1,p0_2) -> (pT_1,pT_2)
    path = [linspace(p0(1), pT(1), max(2,ceil(abs(pT(1)-p0(1))/stepMax)))'  repmat(p0(2),[],1)];
    path = [path;  [repmat(pT(1), max(2,ceil(abs(pT(2)-p0(2))/stepMax)))'  linspace(p0(2), pT(2), max(2,ceil(abs(pT(2)-p0(2))/stepMax)))']];

    sol = sol0;
    Par = S.params;
    baseDelta = S.solver.delta;

    for k = 2:size(path,1)
        target = path(k,:);
        % adaptive attempt: try up to 5 backtracks with smaller step and larger delta
        ok = false; localStep = norm(target - path(k-1,:));
        for attempt = 1:5
            % build handles with a slightly larger delta near failures
            bump = min(3, 1 + 0.5*(attempt-1));     % 1, 1.5, 2, 2.5, 3 x
            Par.H0    = target;
            Par.delta = max(1e-3, min(0.03, baseDelta*bump));

            odefun = @(s,y,lam) BendV_Lag_EIGp_DE_impl(s,y,lam,Par);
            bcfun  = @(ya,yb,lam) BendV_Lag_EIGp_BC_impl(ya,yb,lam,Par);

            try
                sol = bvp6c(odefun, bcfun, sol, S.solver.opts);
                [BCmax,BCrep] = bc_diagnostics(sol,bcfun);
                if BCmax <= 1e-6
                    ok = true;
                    break;
                else
                    % reject and backtrack
                    target = path(k-1,:) + 0.5*(target - path(k-1,:)); % halve step
                end
            catch
                % backtrack on failure
                target = path(k-1,:) + 0.5*(target - path(k-1,:));
            end
        end
        if ~ok
            error('continuationSolve: could not reach waypoint (%.3f,%.3f)', target(1),target(2));
        end
        if S.plot.enable
            plot(S.plot.ax, Par.H0(1), Par.H0(2), '.', 'MarkerSize',6); drawnow limitrate;
        end
    end
end

function [C1,C2,C3,C4] = subdivideCell(C)
    % Split into 4 children (SW, SE, NE, NW), inherit depth+1.
    am = 0.5*(C.a1 + C.a2);
    bm = 0.5*(C.b1 + C.b2);
    d1 = C.depth + 1;

    C1 = makeCell(C.a1, am, C.b1, bm, d1); % SW
    C2 = makeCell(am,  C.a2, C.b1, bm, d1); % SE
    C3 = makeCell(am,  C.a2, bm,  C.b2, d1); % NE
    C4 = makeCell(C.a1, am, bm,  C.b2, d1); % NW
end

function [uniform, mixedEdges] = uniformTest(C, eTol, pTol, tau)
    % Decide whether a cell is "uniform" (same morphology) or "mixed".
    % Heuristics:
    %  1) all corner labels equal --> uniform
    %  2) else, if max energy/pressure difference small AND shape distances small --> uniform
    %  3) track which edges are mixed to seed boundary tracing

    labs = C.cornerLabel;
    if all(labs == labs(1))
        uniform = true; mixedEdges = zeros(0,2); return;
    end

    eSpread = max(C.cornerEnergy)   - min(C.cornerEnergy);
    pSpread = max(C.cornerPressure) - min(C.cornerPressure);

    % Optional: include shape descriptor distances if you store them in meta
    % Here, we only use labels + E/P spreads
    shapeSmall = true; % replace with real test if you export descriptors

    uniform = (eSpread < eTol) && (pSpread < pTol) && shapeSmall;

    % Which edges have differing labels?
    edges = [1 2; 2 3; 3 4; 4 1]; % SW-SE, SE-NE, NE-NW, NW-SW
    mixedEdges = [];
    if ~uniform
        for k = 1:4
            if labs(edges(k,1)) ~= labs(edges(k,2))
                mixedEdges(end+1,:) = edges(k,:); %#ok<AGROW>
            end
        end
    end
end

function v = defaultArg(s, field, vDefault)
    if isfield(s, field), v = s.(field); else, v = vDefault; end
end

% ==========================================
% Cached solves & labeling at (H0_1, H0_2)
% ==========================================

function cache = initSolveCache()
    % Simple cache for solved parameter points and their solutions/labels.
    cache = struct();
    cache.keys   = [];        % Nx2 matrix of params
    cache.labels = [];        % Nx1 labels
    cache.E      = [];        % Nx1 energy
    cache.P      = [];        % Nx1 pressure
    cache.sols   = {};        % cell array of bvp solutions (for reuse)
end

function [sol, meta] = solveAtParams(S, cache)
% Robust solve with accept/retry gates, disk cache reuse, and warm-starts.

    % ---------- thresholds ----------
    TH.BCmax      = 1e-6;     % accept gate on BC
    TH.DEmaxHard  = 2e-1;     % hard DE residual gate
    TH.rMin       = 1e-3;     % min radius away from poles (geometry sanity)
    TH.maxResBVP  = 5e-5;     % informational (bvp6c prints this)

    % ---------- fast path: exact on-disk cache ----------
    [resultPath, resultHash] = resultFileFor(S);
    initSol = [];
    if exist(resultPath,'file') == 2
        try
            L = load(resultPath,'sol','label','E_total','P_osm');  % minimal fields required
            % Re-check with current handles (delta/tols can differ run-to-run)
            Par = S.params; Par.H0 = S.H0; Par.delta = S.solver.delta;
            odefun = @(s,y,lam) BendV_Lag_EIGp_DE_impl(s,y,lam,Par);
            bcfun  = @(ya,yb,lam) BendV_Lag_EIGp_BC_impl(ya,yb,lam,Par);

            [BCnow,~]   = bc_diagnostics(L.sol, bcfun);
            [DEnow,~,~] = de_residual(L.sol, odefun);
            rmin        = local_min_radius_interior(L.sol);  % ignore poles

            if (BCnow <= TH.BCmax) && (DEnow <= TH.DEmaxHard) && (rmin >= TH.rMin)
                logmsg(S,' ... HIT disk cache %s | BC=%.2e | DE=%.2e | rmin=%.2e', ...
                    resultHash, BCnow, DEnow, rmin);

                sol = L.sol;
                % Bookkeeping into cache
                cache.keys   = [cache.keys; S.H0];
                cache.labels = [cache.labels; L.label];
                cache.E      = [cache.E; L.E_total];
                cache.P      = [cache.P; L.P_osm];
                cache.sols{end+1} = sol;

                meta = struct('key',resultHash,'label',L.label,'E_total',L.E_total,'P_osm',L.P_osm, ...
                              'cache',cache,'odefun',odefun,'bcfun',bcfun);
                return
            else
                initSol = L.sol;   % warm-start from disk even if we didn't accept directly
                logmsg(S,' ... disk cache warm-start %s | BC=%.2e | DE=%.2e | rmin=%.2e', ...
                    resultHash, BCnow, DEnow, rmin);
            end
        catch ME
            logmsg(S,' ... WARN: cache load check failed (%s); continuing without hit.', ME.message);
        end
    end

    % ---------- neighbor-based warm start (sign-aware) ----------
    p0 = []; solNN = [];
    if isempty(initSol) && ~isempty(cache.keys)
        idx   = nearest_by_sign(cache.keys, S.H0);
        p0    = cache.keys(idx,:);
        solNN = cache.sols{idx};
        initSol = solNN;
    end

    % Coordinate-wise short continuation if far from neighbor
    if ~isempty(p0) && norm(S.H0 - p0) > 0.3
        try
            initSol = continuationSolve(S, p0, solNN, S.H0, 0.15);  % stepMax=0.15
        catch
            initSol = solNN;    % fallback; attempt ladder will handle the rest
        end
    end

    % Final fallback: packaged seed
    if isempty(initSol)
        initSol = initialGuessFromFile(S, S.H0);
    end

    % ---------- attempt ladder (delta tweaks + optional looser tol) ----------
    function [odefun, bcfun] = makeHandles(delta)
        Par = S.params; Par.H0 = S.H0; Par.delta = delta;
        odefun = @(s,y,lam) BendV_Lag_EIGp_DE_impl(s,y,lam,Par);
        bcfun  = @(ya,yb,lam) BendV_Lag_EIGp_BC_impl(ya,yb,lam,Par);
    end

    attempts = [ ...
      struct('delta', S.solver.delta,                        'opts', S.solver.opts, 'label','baseline'); ...
      struct('delta', max(0.5*S.solver.delta, 1e-3),         'opts', S.solver.opts, 'label','delta/2'); ...
      struct('delta', max(0.25*S.solver.delta, 1e-3),        'opts', S.solver.opts, 'label','delta/4'); ...
      struct('delta', S.solver.delta, 'opts', bvpset(S.solver.opts,'RelTol',1e-5,'AbsTol',1e-7), 'label','looserTol') ...
    ];

    solved = false; lastErr = []; odefun = []; bcfun = [];
    key = sprintf('%+0.4f_%+0.4f', S.H0(1), S.H0(2)); %#ok<NASGU>

    for a = 1:numel(attempts)
        [odefun, bcfun] = makeHandles(attempts(a).delta);
        logmsg(S, 'Solve @ H0=(%+.3f,%+.3f) | attempt=%s | delta=%.4g | opts.RT=%.1e', ...
               S.H0(1), S.H0(2), attempts(a).label, attempts(a).delta, bvpget(attempts(a).opts,'RelTol'));
        try
            sol = bvp6c(odefun, bcfun, initSol, attempts(a).opts);

            % ---------- gates & diagnostics ----------
            [BCmax, BCrep]       = bc_diagnostics(sol, bcfun);
            [DEmax, ~, worstIdx] = de_residual(sol, odefun);
            rmin                  = local_min_radius_interior(sol);  % ignore poles

            logmsg(S,' ... mesh=%d | BC=%.2e (i=%d) | DE=%.2e (state %d) | rmin=%.2e', ...
                   numel(sol.x), BCmax, BCrep.idx, DEmax, worstIdx, rmin);

            % Optional quick-stop for debug smoke tests
            if isfield(S,'debug') && isfield(S.debug,'short') && S.debug.short
                if DEmax > 1e3
                    logmsg(S,' ... DEBUG-GATE: DE too large (%.2e); rejecting early.', DEmax);
                    error('debug_gate:DElarge','DE too large');
                end
            end

            isBCOK  = (BCmax <= TH.BCmax);
            isDEOK  = (DEmax <= TH.DEmaxHard);
            isGeomOK= (rmin  >= TH.rMin);

            if isBCOK && isDEOK && isGeomOK
                logmsg(S,' ... ACCEPT | mesh=%d | BC=%.2e | DE=%.2e (state %d) | rmin=%.2e', ...
                       numel(sol.x), BCmax, DEmax, worstIdx, rmin);
                solved = true;

                % ---------- persist accepted result ----------
                [label, E_total, P_osm] = labelFromSolution(sol);
                try
                    [savePath, saveHash] = resultFileFor(S);
                    delta_used  = attempts(a).delta;
                    solver_opts = attempts(a).opts;
                    BCmax_save  = BCmax;
                    DEmax_save  = DEmax;
                    save(savePath, 'sol','label','E_total','P_osm', ...
                         'BCmax_save','DEmax_save','delta_used','solver_opts','-v7.3');
                    writeCatalogRow(S, saveHash, label, E_total, P_osm, BCmax, DEmax, numel(sol.x));
                catch ME
                    logmsg(S,' ... WARN: could not save result/csv (%s)', ME.message);
                end
                break;
            else
                logmsg(S,' ... REJECT | BC=%.2e | DE=%.2e | rmin=%.2e', BCmax, DEmax, rmin);
                initSol = sol;  % warm-start next attempt from current iterate
            end
        catch ME
            logmsg(S, ' ... FAIL: %s', ME.message);
            lastErr = ME;
            % keep initSol; next attempt tweaks delta/tols
        end
    end

    if ~solved
        if ~isempty(lastErr), rethrow(lastErr); else, error('solveAtParams: failed to converge'); end
    end

    % ---------- bookkeeping into in-memory cache ----------
    [label, E_total, P_osm] = labelFromSolution(sol);
    cache.keys   = [cache.keys; S.H0];
    cache.labels = [cache.labels; label];
    cache.E      = [cache.E; E_total];
    cache.P      = [cache.P; P_osm];
    cache.sols{end+1} = sol;

    meta = struct('key',resultHash,'label',label,'E_total',E_total,'P_osm',P_osm, ...
                  'cache',cache,'odefun',odefun,'bcfun',bcfun);
end

function out = getOrSolve(params, sim, cache)
    h    = simpleDataHash(params, sim.model_version);
    fmat = resultFileFor(h);

    if exist(fmat,'file')
        % Already solved — fast path
        out = load(fmat,'result','meta');   % define fields below
        if ~isCataloged(h), appendCatalogRow(h, params, out.meta); end
        return
    end

    % Try warm-start from cache/catalog (nearest neighbor in (H01,H02))
    warm = pickWarmStart(params, cache);

    % Solve (your bvp6c wrapper), returning result + meta
    [result, meta] = solveAtParams(params, warm, sim.solver);

    % Write atomically
    tmp = [fmat '.tmp'];
    save(tmp, 'result','meta','-v7.3');
    movefile(tmp, fmat, 'f');

    appendCatalogRow(h, params, meta);
    out.result = result; out.meta = meta;
end

% --- helpers local to this file ------------------------------------------

function rmin = local_min_radius_interior(sol)
% Min radius away from poles (ignore s=0 and s=pi so r=0 at poles is allowed)
    rA = sol.y(4,:);   rB = sol.y(13,:);
    if numel(rA) >= 3
        rA = rA(2:end-1);
    end
    if numel(rB) >= 3
        rB = rB(2:end-1);
    end
    rmin = min([rA(:); rB(:)]);
    if isempty(rmin)
        rmin = Inf;  % degenerate mesh; don't fail geometry gate
    end
end

function saveCache(cache)
    f = fullfile('SimResults','cache.mat');
    tmp = [f '.tmp'];
    save(tmp, 'cache','-v7');
    movefile(tmp,f,'f');
end

function cache = loadCache()
    f = fullfile('SimResults','cache.mat');
    if exist(f,'file'), S = load(f); cache = S.cache; else, cache = struct(); end
end

function C = mergeCaches(C, D)
    if isempty(D) || isempty(D.keys), return; end
    for i=1:size(D.keys,1)
        k = D.keys(i,:);
        if isempty(C.keys) || ~ismembertol(k, C.keys, 1e-12, 'ByRows',true)
            C.keys   = [C.keys; k];
            C.labels = [C.labels; D.labels(i)];
            C.E      = [C.E; D.E(i)];
            C.P      = [C.P; D.P(i)];
            C.sols{end+1} = D.sols{i};
        end
    end
end

function key = makeSolveKey(S)
    % Only include things that define the physics; omit delta/tols.
    key = struct( ...
        'H0', S.H0, ...
        'A',  S.params.A, 'V',  S.params.V, ...
        'KA', S.params.KA,'KB', S.params.KB,'KG', S.params.KG, ...
        'aS', S.params.aS,'bS', S.params.bS, ...
        'codever','qtree-2025-02-15');   % bump when DE/BC change
end

function f = resultFileFor(hash)
    f = fullfile('SimResults','solutions',[hash '.mat']);
end

function writeCatalogRow(S, hash, label, E, P, bcmax, demax, mesh)
    % Assemble a single logical row (cell array)
    ts  = datestr(now,'yyyy-mm-dd HH:MM:SS');
    row = {ts, hash, ...
           S.H0(1), S.H0(2), S.params.A, S.params.V, ...
           S.params.KA, S.params.KB, S.params.KG, ...
           label, E, P, bcmax, demax, mesh};

    % ---- 1) Append to MAT table (catalog.mat)
    if isfield(S.paths,'catalog_mat') && exist(S.paths.catalog_mat,'file')==2
        L = load(S.paths.catalog_mat,'T'); 
        T = L.T;
        % Ensure same variable order
        newRow = cell2table(row, 'VariableNames', T.Properties.VariableNames);
        T = [T; newRow]; %#ok<AGROW>
        save(S.paths.catalog_mat,'T','-v7.3');
    end

    % ---- 2) Append to CSV (catalog.csv)
    if isfield(S.paths,'catalog_csv')
        fid = fopen(S.paths.catalog_csv,'at');   % creates if missing
        fmt = '%s,%s,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%d,%.16g,%.16g,%.3e,%.3e,%d\n';
        fprintf(fid, fmt, row{:});
        fclose(fid);
    end
end

function appendCatalogRow(h, p, m)
    % m: meta collected post-solve
    row = table( ...
        string(h), datetime('now','Format','yyyy-MM-dd HH:mm:ss'), ...
        p.H0_1, p.H0_2, p.x1, p.v, p.k1, p.k2, p.kG, ...
        m.energy, m.pressure, m.residualMax, m.meshN, ...
        m.shapeLabel, m.version, ...
        'VariableNames', ["hash","timestamp","H0_1","H0_2","x1","v","k1","k2","kG", ...
                          "energy","pressure","residualMax","meshN","shape","model_version"]);

    catCsv = fullfile('SimResults','catalog.csv');
    if ~exist(catCsv,'file')
        writetable(row, catCsv);
    else
        % Append if this hash is new
        if ~isCataloged(h)
            writetable(row, catCsv, 'WriteMode','Append');
        end
    end

    % Keep a MATLAB mirror (optional, fast loads)
    rebuildCatalogMatIfStale();
end

function S = initResultsIO(S)
    run_id = datestr(now,'yyyy-mm-dd_HHMMSS');
    S.paths.run    = fullfile(S.paths.root, ['run_', run_id], filesep);
    S.paths.sol_dir= fullfile(S.paths.run, 'solutions');
    if ~exist(S.paths.run,'dir'), mkdir(S.paths.run); end
    if ~exist(S.paths.sol_dir,'dir'), mkdir(S.paths.sol_dir); end

    % --- define BOTH catalogs explicitly
    S.paths.catalog_mat = fullfile(S.paths.run, 'catalog.mat');
    S.paths.catalog_csv = fullfile(S.paths.run, 'catalog.csv');

    % --- initialize MAT table with a canonical schema
    if ~exist(S.paths.catalog_mat,'file')
        T = table( ...
            strings(0,1), strings(0,1), ...           % timestamp, hash
            zeros(0,1), zeros(0,1), ...               % H0_1, H0_2
            zeros(0,1), zeros(0,1), ...               % A, V
            zeros(0,1), zeros(0,1), zeros(0,1), ...   % KA, KB, KG
            zeros(0,1,'int16'), ...                   % label
            zeros(0,1), zeros(0,1), ...               % E, P
            zeros(0,1), zeros(0,1), ...               % BCmax, DEmax
            zeros(0,1,'uint16'), ...                  % mesh
            'VariableNames', {'timestamp','hash','H0_1','H0_2','A','V','KA','KB','KG', ...
                              'label','E','P','BCmax','DEmax','mesh'});
        save(S.paths.catalog_mat,'T','-v7.3');
    end

    % --- initialize CSV with a header
    if ~exist(S.paths.catalog_csv,'file')
        hdr = {'timestamp','hash','H0_1','H0_2','A','V','KA','KB','KG','label','E','P','BCmax','DEmax','mesh'};
        fid = fopen(S.paths.catalog_csv,'wt');
        fprintf(fid, '%s\n', strjoin(hdr,','));
        fclose(fid);
    end
end

function logmsg(S, fmt, varargin)
    msg = sprintf(fmt, varargin{:});
    if S.log.fid > 0
        fprintf(S.log.fid, '%s\n', msg);
    end
    if S.print2cmdwin
        fprintf('%s\n', msg);
    end
end

function h = simpleDataHash(params, model_version)
    % Only include physically-meaningful knobs and the model ABI
    key = struct( ...
        'H01', params.H0_1, ...
        'H02', params.H0_2, ...
        'x1',  params.x1, ...
        'v',   params.v, ...
        'k1',  params.k1, ...
        'k2',  params.k2, ...
        'kG',  params.kG, ...
        'solverABI', model_version);  % bump when you change BC/ODE layout
    bytes = getByteStreamFromArray(key);
    h = lower(matlab.net.internal.hash.MD5(bytes));  % or SHA-256 if you have a helper
end

function tf = isCataloged(h)
    persistent cacheSet lastCheck;
    catCsv = fullfile('SimResults','catalog.csv');
    if isempty(cacheSet) || ~isequal(lastCheck, dir(catCsv).datenum)
        T = readtable(catCsv, 'TextType','string');
        cacheSet = containers.Map(T.hash, true(1,height(T)));
        lastCheck = dir(catCsv).datenum;
    end
    tf = isKey(cacheSet, string(h));
end

function rebuildCatalogMatIfStale()
    csv = fullfile('SimResults','catalog.csv');
    mat = fullfile('SimResults','catalog.mat');
    if ~exist(mat,'file') || dir(mat).datenum < dir(csv).datenum
        T = readtable(csv, 'TextType','string');
        save(mat, 'T');
    end
end

function y = normalizeForHash(x)
    if istable(x)
        x = table2struct(x, 'ToScalar', true);
    end
    if isstruct(x)
        if numel(x) > 1
            y = arrayfun(@normalizeForHash,x,'UniformOutput',false); return;
        end
        f = sort(fieldnames(x));
        s = struct();
        for i=1:numel(f), s.(f{i}) = normalizeForHash(x.(f{i})); end
        y = s; return
    elseif iscell(x)
        y = cellfun(@normalizeForHash,x,'UniformOutput',false); return
    elseif isstring(x)
        y = cellstr(x); return
    elseif isa(x,'categorical')
        y = struct('codes',double(x),'cats',{categories(x)},'ord',isordinal(x)); return
    else
        y = x;
    end
end

function [rMax, rComp, worstIdx] = de_residual(sol, odefun)
    % Second-order central difference on nonuniform grid, ignoring pole buffers.

    lam = sol.parameters(:);
    s   = sol.x;             % s in [0, pi]
    Y   = sol.y;
    n   = numel(s);
    if n < 3, error('de_residual: need at least 3 mesh points'); end

    % interior nodes (2..n-1)
    h   = diff(s);
    dY  = zeros(size(Y,1), n-2);
    for j = 2:n-1
        h0 = s(j)   - s(j-1);
        h1 = s(j+1) - s(j);
        dY(:,j-1) = ( -(h1/(h0*(h0+h1)))*Y(:,j-1) ...
                      + ((h1-h0)/(h0*h1))*Y(:,j) ...
                      + (h0/(h1*(h0+h1)))*Y(:,j+1) );
    end

    sC = s(2:end-1);
    YC = Y(:,2:end-1);
    
    % --- mask out pole buffers ---
    % estimate delta from first and last spacing if not known here
    % but better: pass it; we can infer from odefun only with extra plumbing.
    % We'll be conservative and drop 2*mean(h) near each pole as a buffer:
    pad = 2*mean(h);
    mask = (sC > pad) & (sC < (pi - pad));
    if ~any(mask)
        % fallback: evaluate everywhere
        mask = true(size(sC));
    end
    sC = sC(mask);
    YC = YC(:,mask);
    dY = dY(:,mask);

    F = zeros(size(YC));
    for j = 1:numel(sC)
        F(:,j) = odefun(sC(j), YC(:,j), lam);
    end

    R = dY - F;
    if isempty(R)
        rComp = zeros(size(Y,1),1);
    else
        rComp = max(abs(R),[],2);
    end
    [rMax, worstIdx] = max(rComp);
end

% ======================================================
% Boundary tracing from a mixed edge pa--pb (labels differ)
% ======================================================

function poly = traceBoundary(S, pa, pb, cache)
    % Return a polyline of boundary points between two morphology classes.
    % pa, pb: 1x2 endpoints with different labels (already solved by caller).
    % cache:  solve cache (updated in place).
    %
    % Heuristics (tune in practice):
    maxSteps = 200;
    ds       = 0.25;     % predictor step in parameter space (H0 units)
    epsBis   = 0.02;     % bisection tolerance along normal
    pad      = 0.05;     % small offset to start bracketing around prediction
    stopBox  = [S.grid.A_lim; S.grid.B_lim]; % bounds

    % 1) Find initial boundary point p* by bisection on segment pa->pb
    [pStar, cache] = edgeBisectionToSwitch(S, pa, pb, cache, epsBis);

    poly = pStar;   % Nx2 list of boundary points
    % Initial normal is along the mixed edge
    n = unit(pb - pa);              % normal ~ direction across boundary
    t = rot90ccw(n);                % tangent = +90deg rotation

    for k = 1:maxSteps
        % 2) Predictor: step along tangent
        pPred = pStar + ds * t;

        % keep inside box
        if ~inBox(pPred, stopBox), break; end

        % 3) Corrector: bracket along normal around pPred, then bisection to label switch
        [pNext, nNew, success, cache] = correctToBoundary(S, pPred, n, pad, epsBis, cache);
        if ~success, break; end

        % 4) Append and update frame
        poly(end+1,:) = pNext; %#ok<AGROW>

        % New tangent: perpendicular to updated normal
        t = rot90ccw(nNew);
        pStar = pNext;
        n = nNew;
    end
end

% edgeBisectionToSwitch:
function [pStar, cache] = edgeBisectionToSwitch(S, pa, pb, cache, epsBis)
    [la, cache] = getLabel(S, pa, cache);
    [lb, cache] = getLabel(S, pb, cache);
    if la == lb, error('edgeBisectionToSwitch: endpoints share label'); end
    a = pa; b = pb;
    for it = 1:40
        m = 0.5*(a+b);
        [lm, cache] = getLabel(S, m, cache);
        if lm == la, a = m; else, b = m; end
        if norm(a-b) < epsBis, break; end
    end
    pStar = 0.5*(a+b);
end

% correctToBoundary:
function [pNext, nOut, success, cache] = correctToBoundary(S, pPred, nIn, pad, epsBis, cache)
    success = false;
    n = unit(nIn);

    pL = pPred - pad*n;  pR = pPred + pad*n;
    if ~inBox(pL, [S.grid.A_lim; S.grid.B_lim]) || ~inBox(pR, [S.grid.A_lim; S.grid.B_lim])
        pNext = pPred; nOut = n; return;
    end

    [lL, cache] = getLabel(S, pL, cache);
    [lR, cache] = getLabel(S, pR, cache);

    if lL == lR
        n = -n;
        pL = pPred - pad*n;  pR = pPred + pad*n;
        [lL, cache] = getLabel(S, pL, cache);
        [lR, cache] = getLabel(S, pR, cache);
        if lL == lR
            pad2 = 2*pad;
            pL = pPred - pad2*n; pR = pPred + pad2*n;
            [lL, cache] = getLabel(S, pL, cache);
            [lR, cache] = getLabel(S, pR, cache);
            if lL == lR, pNext = pPred; nOut = n; return; end
        end
    end

    a = pL; b = pR; la = lL; lb = lR;
    for it = 1:40
        m = 0.5*(a+b);
        [lm, cache] = getLabel(S, m, cache);
        if lm == la, a = m; else, b = m; end
        if norm(a-b) < epsBis, break; end
    end
    pNext = 0.5*(a+b);
    nOut  = unit(b - a);
    success = true;
end

function tf = inBox(p, box)
% box = [ [Amin Amax]; [Bmin Bmax] ]
    tf = (p(1)>=box(1,1) && p(1)<=box(1,2) && p(2)>=box(2,1) && p(2)<=box(2,2));
end

function v = unit(v), n = norm(v); if n==0, return; end; v = v./n; end

function r = rot90ccw(v), r = [ -v(2), v(1) ]; end

% --- change signature to return cache too
function [lbl, cache] = getLabel(S, p, cache)
    S.H0 = p;
    [~, meta] = solveAtParams(S, cache);
    lbl   = meta.label;
    cache = meta.cache;   % propagate
end

function initSol = initialGuessFromFile(S, ~)
    p1 = fullfile(S.paths.mats, S.initialGuess);
    if exist(p1,'file')~=2
        p1 = fullfile('InitialShapes', S.initialGuess);
    end
    tmp = load(p1);
    initSol = tmp.Version(1).Solution;
end

function [label, E_total, P_osm] = labelFromSolution(sol)
    E_total = sol.y(9,end) - sol.y(18,end) + 0.5*sol.y(4,end);
    P_osm   = sol.parameters(1);
    rA = sol.y(4,:); rB = sol.y(13,:);
    rr = [rA rB];

    if exist('findpeaks','file')
        np = numel(findpeaks(rr));
    else
        % very light peak count without toolbox
        np = sum( rr(2:end-1) > rr(1:end-2) & rr(2:end-1) >= rr(3:end) );
    end

    if     np <= 2, label = 1;
    elseif np <= 4, label = 2;
    else            label = 3;
    end
end

% Sets the Boundary Conditions for the two-phase vesicle shape equations.
%
% Energy = Bending + Gauss Bending + Line Tension
%
% Constraints: Volume
% Eigenvalues: Pressure
%
% Var:   Q      H     Psi     r      z     LAM     s      V      B
%   α: y( 1)  y( 2)  y( 3)  y( 4)  y( 5)  y( 6)  y( 7)  y( 8)  y( 9)  
%   β: y(10)  y(11)  y(12)  y(13)  y(14)  y(15)  y(16)  y(17)  y(18)

function res = BendV_Lag_EIGp_BC_impl(y_poles, y_neck, lam, par)
    % -------- Simulation Parameters --------
    kA = par.KA;
    kB = par.KB;
    kG = par.KG;
    Vf = par.V;
    H0 = par.H0;

    % α-phase (s)outh pole conditions
    south_pole = num2cell(y_poles(1:9));
    [QAs, HAs, PAs, rAs, zAs, LAs, sAs, VAs, EAs] = deal(south_pole{:});
    res_south = [
        QAs
        % HAs
        PAs
        rAs
        zAs
        % LAs
        sAs
        VAs
        EAs
    ];

    % β-phase (n)orth pole conditions
    north_pole = num2cell(y_poles(10:18));
    [QBn, HBn, PBn, rBn, zBn, LBn, sBn, VBn, EBn] = deal(north_pole{:});
    res_north = [
        QBn
        % HBn
        PBn - pi
        % rBn
        % zBn
        % LBn
        sBn
        VBn
        EBn
    ];

    % α-phase nec(k) conditions
    alpha_neck = num2cell(y_neck(1:9));
    [QAk, HAk, PAk, rAk, zAk, LAk, sAk, VAk, EAk] = deal(alpha_neck{:});

    % β-phase nec(k) conditions
    beta_neck = num2cell(y_neck(10:18));
    [QBk, HBk, PBk, rBk, zBk, LBk, sBk, VBk, EBk] = deal(beta_neck{:});

    res_neck = [
        PBk - PAk
        rBk - rAk
        zBk - zAk
        (VAk - VBk) - Vf;
        (QBk - QAk) - sin(PAk)/rAk;
        kB*(2*HBk - H0(2)) ...
            - kA*(2*HAk - H0(1)) ...
            + kG*(sin(PAk)/rAk);
        kA*(2*HAk - H0(1))*(HAk - sin(PAk)/rAk + H0(1)/2) ...
        - kB*(2*HBk - H0(2))*(HBk - sin(PAk)/rAk + H0(2)/2) ...
        - (cos(PAk)/rAk) + LBk - LAk
    ];

    res = [res_south; res_north; res_neck];
end

function dyds = BendV_Lag_EIGp_DE_impl(S, y, lam, par)
    % -------- Simulation Parameters --------
    kA = par.KA;
    kB = par.KB;
    H0 = par.H0;
    aS = par.aS;
    bS = par.bS;
    delta= par.delta;

    % Scale S to each region
    SA = aS*S;
    SB = bS*S + pi;

    % α-phase variables
    alpha_vars = num2cell(y(1:9));

    % β-phase variables
    beta_vars = num2cell(y(10:18));

    RHS_pole = @(Q, H, P, r, z, L, s, V, B, S, k, H0, phase) [ ...
        H*L + 0.5*lam - k*H0*H^2 + 0.5*k*H*H0^2;
        0;
        H;
        phase;
        0;
        0;
        1;
        0.75*r*sin(P)*sin(S);
        0.25*k*(2*H - H0)^2 * sin(S);
    ];

    RHS = @(Q, H, P, r, z, L, s, V, B, S, k, H0) [ ...
        (-Q*cos(P)/r - k*(2*H - H0)*(H*H0 + 2*(H - sin(P)/r)^2) + 2*H*L + lam)*sin(S)/r;
        Q/(2*k)*sin(S)/r;
        (2*H - sin(P)/r)*sin(S)/r;
        cos(P)*sin(S)/r;
        sin(P)*sin(S)/r;
        0;
        sin(S)/r;
        0.75*r*sin(P)*sin(S);
        0.25*k*(2*H - H0)^2 * sin(S);
    ];

    if S < delta*pi
        % Taylor Approximaton of the ODEs at s = 0 and pi
        RegionA = RHS_pole(alpha_vars{:}, SA, kA, H0(1), 1);
        RegionB = RHS_pole(beta_vars{:},  SB, kB, H0(2),-1);
    else
        % Bulk ODEs
        RegionA = RHS(alpha_vars{:}, SA, kA, H0(1));
        RegionB = RHS(beta_vars{:},  SB, kB, H0(2));
    end
    dyds = [RegionA*aS; RegionB*bS];
end

function [BCmax, report] = bc_diagnostics(sol, bcfun)
    Ya = sol.y(:,1); Yb = sol.y(:,end); lam = sol.parameters(:);
    res = bcfun(Ya, Yb, lam);
    [BCmax, iMax] = max(abs(res));
    report.idx = iMax;
    report.res = res;

    % --- robust min radius AWAY FROM POLES ---
    s  = sol.x;                    % s in [0, pi]
    rA = abs(sol.y(4,:));          % alpha radius
    rB = abs(sol.y(13,:));         % beta radius

    % choose a small buffer near each pole; purely diagnostic
    hmean = mean(diff(s));
    % heuristics: at least a few mesh spacings, and ~1% of the domain
    buf = max(5*hmean, 0.01*pi);

    mask = (s > buf) & (s < (pi - buf));
    if any(mask)
        rMinAway = min( [ min(rA(mask)), min(rB(mask)) ] );
    else
        % if mesh is too coarse, fall back to whole domain (rare in practice)
        rMinAway = min( [ min(rA), min(rB) ] );
    end

    report.min_r = rMinAway;       % <-- use this in your acceptance gate
end

function idx = nearest_by_sign(cacheKeys, target)
    tgtSign = sign(target);
    same = all(sign(cacheKeys) == tgtSign, 2);
    pool = find(same);
    if isempty(pool), pool = 1:size(cacheKeys,1); end
    [~,j] = min(sum((cacheKeys(pool,:) - target).^2,2));
    idx = pool(j);
end