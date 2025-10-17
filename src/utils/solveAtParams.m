% -------------------------------------------------------------------------
% EXTRACTED HELPER for "solveAtParams"
%   - Source: sim_driver_quad_tree_full.m
%   - Extracted: 2025-10-11 11:50:11
%   - Sub-helpers (nested functions) are retained in this file.
% -------------------------------------------------------------------------

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