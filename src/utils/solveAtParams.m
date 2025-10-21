function [result, meta] = solveAtParams(params, sim, warm)
% SOLVEATPARAMS  Compute a vesicle solution at (H0_1,H0_2) with robust gates.
% Contract (refactor):
%   INPUT:
%     params : struct with fields at least H0_1, H0_2.
%              Optional physics knobs: A, V, KA, KB, KG, etc.
%              Optional solver knobs:  delta, opts (bvpset struct)
%     warm   : optional struct for warm-starts. Recognized fields:
%              - warm.result.sol  or warm.sol   : a prior bvpinit/bvp6c solution struct
%              - warm.meta        : prior meta (ignored here except maybe for seeds)
%   OUTPUT:
%     result : struct with fields {s,y,bcResidual,deResidual} OR the raw bvp6c 'sol'
%              (below we return the raw bvp6c 'sol' as 'result.sol' and fill quick diags)
%     meta   : struct with quick scalars for catalog:
%              label (string), E (double), P (double),
%              BCmax (double), DEmax (double), mesh (double)

    arguments
        params (1,1) struct
        sim    (1,1) struct
        warm   (1,1) struct = struct()
    end

    % --- thresholds & numerical knobs ---
    TH     = sim.TH;
    delta0 = TH.delta;
    opts0  = TH.opts;

    % --- physics bundle (fixed for this run) ---
    MP = sim.MP;
    H0 = [params.H0_1, params.H0_2];
    [aS,bS] = computePhaseScales(MP.A);
    Par = struct('H0',H0, 'A',MP.A,'V',MP.V,'KA',MP.KA,'KB',MP.KB,'KG',MP.KG, ...
                 'aS',aS,'bS',bS, 'delta',delta0);

    % --- initial guess / warm start ---
    initSol = [];
    if isfield(warm,'result') && isstruct(warm.result) && isfield(warm.result,'sol')
        initSol = warm.result.sol;
    elseif isfield(warm,'sol') && isstruct(warm.sol)
        initSol = warm.sol;
    end
    if isempty(initSol)
        % fallback seed from InitialShapes matching physics
        seedParams = struct('A',MP.A,'V',MP.V,'KG',MP.KG,'KA',MP.KA,'KB',MP.KB);
        initSol = initialGuessFromFile(seedParams, H0);
    end

    % --- (optional) continuation from neighbor toward target H0 ---
    if ~isempty(initSol) && isfield(warm,'fromParams') && all(isfield(warm.fromParams, {'H0_1','H0_2'}))
        Hfrom = [warm.fromParams.H0_1, warm.fromParams.H0_2];
        if norm(H0 - Hfrom) > 0.15
            try
                initSol = continuation_towards_H0(initSol, Hfrom, H0, sim, 0.15);
            catch
                % continuation is a best-effort preconditioner; if it fails, keep initSol as-is
            end
        end
    end

    % --- attempt ladder (delta / tolerances) ---
    attempts = [ ...
      struct('delta', delta0,                        'opts', opts0); ...
      struct('delta', max(0.5*delta0, 1e-3),         'opts', opts0); ...
      struct('delta', max(0.25*delta0, 1e-3),        'opts', opts0); ...
      struct('delta', delta0, 'opts', bvpset(opts0,'RelTol',1e-5,'AbsTol',1e-7)) ...
    ];

    solved = false; sol = [];
    for a = 1:numel(attempts)
        Par.delta = attempts(a).delta;          % numerical, not in hash
        Par.H0    = H0;                         % ensure H0 is current target
        odefun = @(s,y,lam) BendV_Lag_EIGp_DE_impl(s,y,lam,Par);
        bcfun  = @(ya,yb,lam) BendV_Lag_EIGp_BC_impl(ya,yb,lam,Par);
        try
            sol = bvp6c(odefun, bcfun, initSol, attempts(a).opts);

            % gates
            [BCmax,~]       = bc_diagnostics(sol, bcfun);
            [DEmax,~,~]     = de_residual(sol, odefun);
            rmin            = local_min_radius_interior(sol);

            if (BCmax <= TH.BCmax) && (DEmax <= TH.DEmaxHard) && (rmin >= TH.rMin)
                solved = true; break;
            else
                initSol = sol;  % warm next rung from current iterate
            end
        catch
            % try next rung
        end
    end
    if ~solved, error('solveAtParams: failed to converge'); end

    % --- label & metrics ---
    [label, E_total, P_osm] = labelFromSolution(sol);
    [BCmax,~] = bc_diagnostics(sol, @(ya,yb,lam) BendV_Lag_EIGp_BC_impl(ya,yb,lam,Par));
    [DEmax,~,~] = de_residual(sol, @(s,y,lam) BendV_Lag_EIGp_DE_impl(s,y,lam,Par));

    result = struct('sol', sol, 'mesh', numel(sol.x));
    meta   = struct('label',label,'E',E_total,'P',P_osm, ...
                    'BCmax',BCmax,'DEmax',DEmax,'mesh',result.mesh);
end
