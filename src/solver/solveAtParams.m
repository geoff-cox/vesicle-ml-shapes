% === FILE: solveAtParams.m ===
function [result, meta] = solveAtParams(params, sim, warm)
% SOLVEATPARAMS  Robust continuation solve at target H0 with deterministic fallback.
% Fixes vs uploaded solveAtParams.m.txt:
%   - accepts SP.SaveHomotopy or SP.saveHomotopy
%   - accepts warm.result.sol OR warm.sol (seed)
%   - default continuation origin is [0,0] when no predecessor
%   - axis-cross reseed does not leave H0_from unchanged
%   - homotopy struct is initialized safely
%   - removes calls to undefined helpers try_axis_paths / coarsen_mesh / local_min_radius_interior
%   - uses bc_diagnostics report.min_r for rMin gate

    arguments
        params (1,1) struct
        sim    (1,1) struct
        warm   (1,1) struct = struct()
    end

    TH     = sim.TH;
    delta0 = TH.delta;
    opts0  = TH.opts;

    verbose      = getFlag(sim.SP, {'Verbose','verbose'}, false);
    saveHomotopy = getFlag(sim.SP, {'SaveHomotopy','saveHomotopy'}, false);

    say = @(fmt,varargin) if_verbose(verbose, fmt, varargin{:});

    MP   = sim.MP;
    Htgt = [params.H0_1, params.H0_2];
    [aS,bS] = computePhaseScales(MP.A);

    basePar = struct('A',MP.A,'V',MP.V,'KA',MP.KA,'KB',MP.KB,'KG',MP.KG, ...
                     'aS',aS,'bS',bS,'delta',delta0);

    say('\n=== solveAtParams: H0=[%+.6g,%+.6g] ===', Htgt);

    % ---- initial guess ----
    initSol = [];
    if isfield(warm,'result') && isstruct(warm.result) && isfield(warm.result,'sol')
        initSol = warm.result.sol;
        say('Warm-start from prior solve.');
    elseif isfield(warm,'sol')
        initSol = warm.sol;
        say('Warm-start from seed (warm.sol).');
    end

    if isempty(initSol)
        seedParams = struct('A',MP.A,'V',MP.V,'KG',MP.KG,'KA',MP.KA,'KB',MP.KB);
        initSol = initialGuessFromFile(seedParams, Htgt);
        say('Initial guess from initial-shapes.');
    end

    if ~isfield(initSol,'parameters'); initSol.parameters = []; end

    % ---- continuation origin ----
    if isfield(warm,'fromParams') && all(isfield(warm.fromParams,{'H0_1','H0_2'}))
        H0_from = [warm.fromParams.H0_1, warm.fromParams.H0_2];
    else
        H0_from = [0,0];
    end

    % ---- axis-cross reseed ----
    h1Crosses = (H0_from(1) <= 0 && Htgt(1) >= 0) || (H0_from(1) >= 0 && Htgt(1) <= 0);
    h2Crosses = (H0_from(2) <= 0 && Htgt(2) >= 0) || (H0_from(2) >= 0 && Htgt(2) <= 0);
    crossesZero = (h1Crosses || h2Crosses) && any(H0_from ~= 0);

    if crossesZero
        say('  ↺ axis-cross: reseed at [0,0] and restart continuation.');
        seedParams = struct('A',MP.A,'V',MP.V,'KG',MP.KG,'KA',MP.KA,'KB',MP.KB);
        seedSol = initialGuessFromFile(seedParams, [0,0]);
        if ~isempty(seedSol), initSol = seedSol; end
        H0_from = [0,0];
    end

    stepCaps = [Inf, 0.50, 0.25, 0.12, 0.06, 0.03];

    if isfield(TH,'delta_list') && ~isempty(TH.delta_list)
        deltas = TH.delta_list(:)';
    else
        deltas = [delta0, min(2*delta0, 0.02), max(0.5*delta0, 5e-3)];
    end

    optSets  = {opts0, bvpset(opts0,'RelTol',1e-5,'AbsTol',1e-7)};  % exploration fallback

    homolog = struct('H0',{},'delta',{},'mesh',{},'BCmax',{},'DEmax',{},'rMin',{}); % safe init

    accepted = false;
    sol = [];
    usedPar = basePar;

    minSeg = defaultMinSeg(TH);

    for rung = 1:numel(stepCaps)
        stepMax = stepCaps(rung);

        for dd = 1:numel(deltas)
            deltaNow = deltas(dd);

            for oo = 1:numel(optSets)
                optsNow = optSets{oo};
                curSol  = initSol;

                Hprev   = H0_from;
                Htarget = Htgt;

                if isfinite(stepMax), maxSeg = stepMax;
                else,                maxSeg = norm(Htarget - Hprev); end

                while true
                    v = Htarget - Hprev;
                    L = norm(v);

                    if L <= maxSeg
                        next = Htarget;
                    else
                        next = Hprev + v * (maxSeg / L);
                    end

                    usedPar = basePar;
                    usedPar.H0    = next;
                    usedPar.delta = deltaNow;

                    odefun = @(s_,y_,lam_) BendV_Lag_EIGp_DE_impl(s_,y_,lam_,usedPar);
                    bcfun  = @(ya,yb,lam_) BendV_Lag_EIGp_BC_impl(ya,yb,lam_,usedPar);

                    try
                        curSol = bvp6c(odefun, bcfun, curSol, optsNow);

                        if saveHomotopy
                            [BCi, rep]  = bc_diagnostics(curSol, bcfun);
                            [DEi,~,~]   = de_residual(curSol, odefun);
                            homolog(end+1) = struct('H0',next,'delta',deltaNow,'mesh',numel(curSol.x), ...
                                                    'BCmax',BCi,'DEmax',DEi,'rMin',rep.min_r); %#ok<AGROW>
                        end

                        if all(next == Htarget)
                            break;
                        else
                            Hprev = next;
                        end

                    catch ME
                        prev = maxSeg;
                        maxSeg = max(maxSeg/2, minSeg);
                        say('  step fail: %s | %.4g → %.4g', ME.message, prev, maxSeg);

                        if maxSeg <= minSeg + eps
                            curSol = [];
                            break;
                        end
                    end
                end

                if ~isempty(curSol)
                    [BCmax, rep]	= bc_diagnostics(curSol, bcfun);
                    [DEmax,~,~]		= de_residual(curSol, odefun);
                    rmin			= rep.min_r;
					rNeck			= min(abs(curSol.y(4,end)), abs(curSol.y(13,end)));
					rNeckMin		= defaultArg(TH,'rNeckMin', 0);  % 0 disables

					if (BCmax <= TH.BCmax) && (DEmax <= TH.DEmaxHard) && (rmin >= TH.rMin) && (rNeck >= rNeckMin)
						sol = curSol; accepted = true;
						break
					else
						initSol = curSol; % warm next attempt
					end
                end

                if accepted, break; end
            end
            if accepted, break; end
        end
        if accepted, break; end
    end

    if ~accepted
        error('solveAtParams: failed to converge after all rungs');
    end

    [label, E_total, P_osm] = labelFromSolution(sol);

    bcfunT  = @(ya,yb,lam_) BendV_Lag_EIGp_BC_impl(ya,yb,lam_,usedPar);
    odefunT = @(s_,y_,lam_) BendV_Lag_EIGp_DE_impl(s_,y_,lam_,usedPar);
    [BCmax, rep] = bc_diagnostics(sol, bcfunT);
    [DEmax,~,~]  = de_residual(sol, odefunT);

    result = struct('sol', sol, 'mesh', numel(sol.x));
    meta   = struct('label',label,'E',E_total,'P',P_osm, ...
                    'BCmax',BCmax,'DEmax',DEmax,'mesh',result.mesh, ...
                    'rMinAway', rep.min_r, 'rNeck', rNeck);

    if saveHomotopy
        meta.homotopy = homolog;
    end

    say('=== accepted: mesh=%d | BC=%.2e | DE=%.2e | rMinAway=%.2e | rNeck=%.2e ===\n', ...
        meta.mesh, meta.BCmax, meta.DEmax, meta.rMinAway, meta.rNeck);
end

% ---- helpers ----
function if_verbose(flag, fmt, varargin)
    if flag, fprintf([fmt '\n'], varargin{:}); end
end

function tf = getFlag(SP, names, default)
    tf = default;
    if ~isstruct(SP), return; end
    for k=1:numel(names)
        nm = names{k};
        if isfield(SP,nm), tf = logical(SP.(nm)); return; end
    end
end

function m = defaultMinSeg(TH)
    if isfield(TH,'minH0Step') && TH.minH0Step > 0
        m = TH.minH0Step;
    else
        m = 0.01;
    end
end

function [BCmax, report] = bc_diagnostics(sol, bcfun)

    lam = [];
    if isfield(sol,'parameters') && ~isempty(sol.parameters)
        lam = sol.parameters(:);
    end
    Ya = sol.y(:,1); Yb = sol.y(:,end);
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

function [ok, initSol, fromH0] = try_axis_paths(warmSol, Hfrom, Htgt, sim, stepMax)

    % TRY_AXIS_PATHS  Axis-aligned fallback when 2D continuation stalls.
    % Try two one-dimensional homotopies: (H0_1 then H0_2) and (H0_2 then H0_1).
    % Returns first that completes at least one accepted step (ok=true) along either axis.

    ok = false; initSol = []; fromH0 = Hfrom;

    % Path A: H0_1 -> H0_2
    [okA, solA] = oned_homotopy(warmSol, Hfrom, [Htgt(1) Hfrom(2)], sim, stepMax);
    if okA
        [okB, solB] = oned_homotopy(solA, [Htgt(1) Hfrom(2)], Htgt, sim, stepMax);
        if okB, ok = true; initSol = solB; fromH0 = Htgt; return; end
    end

    % Path B: H0_2 -> H0_1
    [okB1, solB1] = oned_homotopy(warmSol, Hfrom, [Hfrom(1) Htgt(2)], sim, stepMax);
    if okB1
        [okB2, solB2] = oned_homotopy(solB1, [Hfrom(1) Htgt(2)], Htgt, sim, stepMax);
        if okB2, ok = true; initSol = solB2; fromH0 = Htgt; return; end
    end
end

function [ok, solOut] = oned_homotopy(solIn, p0, p1, sim, stepMax)
    % March in fixed-size steps along one axis; accept first step that passes gates.
    ok = false; solOut = solIn;
    if all(p0==p1), ok = true; return; end
    nSteps = max(1, ceil(norm(p1-p0)/stepMax));
    for k=1:nSteps
        alpha = k/nSteps;
        Hk = (1-alpha)*p0 + alpha*p1;
        [accept, solK] = try_one_shot(solOut, Hk, sim);
        if accept
            ok = true; solOut = solK; return;  % return as soon as we get one accepted step
        else
            % If Newton failed hard, coarsen the mesh and retry once at same Hk
            [accept2, solK2] = try_one_shot(coarsen_mesh(solOut, 0.5), Hk, sim);
            if accept2
                ok = true; solOut = solK2; return;
            end
        end
    end
end

function [accept, sol] = try_one_shot(solInit, H0, sim)
    % Single attempt at a given H0 using the current ladders deltas/tols.
    accept = false; sol = [];
    TH = sim.TH; MP = sim.MP;
    [aS,bS] = computePhaseScales(MP.A);
    Par0 = struct('H0',H0,'A',MP.A,'V',MP.V,'KA',MP.KA,'KB',MP.KB,'KG',MP.KG, ...
                  'aS',aS,'bS',bS,'delta',sim.TH.delta);

    deltas = sim.TH.delta_list;  % e.g., [0.01, 0.02, 0.005]
    opts   = {sim.TH.opts, bvpset(sim.TH.opts,'RelTol',1e-5,'AbsTol',1e-7)};

    for d = 1:numel(deltas)
      for o = 1:numel(opts)
          Par = Par0; Par.delta = deltas(d);
          odefun = @(s,y,lam) BendV_Lag_EIGp_DE_impl(s,y,lam,Par);
          bcfun  = @(ya,yb,lam) BendV_Lag_EIGp_BC_impl(ya,yb,lam,Par);
          try
              sol1 = bvp6c(odefun, bcfun, solInit, opts{o});
              [BCmax,~] = bc_diagnostics(sol1, bcfun);
              [DEmax,~,~] = de_residual(sol1, odefun);
              rmin = local_min_radius_interior(sol1);
              if (BCmax <= TH.BCmax) && (DEmax <= TH.DEmaxHard) && (rmin >= TH.rMin)
                  accept = true; sol = sol1; return;
              else
                  solInit = sol1; % warm next try
              end
          catch
              % continue
          end
      end
    end
end

function sol2 = coarsen_mesh(sol1, ratio)
    % Simple resample to a coarser grid to regularize Jacobian.
    if nargin<2, ratio = 0.5; end
    x = sol1.x;  m = max(6, ceil(numel(x)*ratio));
    x2 = linspace(x(1), x(end), m);
    y2 = deval(sol1, x2);
    sol2 = bvpinit(x2, y2, sol1.parameters);
end

function [rMax, rComp, worstIdx] = de_residual(sol, odefun)
    % DE_RESIDUAL  Max residual of ODE on a nonuniform mesh (central difference).
    % Uses second-order nonuniform central differences on interior nodes only.

    lam = [];
    if isfield(sol,'parameters') && ~isempty(sol.parameters)
        lam = sol.parameters(:);
    end
    s   = sol.x;             % s in [0, pi]
    Y   = sol.y;
    n   = numel(s);
    if n < 3, error('de_residual: need at least 3 mesh points'); end

    % interior nodes (2..n-1) using nonuniform central difference weights
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
    % We will be conservative and drop 2*mean(h) near each pole as a buffer:
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

function initSol = initialGuessFromFile(params, H0)

    here    = fileparts(mfilename('fullpath'));   % src/solver
    srcRoot = fileparts(here);                     % src
    ishapes = fullfile(srcRoot,'initial-shapes');

    initSol = [];

    if isfield(params,'initialGuess') && ~isempty(params.initialGuess)
        cand = params.initialGuess;
        if exist(cand,'file')~=2
            cand = fullfile(ishapes, cand);
        end
        assert(exist(cand,'file')==2, 'initialGuessFromFile: file not found: %s', cand);
        tmp = load(cand);
        initSol = tmp.Version(1).Solution;
        if ~isfield(initSol,'parameters'), initSol.parameters = []; end
        return
    end

    d = dir(fullfile(ishapes, 'SIM_Node_*.mat'));
    assert(~isempty(d), 'initialGuessFromFile: no SIM_Node_*.mat found in %s', ishapes);

    % If physics is provided, try exact seed file first:
    if all(isfield(params,{'A','V','KG','KA','KB'}))
        base = sprintf('SIM_Node_%d_%d_%d_%d_%d_+00_+00.mat', ...
            round(params.A*100), round(params.V*100), params.KG, params.KA, params.KB);
        cand = fullfile(ishapes, base);
        if exist(cand,'file')==2
            tmp = load(cand);
            initSol = tmp.Version(1).Solution;
            if ~isfield(initSol,'parameters'), initSol.parameters = []; end
            return
        end
    end

    % Otherwise choose the first file (safe fallback).
    % NOTE: Removed the incorrect "nearest by first token" heuristic. fileciteturn25file8L47-L56
    f = fullfile(d(1).folder, d(1).name);
    tmp = load(f);
    initSol = tmp.Version(1).Solution;
    if ~isfield(initSol,'parameters'), initSol.parameters = []; end
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

function rmin = local_min_radius_interior(sol)
    % Min radius away from poles (ignore s=0 and s=pi so r=0 at poles is allowed)
    rA = abs(sol.y(4,:));
    rB = abs(sol.y(13,:));
    if numel(rA) >= 3, rA = rA(2:end-1); end
    if numel(rB) >= 3, rB = rB(2:end-1); end
    rmin = min([rA(:); rB(:)]);
    if isempty(rmin), rmin = Inf; end
end

function res = BendV_Lag_EIGp_BC_impl(y_poles, y_neck, lam, par)
    % BENDV_LAG_EIGP_BC_IMPL  Boundary conditions at poles and neck junction.
    % South pole (α phase): regularity in Q,P,r,z plus integral constraints.
    % North pole (β phase): regularity in Q,P plus integral constraints.
    % Neck junction: continuity of geometry and force balance between phases.

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
        rBn
        % zBn
        % LBn
        sBn
        % VBn
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
    % BENDV_LAG_EIGP_DE_IMPL  ODE system for two-phase vesicle equilibrium.
    % Evaluates d/ds of 18-component state vector (9 α-phase, 9 β-phase).
    % Uses a Taylor-expanded RHS near poles (S < delta*pi), otherwise bulk RHS.

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

    % α-phase variables: [Q,H,P,r,z,L,s,V,E]
    alpha_vars = num2cell(y(1:9));

    % β-phase variables: [Q,H,P,r,z,L,s,V,E]
    beta_vars = num2cell(y(10:18));

    % RHS_pole handles singular pole expansion; RHS is bulk form.
    % FIXED: Pole expansion requires sin(S)/r = 1 (from ds/dS consistency)
    % Previous version incorrectly assumed sin(S)/r = 1/2, causing factor-of-2 errors
    % Corrected: February 2026 - Issue 1 from code audit
    RHS_pole = @(Q, H, P, r, z, L, s, V, B, S, k, H0, phase) [ ...
        2*H*L + lam - 2*k*H0*H^2 + k*H*H0^2;
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
