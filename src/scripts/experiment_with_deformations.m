% test_multistart_deform_at_fixed_H0.m
%
% Minimal "multi-start via deformation" test harness:
%   (1) pick an (H0_1,H0_2), load an existing solution from the catalog
%   (2) deform the solution (various strategies)
%   (3) re-solve at the SAME H0 using bvp6c (or a solveAtParams fallback)
%   (4) plot original vs deformed guess vs converged solution
%
% Run from anywhere (it finds src/ and sim-results/ relative to this file).

%% --- project paths ---
here    = fileparts(mfilename('fullpath'));     % put this file anywhere inside repo
% Heuristic: walk upward until we find "src"
p = here;
for k = 1:6
    if exist(fullfile(p,'src'),'dir')
        projRoot = p;
        break
    end
    p = fileparts(p);
end
assert(exist('projRoot','var')==1, ...
    'Could not locate repo root containing src/. Move this script inside the repo.');

srcRoot = fullfile(projRoot,'src');
addpath(genpath(srcRoot));

simDir  = fullfile(projRoot,'sim-results');
resDir  = fullfile(simDir,'hashed_results');

%% --- choose target H0 and (optionally) which branch tag to load ---
targetH0   = [0, 0];   % <-- CHANGE ME
branchTag  = "";              % "" = single-branch results; or e.g. "NE", "SW", etc.
H0tol      = 1e-10;           % exact-match tolerance (fallback is nearest)

%% --- solver knobs for the re-solves ---
TH = struct();
TH.delta = 0.01;              % use your usual delta (can sweep this too)
TH.opts  = bvpset('RelTol',1e-6,'AbsTol',1e-8,'NMax',1500);

% Optional "sanity" gates (loose; this is a testing script)
gateBCmax = 1e-4;             % boundary residual max norm threshold
gateMinR  = 0;                % set >0 if you want to reject self-intersections

%% --- load solution at (or near) targetH0 via the catalog ---
[hash0, entry0] = find_hash_for_H0(simDir, targetH0, branchTag, H0tol);
S0 = load(fullfile(resDir, hash0 + ".mat"), 'result', 'meta');
sol0 = S0.result.sol;

fprintf('Loaded hash %s for H0=[%+.6g,%+.6g]\n', hash0, targetH0(1), targetH0(2));
fprintf('  E=%.6g, P=%.6g, label=%s\n', ...
    defaultArg(S0.meta,'E',NaN), defaultArg(S0.meta,'P',NaN), string(defaultArg(S0.meta,'label',"")));

%% --- build Par struct for bvp6c (physics + H0 + delta) from catalog entry ---
P = entry0.params;  % contains A,V,KA,KB,KG,H0_1,H0_2, aS,bS (usually)
if ~isfield(P,'aS') || ~isfield(P,'bS')
    [P.aS, P.bS] = computePhaseScales(P.A);
end
Par = struct('A',P.A,'V',P.V,'KA',P.KA,'KB',P.KB,'KG',P.KG, ...
             'aS',P.aS,'bS',P.bS,'delta',TH.delta,'H0',targetH0);

% Decide whether we can call the model ODE/BC directly (preferred)
haveModelFns = (exist('BendV_Lag_EIGp_DE_impl','file')==2) && (exist('BendV_Lag_EIGp_BC_impl','file')==2);
if haveModelFns
    odefun = @(s,y,lam) BendV_Lag_EIGp_DE_impl(s,y,lam,Par);
    bcfun  = @(ya,yb,lam) BendV_Lag_EIGp_BC_impl(ya,yb,lam,Par);
else
    warning(['Could not find BendV_Lag_EIGp_DE_impl / _BC_impl as standalone functions.\n' ...
             'Falling back to solveAtParams as the re-solver (still uses your deformed guess).']);
end

%% --- deformation strategies (add/remove freely) ---
strategies = {};

% (A) Small sinusoidal perturbation to P(s) (both phases), vanishing at poles
strategies{end+1} = struct( ...
    'name',"P-mode sin(3s), eps sweep", ...
    'maker',@(sol,eps) deform_P_mode(sol, eps, 3, "sin_integer"));

% (B) Neck-biased perturbation (nonzero at s=pi), still zero at s=0
strategies{end+1} = struct( ...
    'name',"P neck-kick sin((2.5)s), eps sweep", ...
    'maker',@(sol,eps) deform_P_mode(sol, eps, 2, "sin_half"));

% (C) Random smooth combination of a few sine modes (deterministic seed)
strategies{end+1} = struct( ...
    'name',"P random smooth (modes 2..5)", ...
    'maker',@(sol,eps) deform_P_random(sol, eps, 2:5, 12345));

% (D) Pressure / lambda jitter only (often surprisingly effective)
strategies{end+1} = struct( ...
    'name',"lambda jitter", ...
    'maker',@(sol,eps) deform_lambda(sol, eps));

% amplitudes to try (bigger first; we'll back off on failure)
epsList = [1e-1, 3e-2, 1e-2, 3e-3, 1e-3];

%% --- run each strategy and plot ---
t = tiledlayout(numel(strategies), 1, 'TileSpacing','compact', 'Padding','compact');
title(t, sprintf('Fixed-H0 multi-start test | H0=[%+.3g,%+.3g] | base hash=%s', ...
      targetH0(1), targetH0(2), hash0), 'Interpreter','none');

for i = 1:numel(strategies)
    st = strategies{i};
    best = struct('ok',false,'sol',[],'guess',[],'eps',NaN,'msg',"");

    for eps = epsList
        guess = st.maker(sol0, eps);

        try
            if haveModelFns
                solNew = bvp6c(odefun, bcfun, guess, TH.opts);
            else
                % Fallback: use solveAtParams at fixed H0 with your deformed guess as warm start
                sim = struct('MP',struct('A',P.A,'V',P.V,'KA',P.KA,'KB',P.KB,'KG',P.KG), ...
                             'TH',struct('delta',TH.delta,'opts',TH.opts,'BCmax',inf,'DEmaxHard',inf,'rMin',-inf,'minH0Step',0.01), ...
                             'SP',struct('Verbose',false,'saveHomotopy',false));
                params = struct('H0_1',targetH0(1),'H0_2',targetH0(2));
                warm   = struct('sol',guess,'fromParams',params);
                out    = solveAtParams(params, sim, warm);
                solNew = out.sol;
            end

            % quick diagnostics
            [BCmax, minR] = bc_and_minr(solNew, haveModelFns, bcfun);
            if BCmax <= gateBCmax && minR >= gateMinR
                best.ok   = true;
                best.sol  = solNew;
                best.guess= guess;
                best.eps  = eps;
                best.msg  = "converged";
                break
            else
                % keep the best "near miss" if you want; for now just report and continue
                best.msg = sprintf("rejected by gates (BCmax=%.2e, minR=%.2e)", BCmax, minR);
            end

        catch ME
            best.msg = string(ME.message);
        end
    end

    % --- plotting ---
    nexttile; hold on; axis equal; grid on;
    plot_shape(sol0, '-', 1.6, 'Original');

    if ~isempty(best.guess)
        plot_shape(best.guess, '--', 1.0, sprintf('Deformed guess (eps=%.1e)', best.eps));
    end
    if best.ok
        plot_shape(best.sol, '-', 1.6, sprintf('New solution (eps=%.1e)', best.eps));
    end

    line([0 0], ylim, 'Color',[0.5 0.5 0.5], 'LineStyle','--');
    xlabel('r'); ylabel('z');
    ttl = sprintf('[%d] %s | %s', i, st.name, best.msg);
    title(ttl, 'Interpreter','none');

    % print a short numeric summary
    if best.ok
        d = shape_distance(sol0, best.sol);
        [Enew, Pnew] = energy_and_pressure(best.sol);
        fprintf('[%d] %s: OK at eps=%.1e | dist=%.3e | E=%.6g | P=%.6g\n', ...
            i, st.name, best.eps, d, Enew, Pnew);
    else
        fprintf('[%d] %s: FAIL | last msg: %s\n', i, st.name, best.msg);
    end
end

%% ================= local helpers =================

function guess = deform_P_mode(sol, eps, k, modeType)
% Perturb tangent angle P(s) for both phases.
% alpha P index = 3, beta P index = 12.
    s = sol.x;
    switch modeType
        case "sin_integer"
            bump = sin(k*s);                 % 0 at s=0 and s=pi
        case "sin_half"
            bump = sin((k+0.5)*s);           % 0 at s=0, nonzero at s=pi (neck kick)
        otherwise
            error('Unknown modeType');
    end

    y = sol.y;
    y(3,:)  = y(3,:)  + eps*bump;            % alpha P
    y(12,:) = y(12,:) + eps*bump;            % beta  P

    guess = bvpinit(sol.x, y, safe_params(sol.parameters));
end

function guess = deform_P_random(sol, eps, modes, seed)
% Random smooth combination of sin(m*s), vanishing at poles.
    rng(seed);
    s = sol.x;
    a = randn(size(modes));
    a = a / norm(a);                         % normalize direction
    bump = zeros(size(s));
    for j = 1:numel(modes)
        bump = bump + a(j)*sin(modes(j)*s);
    end
    y = sol.y;
    y(3,:)  = y(3,:)  + eps*bump;
    y(12,:) = y(12,:) + eps*bump;

    guess = bvpinit(sol.x, y, safe_params(sol.parameters));
end

function guess = deform_lambda(sol, eps)
% Jitter the free parameter lambda (osmotic pressure) only.
% eps is relative perturbation scale; e.g. eps=1e-2 means ±1% (here: deterministic sign).
    lam = safe_params(sol.parameters);
    if isempty(lam), lam = 0; end
    lam(1) = lam(1) * (1 + eps);  % you can also try +/- or additive
    guess = bvpinit(sol.x, sol.y, lam);
end

function plot_shape(sol, ls, lw, lbl)
% 2D profile overlay for alpha and beta phases.
    rA = sol.y(4,:);  zA = sol.y(5,:);
    rB = sol.y(13,:); zB = sol.y(14,:);

    plot(rA, zA, ls, 'LineWidth', lw);
    plot(rB, zB, ls, 'LineWidth', lw);

    % annotate lightly (legend is too noisy across many overlays)
    if nargin>=4 && strlength(string(lbl))>0
        text(rA(round(end*0.7)), zA(round(end*0.7)), " " + string(lbl), ...
            'FontSize',9, 'Interpreter','none');
    end
end

function [E_total, P_osm] = energy_and_pressure(sol)
% Mirrors the quick diagnostics used elsewhere in your codebase.
    E_total = sol.y(9,end) - sol.y(18,end) + 0.5*sol.y(4,end);
    if isfield(sol,'parameters') && ~isempty(sol.parameters)
        P_osm = sol.parameters(1);
    else
        P_osm = NaN;
    end
end

function d = shape_distance(solA, solB)
% Simple L2 distance between (r,z) profiles on a common grid.
    sRef = linspace(0, pi, 300);
    YA = deval(solA, sRef);
    YB = deval(solB, sRef);

    % compare both phases' (r,z)
    vA = [YA(4,:)  YA(5,:)  YA(13,:) YA(14,:)];
    vB = [YB(4,:)  YB(5,:)  YB(13,:) YB(14,:)];
    d = sqrt(mean((vA(:)-vB(:)).^2));
end

function [BCmax, minRaway] = bc_and_minr(sol, haveModelFns, bcfun)
% Boundary residual at endpoints + rough min radius away from poles.
    BCmax = NaN;
    if haveModelFns
        lam = safe_params(sol.parameters);
        res = bcfun(sol.y(:,1), sol.y(:,end), lam);
        BCmax = max(abs(res));
    end

    s  = sol.x;
    rA = abs(sol.y(4,:));
    rB = abs(sol.y(13,:));

    if numel(s) >= 10
        buf = max(5*mean(diff(s)), 0.01*pi);
        mask = (s > buf) & (s < (pi-buf));
        if any(mask)
            minRaway = min([min(rA(mask)), min(rB(mask))]);
        else
            minRaway = min([min(rA), min(rB)]);
        end
    else
        minRaway = min([min(rA), min(rB)]);
    end
end

function lam = safe_params(p)
    if nargin==0 || isempty(p)
        lam = [];
    else
        lam = p(:);
    end
end

function val = defaultArg(S, field, fallback)
% Safe meta-field getter (local fallback if you don't have defaultArg.m on path).
    if isstruct(S) && isfield(S,field) && ~isempty(S.(field))
        val = S.(field);
    else
        val = fallback;
    end
end

function [hash, entry] = find_hash_for_H0(simDir, targetH0, branchTag, tol)
% Use catalog to find the hash for a given H0 (exact within tol, else nearest).
    T = catalog_load(simDir);
    assert(istable(T) && any(strcmp(T.Properties.VariableNames,'entry')), ...
        'catalog_load(simDir) must return a table with variable "entry".');

    H = nan(height(T),2);
    ok = false(height(T),1);
    tag = strings(height(T),1);

    for i=1:height(T)
        e = T.entry{i};
        if ~isstruct(e) || ~isfield(e,'params'), continue; end
        p = e.params;
        if isfield(p,'H0_1') && isfield(p,'H0_2')
            H(i,:) = [p.H0_1, p.H0_2];
            ok(i) = true;
        end
        if isfield(e,'meta') && isfield(e.meta,'branch_tag')
            tag(i) = string(e.meta.branch_tag);
        else
            tag(i) = "";
        end
    end

    % keep only actual simulation entries (not seeds)
    idx = find(ok);
    H = H(idx,:);
    tag = tag(idx);
    E = T.entry(idx);

    % optional branch tag filter
    if strlength(string(branchTag))>0
        keep = (tag == string(branchTag));
        idx  = idx(keep);
        H    = H(keep,:);
        E    = E(keep);
        assert(~isempty(idx), 'No catalog entries found with branchTag="%s".', string(branchTag));
    end

    d = vecnorm(H - targetH0(:).', 2, 2);
    [dmin, j] = min(d);

    if dmin <= tol
        entry = E{j};
        hash  = string(entry.meta.hash);
    else
        fprintf('No exact match within tol=%.1e; using nearest H0=[%+.6g,%+.6g] (dist=%.3g)\n', ...
            tol, H(j,1), H(j,2), dmin);
        entry = E{j};
        hash  = string(entry.meta.hash);
    end
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
        - kB*(2*HBk - H0(2))*(HBk - sin(PBk)/rBk + H0(2)/2) ...
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
        %2*H*L + lam - 2*k*H0*H^2 + k*H*H0^2;
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
