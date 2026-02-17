function bvp = vesicleBVP_bundle(par)
%VESICLEBVP_BUNDLE  Factory returning odefun/bcfun handles for bvp4c/bvp6c.
%
% Usage:
%   bvp = vesicleBVP_bundle(par);
%   sol = bvp6c(bvp.odefun, bvp.bcfun, initSol, opts);   % or bvp4c
%
% Required par fields (as used by the DE/BC):
%   KA, KB, KG, H0 (1x2), V, aS, bS, delta, poleDeg
%
% Optional par fields:
%   rFloor (default 1e-12), debugInterface (default false)

    par = vesicleBVP_defaults(par);

    % Solver-facing handles (capture par in closure)
    bvp.par   = par;
    bvp.odefun = @(t,y,P) BendV_Lag_EIGp_DE_impl(t,y,P,par);
    bvp.bcfun  = @(ya,yb,P) BendV_Lag_EIGp_BC_impl(ya,yb,P,par);

    % Optional: expose shared helpers for diagnostics/tests
    bvp.interface = @(YA,YB) interface_conditions(YA,YB,par);
    bvp.unpack9   = @unpack9;
end

% =========================
% Defaults / validation
% =========================
function par = vesicleBVP_defaults(par)
    if ~isfield(par,'rFloor') || isempty(par.rFloor)
        par.rFloor = 1e-12;
    end
    if ~isfield(par,'poleDeg') || isempty(par.poleDeg)
        par.poleDeg = 2;
    end
    if ~isfield(par,'debugInterface') || isempty(par.debugInterface), par.debugInterface = false; end
end

% ============================================================
% LOCAL HELPERS.
% ============================================================

function res = BendV_Lag_EIGp_BC_impl(y_poles, y_neck, P, par)
    % BENDV_LAG_EIGP_BC_IMPL  Boundary conditions for the two-phase 
    % vesicle BVP.
    %
    % Interpretation of inputs (consistent with the t-parameterization):
    %   y_poles = Y(t=0)  contains:
    %       alpha at    S_alpha(0) = 0      (south pole) 
    %       beta at     S_beta(0)  = pi     (north pole)
    %   y_neck  = Y(t=pi) contains:
    %       alpha at    S_alpha(pi) = S*    (neck) 
    %       beta at     S_beta(pi)  = S*    (neck)
    %
    % State blocks:
    %   Y_alpha = [Q; H; psi; r; z; L; s; V; E]
    %   Y_beta  = [Q; H; psi; r; z; L; s; V; E]
    %
    % BC count: 19 = 18 states + 1 unknown parameter P.

    % -------- Parameters --------
    Vf = par.V;         % prescribed reduced volume (nondimensional)

    % ---- unpack at poles (t=0) ----
    YA0 = y_poles(1:9);
    YB0 = y_poles(10:18);

    [QA0, HA0, psiA0, rA0, zA0, LA0, sA0, VA0, EA0] = unpack9(YA0); %#ok<ASGLU>
    [QB0, HB0, psiB0, rB0, zB0, LB0, sB0, VB0, EB0] = unpack9(YB0); %#ok<ASGLU>

    % ---- unpack at neck/interface (t=pi) ----
    YAk = y_neck(1:9);
    YBk = y_neck(10:18);

    [QAk, HAk, psiAk, rAk, zAk, LAk, sAk, VAk, EAk] = unpack9(YAk); %#ok<ASGLU>
    [QBk, HBk, psiBk, rBk, zBk, LBk, sBk, VBk, EBk] = unpack9(YBk); %#ok<ASGLU>

    % ============================================================
    % (1) South pole (alpha) regularity + anchoring integrals
    % ============================================================
    % South pole (alpha): smooth axisymmetry requires r=0 and psi=0,
    % and regularity of the force variable requires Q=0.
    % Anchor z to remove vertical translation.
    % Anchor s, V, E to fix integration constants.
    res_south = [
        QA0          % force regularity
        psiA0        % psi=0 at south pole
        rA0          % r=0 at axis
        zA0          % anchor vertical translation
        sA0          % anchor arclength accumulator
        VA0          % anchor volume accumulator
        EA0          % anchor energy accumulator
    ];

    % ============================================================
    % (2) North pole (beta) regularity + anchoring integrals
    % ============================================================
    % We anchor VB0 to prevent the volume constraint from being satisfied
    % by an arbitrary integration constant in V_beta.
    % North pole (beta): smooth axisymmetry requires r=0 and psi=pi,
    % and regularity requires Q=0.
    % Anchor s and V to fix integration constants (important for the 
    % volume constraint).
    %
    % NOTE: We intentionally do NOT constrain EB0 here to keep the 
    % BC count at 19.
    res_north = [
        QB0          % force regularity
        psiB0 - pi   % psi=pi at north pole
        rB0          % r=0 at axis
        sB0          % anchor arclength accumulator
        VB0          % anchor volume accumulator (important)
    ];

    % =====================================================================
    % (3) Interface continuity + matching conditions via shared helper
    % =====================================================================
    ic = interface_conditions(YAk, YBk, par);

    % Volume constraint: with VB0 fixed, this enforces the 
    % shape (not a free offset)
    res_volume = (VAk - VBk) - Vf;

    res_neck = [
        ic.cont
        res_volume
        ic.match
    ];

    % =====================================================================
    % Collect all residuals (7 + 5 + 7 = 19 total)
    % =====================================================================
    res = [res_south; res_north; res_neck];

end

function dyds = BendV_Lag_EIGp_DE_impl(S, y, P, par)
    % BENDV_LAG_EIGP_DE_IMPL  ODE system for two-phase vesicle equilibrium 
    % on t in [0,pi].
    % 
    % Evaluates d/ds of 18-component state vector (9 α-phase, 9 β-phase).
    %
    % Uses pole Taylor truncations for:
    %   sin(psi)/r, sin(S)/r, and the product Q*cos(psi)/r
    % near t=0 (both poles), and bulk formulas elsewhere.
    %
    % phases: A = α, B = β

    % t is the solver parameter in [0,pi]
    t = S;

    % -------- Parameters --------
    kA    = par.KA;
    kB    = par.KB;
    H0    = par.H0;       % [H0_alpha, H0_beta]
    aS    = par.aS;       % S_alpha(t) = aS*t
    bS    = par.bS;       % S_beta(t)  = bS*t + pi  (bS < 0)
    delta = par.delta;    % pole-buffer width in units of pi
    deg   = par.poleDeg;  % (recommended) add this: 1,2,or 3. Default 2.

    if isempty(deg); deg = 2; end

    deltaS = delta * pi;  % interpret delta as an S-buffer fraction of pi

    % Mapped S-variables on each phase
    SA = aS*t;         % in [0, S*]
    SB = bS*t + pi;    % in [pi, S*] decreasing
    
    % Pole distances in S-variable
    epsA = SA;         % south pole ε
    epsB = pi - SB;    % north pole ε (positive; since bS<0, epsB = -bS*t)
    
    % Decide whether we are in the pole buffer (t near 0)
    inPoleA = (epsA < deltaS);
    inPoleB = (epsB < deltaS);

    % Unpack states (A = α, B = β)
    YA = num2cell(y(1:9));    % [Q,H,psi,r,z,L,s,V,E]
    YB = num2cell(y(10:18));

    % Pole signs: +1 at south, -1 at north
    if inPoleA
        RegionA = rhs_pole(YA{:}, SA, kA, H0(1), P, +1, deg, epsA);
    else
        RegionA = rhs_bulk(YA{:}, SA, kA, H0(1), P);
    end
    if inPoleB
        RegionB = rhs_pole(YB{:}, SB, kB, H0(2), P, -1, deg, epsB);
    else
        RegionB = rhs_bulk(YB{:}, SB, kB, H0(2), P);
    end

    % Chain rule for the affine maps S_alpha(t), S_beta(t)
    dyds = [RegionA * aS; RegionB * bS];

end

function pole = pole_asymptotics(eps, Hp, k, H0, L, P, poleSign, deg)
    % POLE_ASYMPTOTICS
    % Pole-safe approximations for singular ratios and products.
    % eps      : distance to pole in S (south: eps=S; north: eps=pi-S)
    % Hp       : H at the pole (use current H in pole region; 
    %            equals H_p asymptotically)
    % k,H0,L,P : phase parameters/states (nondimensional)
    % poleSign : +1 south, -1 north (affects cos(psi)/r expansion sign)
    % deg      : 1,2 (Taylor degree in eps)
    %
    % Returns:
    %   sinpsi_over_r  ~ sin(psi)/r
    %   sinS_over_r    ~ sin(S)/r
    %   Qcos_over_r    ~ Q * cos(psi)/r  (finite limit used in Q equation)
    %   Hpp            = H_p'' at the pole
    %   Qp             = Q'(0) at the pole (in S-variable)

    % Regularity-based Q'(0) at the pole (finite)
    Qp = Hp*L + 0.5*P - k*H0*Hp^2 + 0.5*k*Hp*H0^2;

    % Corresponding H'' at the pole
    Hpp = Qp/(2*k);

    % --- sin(psi)/r ---
    sinpsi_over_r = Hp;
    if deg >= 2
        sinpsi_over_r = sinpsi_over_r + 0.25*Hpp*eps^2;
    end

    % --- sin(S)/r ---
    sinS_over_r = 1;
    if deg >= 2
        sinS_over_r = sinS_over_r - (1 - Hp^2)*eps^2/8;
    end

    % --- Q*cos(psi)/r ---
    % Use a pole-safe evaluation: Q ~ Qp*eps, cos(psi)/r ~ poleSign*(1/eps + c1*eps + c3*eps^3 + ...)
    % => Q*cos/r ~ poleSign*(Qp + Qp*c1*eps^2 + Qp*c3*eps^4 + ...)
    Qcos_over_r = Qp;  % degree 1

    if deg >= 2
        c1 = (1 - 9*Hp^2)/24;      % coefficient of eps in cos/r
        Qcos_over_r = Qcos_over_r + poleSign*(Qp*c1*eps^2);  % contributes at eps^2
    end

    pole.sinpsi_over_r = sinpsi_over_r;
    pole.sinS_over_r   = sinS_over_r;
    pole.Qcos_over_r   = Qcos_over_r;
    pole.Hpp           = Hpp;
    pole.Qp            = Qp;
end

function dy = rhs_bulk(Q, H, psi, r, z, L, s, V, E, S, k, H0, P)
    % Bulk RHS
    sinpsi_over_r = sin(psi)/r;
    sinS_over_r   = sin(S)/r;

    dy = [ ...
        (-Q*cos(psi)/r ...
         - k*(2*H - H0)*(H*H0 + 2*(H - sinpsi_over_r)^2) ...
         + 2*H*L + P ...
        )*sinS_over_r;                          % dQ/dS
        0.5*Q/k * sinS_over_r;                  % dH/dS
        (2*H - sinpsi_over_r) * sinS_over_r;    % dpsi/dS
        cos(psi) * sinS_over_r;                 % dr/dS
        sin(psi) * sinS_over_r;                 % dz/dS
        0;                                      % dL/dS
        sinS_over_r;                            % ds/dS
        0.75*r*sin(psi)*sin(S);                 % dV/dS  (pole-safe)
        0.5*k*(2*H - H0)^2 * sin(S);            % dE/dS  (pole-safe)
    ];
end

function dy = rhs_pole( ...
    Q, H, psi, r, z, L, s, V, E, S, k, H0, P, poleSign, deg, eps)
    % Pole RHS derived from Taylor truncations in eps 
    % (eps = S near south, eps = pi-S near north)

    pole = pole_asymptotics(eps, H, k, H0, L, P, poleSign, deg);

    sinpsi_over_r = pole.sinpsi_over_r;
    sinS_over_r   = pole.sinS_over_r;

    % Replace the singular product Q*cos(psi)/r by its finite asymptotic value
    Qp   = pole.Qp;

    dy = [ ...
        Qp;                                    % dQ/dS
        0.5*Q/k * sinS_over_r;                 % dH/dS (Q itself is O(eps), safe)
        (2*H - sinpsi_over_r) * sinS_over_r;   % dpsi/dS
        cos(psi) * sinS_over_r;                % dr/dS  (finite: sinS_over_r ~ 1)
        sin(psi) * sinS_over_r;                % dz/dS
        0;                                     % dL/dS
        sinS_over_r;                           % ds/dS
        0.75*r*sin(psi)*sin(S);                % dV/dS  (no 1/r)
        0.5*k*(2*H - H0)^2 * sin(S);           % dE/dS  (no 1/r)
    ];
end

function out = interface_conditions(YAk, YBk, par)
%INTERFACE_CONDITIONS  Continuity + matching residuals at the phase interface.
%
% Inputs:
%   YAk, YBk : 9x1 state vectors at the interface for alpha and beta phases
%              [Q; H; psi; r; z; L; s; V; E]
%   par      : struct with fields KA, KB, KG, H0 (H0 is length-2)
%
% Outputs (struct out):
%   out.cont   : 3x1 continuity residuals [psiB-psiA; rB-rA; zB-zA]
%   out.match  : 3x1 matching residuals  [match1; match2; match3]
%   out.ratios : struct with safe ratios:
%                ratios.sinpsi_over_r, ratios.cospsi_over_r, ratios.r_safe
%
% Notes:
%   - Uses alpha-side geometry (psiA, rA) for ratios; continuity enforces equivalence.
%   - Includes a small r-floor to prevent NaNs for poor initial guesses.

    % -------- Parameters --------
    kA = par.KA;        % alpha-phase bending rigidity (nondimensional)
    kB = par.KB;        % beta-phase  bending rigidity (nondimensional)
    kG = par.KG;        % Gaussian rigidity jump parameter (nondimensional)
    H0 = par.H0(:);     % [H0_alpha; H0_beta]

    % optional safety floor
    if isfield(par,'rFloor') && ~isempty(par.rFloor)
        rFloor = par.rFloor;
    else
        rFloor = 1e-12;
    end

    % unpack (alpha neck)
    [QA, HA, psiA, rA, zA, LA, sA, VA, EA] = unpack9(YAk); %#ok<ASGLU>

    % unpack (beta neck)
    [QB, HB, psiB, rB, zB, LB, sB, VB, EB] = unpack9(YBk); %#ok<ASGLU>

    % continuity residuals
    cont = [
        psiB - psiA
        rB   - rA
        zB   - zA
    ];

    % safe ratios (use alpha-side)
    r_safe = rA;
    if abs(r_safe) < rFloor
        r_safe = sign(r_safe + (r_safe==0))*rFloor; % preserve sign if any
    end

    sinpsi_over_r = sin(psiA) / r_safe;
    cospsi_over_r = cos(psiA) / r_safe;

    % matching residuals (nondimensional form consistent with your BC)
    match1 = (QB - QA) - sinpsi_over_r;

    match2 = kB*(2*HB - H0(2)) ...
           - kA*(2*HA - H0(1)) ...
           + kG*sinpsi_over_r;

    match3 = ...
        kA*(2*HA - H0(1))*(HA - sinpsi_over_r + 0.5*H0(1)) ...
      - kB*(2*HB - H0(2))*(HB - sinpsi_over_r + 0.5*H0(2)) ...
      - cospsi_over_r ...
      + (LB - LA);

    out.cont = cont;
    out.match = [match1; match2; match3];

    out.ratios.sinpsi_over_r = sinpsi_over_r;
    out.ratios.cospsi_over_r = cospsi_over_r;
    out.ratios.r_safe        = r_safe;
end

function [Q, H, psi, r, z, L, s, V, E] = unpack9(Y)
    Q   = Y(1);  H   = Y(2);  psi = Y(3);
    r   = Y(4);  z   = Y(5);  L   = Y(6);
    s   = Y(7);  V   = Y(8);  E   = Y(9);
end
