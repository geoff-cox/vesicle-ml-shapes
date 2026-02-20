clc; clear;
%% --- solver knobs for the re-solves ---
delta = 1e-4;              % use your usual delta (can sweep this too)
opts  = bvpset('RelTol',1e-6,'AbsTol',1e-8,'NMax',1500);
useLegacy = false;
poleDeg = 2;
minStep = 0.0001;
defaultStep = 0.1;

init_sols = load("src\initial-shapes\SIM_Node_50_72_0_1_1_+00_+00.mat");
init_sol = init_sols.Version(1).Solution;
sol = bvpinit(init_sol.x, init_sol.y, init_sol.parameters);

A = 0.5000;
V = 0.7200;
KA = 1;
KB = 1;
KG = 0;
aS = 0.5000;
bS = -0.5000;

startH0 = [0 0];
targetH0 = [0 -10];

direction = targetH0 - startH0;
dH0 = direction/norm(direction);

Par = struct( ...
    'A',A,'V',V,'KA',KA,'KB',KB,'KG',KG, ...
    'aS',aS,'bS',bS,'delta',delta,'H0',startH0, ...
    'poleDeg',poleDeg);

step = defaultStep;
checkpoint = [0,-1];
chkptSol = sol;
chkptH0 = Par.H0;
lastZ = sol.y(14,1);
while step > minStep
    bvp = vesicleBVP_bundle(Par);
    try
        sol = bvp6c(bvp.odefun, bvp.bcfun, sol, opts);
    catch
        sol = chkptSol;
        Par.H0 = chkptH0;
        step = step*0.2;
        checkpoint = chkptH0;
    end

    thisZ = sol.y(14,1);
    zJump = abs(thisZ - lastZ);
    fprintf("z-jump = %g -- ",zJump)
    lastZ = thisZ;

    if zJump > 0.2
        sol = chkptSol;
        Par.H0 = chkptH0;
        step = step*0.5;
        checkpoint = chkptH0;
    elseif zJump < 0.02
        step = step*2;
    end

    if norm(Par.H0 - startH0) > norm(checkpoint - startH0)
        chkptSol = sol;
        checkpoint = checkpoint + 1;
        chkptH0 = Par.H0;
    end

    rA = sol.y(4,:);  zA = sol.y(5,:);
    rB = sol.y(13,:); zB = sol.y(14,:);

    plot(rA, zA, 'LineWidth', 3); hold on;
    plot(rB, zB, 'LineWidth', 3); hold off;
    % axis image
    drawnow
    % if sol.fevals > 150000
    %     step = step*0.5;
    % elseif sol.fevals < 100000
    %     step = min(step*2,defaultStep);
    % end
    Par.H0 = Par.H0 + step*dH0;
    fprintf("H0 = [%g, %g]\n",Par.H0)
    if norm(Par.H0 - startH0) > norm(direction)
        break
    end
end

startH0 = [0 -4];
targetH0 = [4 -4];

direction = targetH0 - startH0;
dH0 = direction/norm(direction);

Par = struct( ...
    'A',A,'V',V,'KA',KA,'KB',KB,'KG',KG, ...
    'aS',aS,'bS',bS,'delta',delta,'H0',startH0, ...
    'poleDeg',poleDeg);

step = defaultStep;
while step > minStep
    bvp = vesicleBVP_bundle(Par);
    try
        sol = bvp6c(bvp.odefun, bvp.bcfun, sol, opts);
        step = defaultStep;
    catch
        step = step*0.1;
    end
    rA = sol.y(4,:);  zA = sol.y(5,:);
    rB = sol.y(13,:); zB = sol.y(14,:);
    pause(0.2)

    plot(rA, zA, 'LineWidth', 3); hold on;
    plot(rB, zB, 'LineWidth', 3); hold off;
    axis image
    drawnow
    Par.H0 = Par.H0 + defaultStep*dH0;
    fprintf("H0 = [%g, %g]\n",Par.H0)
    if norm(Par.H0 - startH0) > norm(direction)
        break
    end
end
% sol = bvp6c(odefun, bcfun, guess, opts)
% 
% rA = sol.y(4,:);  zA = sol.y(5,:); hold on;
% rB = sol.y(13,:); zB = sol.y(14,:);
% 
% plot(rA, zA, 'LineWidth', 3);
% plot(rB, zB, 'LineWidth', 3);
% axis image

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
        % VAs
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
        % 2*H*L + lam - 2*k*H0*H^2 + k*H*H0^2;
        H*L  + 0.5*lam - k*H0*H^2  + 0.5*k*H*H0^2;
        0;
        H;
        phase;
        0;
        0;
        1;
        0.75*r*sin(P)*sin(S);
        0.5*k*(2*H - H0)^2 * sin(S);
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
        0.5*k*(2*H - H0)^2 * sin(S);
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