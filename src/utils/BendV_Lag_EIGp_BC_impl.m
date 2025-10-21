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