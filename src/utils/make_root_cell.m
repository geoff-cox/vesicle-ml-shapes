function C = make_root_cell(bounds)
    % bounds = [H0_1_min H0_1_max; H0_2_min H0_2_max]
    H1 = bounds(1,:); H2 = bounds(2,:);
    cLL = [H1(1) H2(1)]; cLR = [H1(2) H2(1)];
    cUL = [H1(1) H2(2)]; cUR = [H1(2) H2(2)];
    C = struct();
    C.corners        = [cLL; cLR; cUL; cUR];
    C.depth          = 0;
    C.cornerSolved   = false(4,1);
    C.cornerLabel    = strings(4,1);
    C.cornerEnergy   = nan(4,1);
    C.cornerPressure = nan(4,1);
    C.isUniform      = false;
    C.mixedEdges     = zeros(0,2);
end