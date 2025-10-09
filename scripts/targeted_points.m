% --- targeted_points.m
run('bootstrap.m');

% Small menu of (H0_1,H0_2) tests:
params = [ ...
    0, 0;     % baseline (compare to Fig. 3(b) families)
    1, 0;    -1, 0; 
    0, 1;     0,-1;
    1, 1;    -1,-1;
    2, -2;   -2, 2   % opposite-sign pair
];

% Helper to call the core solver at a point using solveAtParams:
% (We reuse sim_driver's plumbing by creating a "sim" and "cache" quickly.)
sim = makeSimConfig();                  % from your driver file (on path)
sim.debug.short = true;                 % keep meshes small for quick checks
cache = bootstrapCache(sim);            % seeds the cache

for k = 1:size(params,1)
    sim.H0 = params(k,:);
    try
        [sol, meta] = solveAtParams(sim, cache);   %#ok<NASGU>
        cache = meta.cache;
        fprintf('OK @ (%.2f, %.2f) | label=%d | E=%.4g | P=%.4g\n', ...
                sim.H0(1), sim.H0(2), meta.label, meta.E_total, meta.P_osm);
    catch ME
        fprintf('FAIL @ (%.2f, %.2f): %s\n', sim.H0(1), sim.H0(2), ME.message);
    end
end
