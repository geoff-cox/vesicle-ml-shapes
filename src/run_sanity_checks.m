% === FILE: run_sanity_checks.m ===
function run_sanity_checks()

    clear; clc;

    warnState = bootstrap();
    cleanupObj = onCleanup(@() cleanup(warnState));

    sim = sim_config_sanity();

    pts = [ ...
        0.00  0.00; ...
        0.05  0.00; ...
        0.00  0.05; ...
        0.05  0.05; ...
       -0.05  0.05; ...
    ];

    warm = struct();

    fprintf('Running %d sanity solves...\n', size(pts,1));

    for i = 1:size(pts,1)
        params = struct('H0_1', pts(i,1), 'H0_2', pts(i,2));

        fprintf('\n[%d/%d] H0=[%+.3f,%+.3f]\n', i, size(pts,1), pts(i,1), pts(i,2));

        if isfield(warm,'result') && isfield(warm.result,'sol')
            warm.fromParams = prevParams;
        else
            warm = struct();
        end

        [result, meta] = solveAtParams(params, sim, warm);

        assert(all(isfinite(result.sol.y(:))), 'NaN/Inf in solution');
        assert(meta.BCmax < 5e-5, 'BCmax too large for sanity: %.2e', meta.BCmax);
        assert(meta.DEmax < 5e-1, 'DEmax too large for sanity: %.2e', meta.DEmax);

        fprintf('  label=%s | mesh=%d | BC=%.2e | DE=%.2e\n', meta.label, meta.mesh, meta.BCmax, meta.DEmax);

        ascii_profile(result.sol);

        warm = struct('result', result);
        prevParams = params;
    end

    fprintf('\nSanity checks complete.\n');
end

function sim = sim_config_sanity()
    sim = struct();

    sim.MP = struct('A',0.75,'V',0.72,'KG',0,'KA',1,'KB',1);
    [sim.MP.aS, sim.MP.bS] = computePhaseScales(sim.MP.A);

    sim.TH = struct();
    sim.TH.delta     = 0.015;
    sim.TH.BCmax     = 5e-5;
    sim.TH.DEmaxHard = 5e-1;
    sim.TH.rMin      = 1e-4;
    sim.TH.minH0Step = 0.02;
    sim.TH.opts = bvpset('RelTol',1e-5,'AbsTol',1e-7,'NMax',900);

    sim.SP = struct();
    sim.SP.ModelVersion   = "sanity";
    sim.SP.Verbose        = true;
    sim.SP.LogToFile      = false;
    sim.SP.SaveHomotopy   = false;
    sim.SP.CoarsenOnStall = true;
    sim.SP.MaxIters       = 1;
end

function ascii_profile(sol)
    x = linspace(sol.x(1), sol.x(end), 11);
    Y = deval(sol, x);

    rA = Y(4,:);  zA = Y(5,:);
    rB = Y(13,:); zB = Y(14,:);

    fprintf('  sample rA(zA) and rB(zB):\n');
    fprintf('    i   x     rA     zA     |    rB     zB\n');
    for i = 1:numel(x)
        fprintf('   %2d  %.3f  %+ .3f  %+ .3f  |  %+ .3f  %+ .3f\n', i, x(i), rA(i), zA(i), rB(i), zB(i));
    end
end

function warnState = bootstrap()
    root = pwd;
    restoredefaultpath;
    addpath(fullfile(root,'src'));
    addpath(fullfile(root,'src','utils'));
    addpath(fullfile(root,'src','quadtree'));
    addpath(fullfile(root,'src','model'));

    warning('off','MATLAB:singularMatrix');
    warning('off','MATLAB:nearlySingularMatrix');
    warning('off','MATLAB:rankDeficientMatrix');
    warnState = warning;
end

function cleanup(warnState)
    warning(warnState);
end