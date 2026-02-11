function script_sanity_checks(varargin)
% SCRIPT_SANITY_CHECKS  Fast validation of the exploration pipeline.
%
% Usage:
%   script_sanity_checks
%   script_sanity_checks('doPhysicsSmoke',true)

    p = inputParser;
    addParameter(p,'doPhysicsSmoke',false,@islogical);
    parse(p,varargin{:});
    doPhysicsSmoke = p.Results.doPhysicsSmoke;

    here = fileparts(mfilename('fullpath'));     % src/scripts
    srcRoot = fileparts(here);                    % src

    % Paths (minimal; avoid restoredefaultpath here)
    addpath(genpath(srcRoot));

    tmpRoot = fullfile(tempdir, "vesicle_sanity_" + string(java.util.UUID.randomUUID));
    simDir  = fullfile(tmpRoot, 'sim-results');
    mkdir(simDir); mkdir(fullfile(simDir,'hashed_results'));

    % --- build sim config ---
    sim = make_sim();
    sim.SP.SimDir        = simDir;
    sim.SP.H0Bounds      = [-0.2 0.2; -0.2 0.2];
    sim.SP.QTmaxDepth    = 2;
    sim.SP.QTmaxCells    = 200;
    sim.SP.MaxIters      = 50;
    sim.SP.CacheSaveEvery= 5;
    sim.SP.UseWarmStart  = false;          % mock run: don’t load seeds
    sim.SP.SolverFcn     = @mockSolve;     % fast solver stub

    fprintf('Sanity temp dir: %s\n', simDir);

    % --- run mock exploration ---
    sim_explore_H0_quad_tree(sim);

    % --- invariants ---
    assert(exist(fullfile(simDir,'cache.mat'),'file')==2, 'cache.mat was not created.');
    assert(exist(fullfile(simDir,'catalog.mat'),'file')==2, 'catalog.mat was not created.');

    T = catalog_load(simDir);
    assert(height(T) > 0, 'Catalog is empty after mock run.');
    assert(all(T.hash ~= ""), 'Some catalog hashes are empty.');
    assert(any(contains(string(T.entry{1}.meta.label), "mock")), 'Mock labels missing; mockSolve not used?');

    fprintf('[OK] mock exploration produced %d catalog rows.\n', height(T));

    % --- optional real solver smoke test (single point) ---
    if doPhysicsSmoke
        sim2 = make_sim();
        sim2.SP.SimDir = simDir;
        sim2.SP.UseWarmStart = true;
        sim2.SP.SolverFcn = @solveAtParams;

        params = struct('H0_1',0,'H0_2',0);
        warm = struct(); % solveAtParams will pull a seed
        [result, meta] = solveAtParams(params, sim2, warm); %#ok<NASGU>

        assert(isfield(meta,'BCmax') && isfield(meta,'DEmax'), 'meta missing diagnostics.');
        fprintf('[OK] physics smoke: BC=%.2e, DE=%.2e, mesh=%d\n', meta.BCmax, meta.DEmax, meta.mesh);
    end

    fprintf('All sanity checks passed.\n');

end

% ---------------- helpers ----------------
function sim = make_sim()
    SP.MaxIters = 1e5;
    SP.ModelVersion = "BVP-v3.1";
    SP.LogToFile = false;
    SP.Verbose = false;
    SP.saveHomotopy = false;
    SP.CoarsenOnStall = true;

    TH.BCmax     = 1e-6;
    TH.DEmaxHard = 2e-1;
    TH.rMin      = 1e-3;

    TH.delta = 0.01;
    TH.opts  = bvpset('RelTol',1e-6,'AbsTol',1e-8,'NMax',1500);
    TH.delta_list = [0.01, 0.015, 0.02, 0.008, 0.005];
    TH.minH0Step = 0.01;

    MP.A  = 0.50;
    MP.V  = 0.72;
    MP.KA = 1.0;
    MP.KB = 1.0;
    MP.KG = 0.0;
    [MP.aS, MP.bS] = computePhaseScales(MP.A);

    sim = struct('SP',SP,'TH',TH,'MP',MP);
end

function [result, meta] = mockSolve(params, sim, warm) %#ok<INUSD>
    % Cheap deterministic fake “solve” so you can validate the pipeline.
    H1 = params.H0_1; H2 = params.H0_2;

    % Minimal sol struct for saving; warm-start disabled in sanity run
    sol = struct();
    sol.x = linspace(0,pi,20);
    sol.y = zeros(18, numel(sol.x));
    sol.parameters = 0;

    result = struct('sol',sol,'mesh',numel(sol.x));

    meta = struct();
    meta.label = "mock_" + sprintf('%+.2f_%+.2f', H1, H2);
    meta.E     = H1^2 + H2^2;
    meta.P     = 0;
    meta.BCmax = 0;
    meta.DEmax = 0;
    meta.mesh  = result.mesh;
    meta.rMinAway = 1;
end
