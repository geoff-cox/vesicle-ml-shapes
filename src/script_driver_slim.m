%% Script Driver — Quadtree Parameter-Space Exploration
%
% Entry point for vesicle shape simulations.
% Bootstraps paths, loads seed shapes, runs the adaptive quadtree sweep,
% then restores the MATLAB environment.
%
% Usage:
%   >> cd src
%   >> script_driver_slim
%

warnState = bootstrap();
sim = sim_config();
sim_explore_H0_quad_tree(sim);
cleanup(warnState);

function S = sim_config()

    % ---- simulation settings ----
    SP.MaxIters = 1e5;
    SP.ModelVersion = "BVP-v3.1";
    SP.LogToFile = true;
    SP.Verbose = true;
    SP.saveHomotopy = true;
    SP.CoarsenOnStall = true;

    % ---- hysteresis / multi-solution ----
    % Set SP.BranchTag to a non-empty string (e.g., "upper", "lower",
    % "path-A") to explore an independent solution branch at every
    % (H0_1, H0_2).  Each branch tag produces a distinct hash, so
    % multiple solutions can coexist in the catalog for the same
    % physics parameters.  Leave empty (default) for single-branch mode.
    SP.BranchTag = "path-C";

    % ---- thresholds (gates) ----
    TH.BCmax     = 1e-6;
    TH.DEmaxHard = 2e-1;
    TH.rMin      = 1e-3;

    % ---- solver knobs (not part of physics hash) ----
    TH.delta = 1e-3;
    TH.opts  = bvpset( ...
        'RelTol',1e-6, ...
        'AbsTol',1e-8, ...
        'NMax',1500);
    TH.delta_list = 1e-3;
    TH.minH0Step = 0.01;
    TH.useLegacy = false;
    TH.poleDeg = 2;

    % ---- physical parameters (global for this run) ----
    MP.A  = 0.50;
    MP.V  = 0.72;
    MP.KA = 1.0;
    MP.KB = 1.0;
    MP.KG = 0.0;

    % derived phase scales from A
    [MP.aS, MP.bS] = computePhaseScales(MP.A);

    % ---- apply default overrides ----
    SP.H0Bounds = [-1 1; -1 1];     % [H0_1min H0_1max; H0_2min H0_2max]
    % SP.maxDepth = 7;
    % SP.maxCells = 4000;
    % SP.eTol     = 5e-3;
    % SP.pTol     = 5e-3;
    % SP.shapeTau = 0.08;

    % pack it into a simulation struct
    S = struct('SP',SP,'TH',TH,'MP',MP);

end

%%
% bootstrap  — one-stop setup for the vesicle simulation project.
% Run this once per MATLAB session before launching the driver.

function warnState = bootstrap()
    % Resolve src/ root from this file's location
    srcRoot = fileparts(mfilename('fullpath'));

    restoredefaultpath; rehash toolboxcache;
    addpath(genpath(srcRoot));

    % sim-results lives next to src/
    projRoot = fileparts(srcRoot);
    simDir   = fullfile(projRoot, 'sim-results');
    if ~exist(simDir,'dir'), mkdir(simDir); end

    ishapesDir = fullfile(srcRoot, 'initial-shapes');
    import_initial_shapes_into_catalog(ishapesDir, simDir);

    disp('MATLAB environment initialized for Vesicle Simulation Project.');

    warnState = warning('off','MATLAB:bvp6c:RelTolNotMet');
end

function import_initial_shapes_into_catalog(ishapesDir, simDir)
    T = catalog_load(simDir);
    files = dir(fullfile(ishapesDir,'SIM_Node_*.mat'));

    for k = 1:numel(files)
        tok = regexp(files(k).name, ...
            'SIM_Node_(\d+)_(\d+)_(-?\d+)_(-?\d+)_(-?\d+)_','tokens','once');
        if isempty(tok), continue; end

        A  = str2double(tok{1})/100;
        V  = str2double(tok{2})/100;
        KG = str2double(tok{3});
        KA = str2double(tok{4});
        KB = str2double(tok{5});

        params = struct('A',A,'V',V,'KG',KG,'KA',KA,'KB',KB);
        entry  = struct('params',params, ...
                        'meta',struct('type',"seed",'source',string(files(k).name)));

        seedKey = struct('model_version',"seed-v1", ...
                         'A',A,'V',V,'KG',KG,'KA',KA,'KB',KB, ...
                         'source', string(files(k).name));
        seedHash = simpleDataHash(seedKey, 'SHA-256');

        T = catalog_append(simDir, seedHash, entry); %#ok<NASGU>
    end
end

% cleanup — restore MATLAB path after a run.

function cleanup(warnState)
    warning(warnState);
    srcRoot = fileparts(mfilename('fullpath'));
    rmpath(genpath(srcRoot));
end