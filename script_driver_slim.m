%% Script Driver
%% run_initial_sweep.m

warnState = bootstrap();
sim = sim_config;
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

    % ---- thresholds (gates) ----
    TH.BCmax     = 1e-6;
    TH.DEmaxHard = 2e-1;
    TH.rMin      = 1e-3;

    % ---- solver knobs (not part of physics hash) ----
    TH.delta = 0.01;
    TH.opts  = bvpset( ...
        'RelTol',1e-6, ...
        'AbsTol',1e-8, ...
        'NMax',1500);
    TH.delta_list = [0.01, 0.015, 0.02, 0.008, 0.005];
    TH.minH0Step = 0.01;

    % ---- physical parameters (global for this run) ----
    MP.A  = 0.50;
    MP.V  = 0.72;
    MP.KA = 1.0;
    MP.KB = 1.0;
    MP.KG = 0.0;

    % derived phase scales from A
    [MP.aS, MP.bS] = computePhaseScales(MP.A);
    
    % pack it into a simulation struct
    S = struct('SP',SP,'TH',TH,'MP',MP);

end

%%
% bootstrap.m
% -------------------------------------------------------------
% One-stop setup for the vesicle simulation project.
% Run this once per MATLAB session before launching the driver.
% -------------------------------------------------------------

function warnState = bootstrap
    % Resolve project root from this file, not pwd
    projRoot = fileparts(mfilename('fullpath'));

    restoredefaultpath; rehash toolboxcache;
    addpath(genpath(fullfile(projRoot,'src')));
    addpath(fullfile(projRoot,'bvp6c-solver'));
    addpath(fullfile(projRoot,'initial-shapes'));

    simDir = fullfile(projRoot,'SimResults');
    if ~exist(simDir,'dir'), mkdir(simDir); end

    import_initial_shapes_into_catalog(fullfile(projRoot,'initial-shapes'), simDir);

    disp('MATLAB environment initialized for Vesicle Simulation Project.');
    ver

    warnState = warning('off','MATLAB:bvp6c:RelTolNotMet');
end

function import_initial_shapes_into_catalog(ishapesDir, simDir)
    T = catalog_load(simDir);
    files = dir(fullfile(ishapesDir,'SIM_Node_*.mat'));

    for k=1:numel(files)
        f = fullfile(files(k).folder, files(k).name);

        tok = regexp(files(k).name,'SIM_Node_(\d+)_(\d+)_(-?\d+)_(-?\d+)_(-?\d+)_','tokens','once');
        if isempty(tok), continue; end

        A  = str2double(tok{1})/100;
        V  = str2double(tok{2})/100;
        KG = str2double(tok{3});
        KA = str2double(tok{4});
        KB = str2double(tok{5});

        params = struct('A',A,'V',V,'KG',KG,'KA',KA,'KB',KB);

        entry = struct('params',params, ...
                       'meta',struct('type',"seed",'source',string(files(k).name)));

        % Unique per seed file (prevents accidental merging)
        seedKey = struct('model_version',"seed-v1", ...
                         'A',A,'V',V,'KG',KG,'KA',KA,'KB',KB, ...
                         'source', string(files(k).name));
        seedHash = simpleDataHash(seedKey, 'SHA-256');

        T = catalog_append(simDir, seedHash, entry); %#ok<NASGU>
    end
end

% cleanup.m
% -------------------------------------------------------------
% removes added paths
% -------------------------------------------------------------

function cleanup(warnState)
    
    warning(warnState);

    % Remove your source folders
    rmpath(genpath(fullfile(pwd,'src')));
    rmpath(fullfile(pwd,'bvp6c-solver'));
    rmpath(fullfile(pwd,'initial-shapes'));
end