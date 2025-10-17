% -------------------------------------------------------------------------
% EXTRACTED HELPER for "makeSimConfig"
%   - Source: sim_driver_quad_tree_full.m
%   - Extracted: 2025-10-11 11:50:11
%   - Sub-helpers (nested functions) are retained in this file.
% -------------------------------------------------------------------------

function sim = makeSimConfig()
    % MAKESIMCONFIG  Build configuration/state struct.

    % ---- High-level SIM params ----
    A  = 0.50;
    V  = 0.72;
    KG = 0;
    KA = 1;
    KB = 1;
    S_star = acos(1-2*A);
    aS = S_star/pi;
    bS = (S_star - pi)/pi;

    sim.params = struct(...
        'A', A, 'V', V, 'KG', KG, 'KA', KA, 'KB', KB, ...
        'aS', aS, 'bS', bS...
        );

    % Starting initial-guess file (lives in InitialShapes/)
    sim.initialGuess = sprintf('SIM_Node_%s_%s_%s_%s_%s_+00_+00.mat', ...
        num2str(round(A*100)), num2str(round(V*100)), ...
        num2str(KG), num2str(KA), num2str(KB));

    ig_src = fullfile('InitialShapes', sim.initialGuess);
    assert(exist(ig_src,'file')==2, 'Initial guess not found: %s', ig_src);

    % Mode: 'new' | 'continue' | 'rescan'
    sim.mode = 'continue';           % change to 'continue' or 'rescan' as needed
    sim.contDate   = '15-Oct-2024';
    sim.rescanDate = '18-Oct-2024';

    sim.debug.short = true;  % short run mode

    sim.print2cmdwin = false;

    % Composite base title for MAT files
    sim.baseTitle = sprintf('SIM_Node_%s_%s_%s_%s_%s', ...
        num2str(round(A*100)), num2str(round(V*100)), ...
        num2str(KG), num2str(KA), num2str(KB));

    % ---- Phase-plane grid ----
    A_lim   = [-4, 4];
    B_lim   = [-4, 4];
    A_nodes = [0:A_lim(2), -1:-1:A_lim(1)];  % same ordering as original
    B_nodes = [0:B_lim(2), -1:-1:B_lim(1)];

    sim.grid = struct(...
        'A_lim', A_lim, 'B_lim', B_lim, ...
        'A_nodes',A_nodes,'B_nodes',B_nodes...
        );

    % ---- Constants / Directions ----
    sim.const.dirNames = [' E';' W';' N';' S';'NE';'SW';'NW';'SE'];
    sim.const.director = [ 1,0; -1,0; 0,1; 0,-1; 1,1; -1,-1; -1,1; 1,-1];

    % ---- Solver options ----
    opts            = bvpset('RelTol',1e-6,'Stats','on', 'NMax', 2000);
    delta           = 0.01;     % Taylor expansion width near poles
    maxArr          = 1000;
    stepTol         = 1/500;    % stopping criterion on step size magnitude
    stepGrow        = 40;       % matches your harmonic update scheme
    saveSolutions   = true;
    goodCountGrowThreshold = 100;

    sim.solver = struct(...
        'opts',                     opts, ...
        'delta',                    delta, ...
        'maxArr',                   maxArr, ...
        'stepTol',                  stepTol, ...
        'stepGrow',                 stepGrow, ...
        'saveSolutions',            saveSolutions, ...
        'goodCountGrowThreshold',   goodCountGrowThreshold ...
        );

    % ---- Plotting ----
    sim.plot.enable = false;

    % ---- Runtime mutable fields (filled later) ----
    sim.paths = struct();
    sim.log   = struct('fid',-1,'path','');
    sim.H0    = [];
end