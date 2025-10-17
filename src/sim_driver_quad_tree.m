% -------------------------------------------------------------------------
% CLEANED MAIN for "sim_driver_quad_tree"
%   - Source: sim_driver_quad_tree_full.m
%   - Extracted: 2025-10-11 11:50:11
%   - Sub-helpers (nested functions) are retained in this file.
% -------------------------------------------------------------------------

function sim_driver_quad_tree()
% PHASEDIAGRAMDRIVER
% New, minimal driver for sampling the (H0_1, H0_2) plane using:
%   1) Quadtree refinement of parameter cells, and
%   2) Lightweight boundary tracing between morphology classes.
%
% Numerics: bvp6c
% ODE/BC:   BendV_Lag_EIGp_DE_lam / BendV_Lag_EIGp_BC_lam   (your existing files)

    % --- Paths (adjust if needed)
    addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
    addpath('bvp6c-solver');
    addpath('InitialShapes');

    % --- Build config/state and output folders/log
    sim = makeSimConfig();           % all knobs in one place (mode, grid limits, solver opts)
    sim = initFoldersAndLogging(sim);  % sets S.paths.* and S.log.fid
    sim = initResultsIO(sim);

    % --- (Optional) initialize a figure for live coverage
    if sim.plot.enable
        sim.plot.fig = figure('Position',[40 40 1100 900],'Color','w');
        sim.plot.ax  = axes('Position',[0.08 0.08 0.86 0.86]); hold(sim.plot.ax,'on'); box(sim.plot.ax,'on');
        xlabel(sim.plot.ax,'H0^{(1)}'); ylabel(sim.plot.ax,'H0^{(2)}');
        title(sim.plot.ax,'Quadtree coverage'); grid(sim.plot.ax,'on');
    end

    cache = bootstrapCache(sim);
    cache = loadCacheIfExists(sim, cache);

    % --- Initialize quadtree over (H0_1, H0_2)
    QT = initQuadtree(sim);

    % --- Quadtree processing (solves corners, refines mixed cells, traces boundaries)
    qParams = struct('maxDepth', 4, ...     % refine up to 2^6 resolution per axis
                     'maxCells', 2000, ...
                     'eTol',     5e-3, ...  % energy uniformity threshold
                     'pTol',     5e-3, ...  % pressure uniformity threshold
                     'shapeTau', 0.08);     % (optional) shape-distance threshold if used
    processQuadtree(sim, QT, qParams, cache);

    % --- Done
    saveCache(sim, cache);
    logmsg(sim,'[DONE] Quadtree pass complete.');
    if sim.log.fid > 0, fclose(sim.log.fid); end
end