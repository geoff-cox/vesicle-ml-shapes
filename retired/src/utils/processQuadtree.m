% -------------------------------------------------------------------------
% EXTRACTED HELPER for "processQuadtree"
%   - Source: sim_driver_quad_tree_full.m
%   - Extracted: 2025-10-11 11:50:11
%   - Sub-helpers (nested functions) are retained in this file.
% -------------------------------------------------------------------------

function processQuadtree(S, QT, params, cache)
    % Main loop: pop a cell, solve corners, test stopping, subdivide if needed.
    
    % Debug overrides
    if isfield(S,'debug') && isfield(S.debug,'short') && S.debug.short
        params.maxDepth = 2;           % very shallow refinement
        params.maxCells = 6;          % small number of cells
        S.solver.opts   = bvpset(S.solver.opts,'NMax', 500,'RelTol',5e-6);
    end

    % params: configuration for stopping criteria & max depth, e.g.:
    %   params.maxDepth, params.maxCells, params.eTol, params.pTol, params.shapeTau
    maxCells = defaultArg(params,'maxCells', 2000);
    eTol     = defaultArg(params,'eTol',     5e-3);
    pTol     = defaultArg(params,'pTol',     5e-3);
    tau      = defaultArg(params,'shapeTau', 0.08);
    maxDepth = defaultArg(params,'maxDepth', 6);

    cellCount = 0;

    progressbar = char('-'*ones(1,100));
    while ~isempty(QT.queue) && cellCount < maxCells
        dx = floor((cellCount+1)/maxCells*100);
        C = QT.queue{1};  QT.queue(1) = [];
        progressbar(1:dx) = '#';
        fprintf('Progress: %s\n', progressbar);
        logmsg(S, 'Progress:\nFinalized Cells: %i (max = %i)\n', cellCount+1, maxCells);

        % 1) Solve corners (if needed)
        failCount = 0;
        for i = 1:4
            if ~C.cornerSolved(i)
                try
                    S.H0 = C.corners(i,:);
                    [sol, meta] = solveAtParams(S, cache);
                    C.cornerSolved(i)   = true;
                    C.cornerLabel(i)    = meta.label;
                    C.cornerEnergy(i)   = meta.E_total;
                    C.cornerPressure(i) = meta.P_osm;
                    cache = meta.cache;

                    if S.plot.enable
                        plot(S.plot.ax, S.H0(1), S.H0(2), '.', 'MarkerSize',8);
                        drawnow limitrate;
                    end
                catch
                    failCount = failCount + 1;
                    C.cornerSolved(i) = false;   % mark unsolved
                    logmsg(S, 'Corner (%+.2f,%+.2f) failed; will try after subdivision.', S.H0(1), S.H0(2));
                end
            end
        end
        
        if failCount >= 2
            % too risky to judge; subdivide immediately
            [C1,C2,C3,C4] = subdivideCell(C);
            QT.queue(end+1:end+4) = {C1,C2,C3,C4};
            continue
        end

        % 2) Stopping test
        [uniform, mixedEdges] = uniformTest(C, eTol, pTol, tau);
        C.isUniform  = uniform;
        C.mixedEdges = mixedEdges;

        % 3) If uniform or max depth reached, finalize cell
        if uniform || C.depth >= maxDepth
            QT.cells = [QT.cells; C]; 
            cellCount = cellCount + 1;
            % Optional: attempt boundary tracing from each mixed edge discovered earlier
            continue
        end

        % 4) Otherwise subdivide into 4 children
        [C1,C2,C3,C4] = subdivideCell(C);
        QT.queue{end+1} = C1; 
        QT.queue{end+1} = C2; 
        QT.queue{end+1} = C3; 
        QT.queue{end+1} = C4; 

        cellCount = cellCount + 1;

        % 5) Kick off boundary tracing from each mixed edge (lightweight)
        for e = 1:size(mixedEdges,1)
            i = mixedEdges(e,1); j = mixedEdges(e,2);
            pa = C.corners(i,:);  pb = C.corners(j,:);
            try
                traceBoundary(S, pa, pb, cache); % see block 2
            catch ME
                fprintf('Boundary trace failed at [%g,%g]-[%g,%g]: %s\n', ...
                        pa(1),pa(2),pb(1),pb(2), ME.message);
            end
        end
    end
end