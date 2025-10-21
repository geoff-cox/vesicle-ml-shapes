% -------------------------------------------------------------------------
% EXTRACTED HELPER for "traceBoundary"
%   - Source: sim_driver_quad_tree_full.m
%   - Extracted: 2025-10-11 11:50:12
%   - Sub-helpers (nested functions) are retained in this file.
% -------------------------------------------------------------------------

function poly = traceBoundary(S, pa, pb, cache)
    % Return a polyline of boundary points between two morphology classes.
    % pa, pb: 1x2 endpoints with different labels (already solved by caller).
    % cache:  solve cache (updated in place).
    %
    % Heuristics (tune in practice):
    maxSteps = 200;
    ds       = 0.25;     % predictor step in parameter space (H0 units)
    epsBis   = 0.02;     % bisection tolerance along normal
    pad      = 0.05;     % small offset to start bracketing around prediction
    stopBox  = [S.grid.A_lim; S.grid.B_lim]; % bounds

    % 1) Find initial boundary point p* by bisection on segment pa->pb
    [pStar, cache] = edgeBisectionToSwitch(S, pa, pb, cache, epsBis);

    poly = pStar;   % Nx2 list of boundary points
    % Initial normal is along the mixed edge
    n = unit(pb - pa);              % normal ~ direction across boundary
    t = rot90ccw(n);                % tangent = +90deg rotation

    for k = 1:maxSteps
        % 2) Predictor: step along tangent
        pPred = pStar + ds * t;

        % keep inside box
        if ~inBox(pPred, stopBox), break; end

        % 3) Corrector: bracket along normal around pPred, then bisection to label switch
        [pNext, nNew, success, cache] = correctToBoundary(S, pPred, n, pad, epsBis, cache);
        if ~success, break; end

        % 4) Append and update frame
        poly(end+1,:) = pNext; %#ok<AGROW>

        % New tangent: perpendicular to updated normal
        t = rot90ccw(nNew);
        pStar = pNext;
        n = nNew;
    end
end