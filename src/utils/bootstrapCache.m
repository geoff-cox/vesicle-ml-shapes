% -------------------------------------------------------------------------
% EXTRACTED HELPER for "bootstrapCache"
%   - Source: sim_driver_quad_tree_full.m
%   - Extracted: 2025-10-11 11:50:11
%   - Sub-helpers (nested functions) are retained in this file.
% -------------------------------------------------------------------------

function cache = bootstrapCache(S)
    cache = initSolveCache();
    seeds = [ 0 0;
              1 0; -1 0; 0 1; 0 -1;
              1 1; -1 1; 1 -1; -1 -1 ];
    for i = 1:size(seeds,1)
        S.H0 = seeds(i,:);
        try
            [~, meta] = solveAtParams(S, cache);
            cache = meta.cache;
            if S.plot.enable
                plot(S.plot.ax, S.H0(1), S.H0(2), '.', 'MarkerSize',10); drawnow limitrate;
            end
        catch ME
            fprintf('Bootstrap skip at (%g,%g): %s\n', S.H0(1), S.H0(2), ME.message);
        end
    end
end