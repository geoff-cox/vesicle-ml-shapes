% -------------------------------------------------------------------------
% EXTRACTED HELPER for "edgeBisectionToSwitch"
%   - Source: sim_driver_quad_tree_full.m
%   - Extracted: 2025-10-11 11:50:12
%   - Sub-helpers (nested functions) are retained in this file.
% -------------------------------------------------------------------------

function [pStar, cache] = edgeBisectionToSwitch(S, pa, pb, cache, epsBis)
    [la, cache] = getLabel(S, pa, cache);
    [lb, cache] = getLabel(S, pb, cache);
    if la == lb, error('edgeBisectionToSwitch: endpoints share label'); end
    a = pa; b = pb;
    for it = 1:40
        m = 0.5*(a+b);
        [lm, cache] = getLabel(S, m, cache);
        if lm == la, a = m; else, b = m; end
        if norm(a-b) < epsBis, break; end
    end
    pStar = 0.5*(a+b);
end