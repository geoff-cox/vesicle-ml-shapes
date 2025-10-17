% -------------------------------------------------------------------------
% EXTRACTED HELPER for "unit"
%   - Source: sim_driver_quad_tree_full.m
%   - Extracted: 2025-10-11 11:50:12
%   - Sub-helpers (nested functions) are retained in this file.
% -------------------------------------------------------------------------

function v = unit(v), n = norm(v); if n==0, return; end; v = v./n; end

function r = rot90ccw(v), r = [ -v(2), v(1) ]; end

% --- change signature to return cache too
function [lbl, cache] = getLabel(S, p, cache)
    S.H0 = p;
    [~, meta] = solveAtParams(S, cache);
    lbl   = meta.label;
    cache = meta.cache;   % propagate
end