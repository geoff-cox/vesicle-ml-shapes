% -------------------------------------------------------------------------
% EXTRACTED HELPER for "inBox"
%   - Source: sim_driver_quad_tree_full.m
%   - Extracted: 2025-10-11 11:50:12
%   - Sub-helpers (nested functions) are retained in this file.
% -------------------------------------------------------------------------

function tf = inBox(p, box)
% box = [ [Amin Amax]; [Bmin Bmax] ]
    tf = (p(1)>=box(1,1) && p(1)<=box(1,2) && p(2)>=box(2,1) && p(2)<=box(2,2));
end