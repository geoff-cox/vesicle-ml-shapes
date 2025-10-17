% -------------------------------------------------------------------------
% EXTRACTED HELPER for "local_min_radius_interior"
%   - Source: sim_driver_quad_tree_full.m
%   - Extracted: 2025-10-11 11:50:11
%   - Sub-helpers (nested functions) are retained in this file.
% -------------------------------------------------------------------------

function rmin = local_min_radius_interior(sol)
% Min radius away from poles (ignore s=0 and s=pi so r=0 at poles is allowed)
    rA = sol.y(4,:);   rB = sol.y(13,:);
    if numel(rA) >= 3
        rA = rA(2:end-1);
    end
    if numel(rB) >= 3
        rB = rB(2:end-1);
    end
    rmin = min([rA(:); rB(:)]);
    if isempty(rmin)
        rmin = Inf;  % degenerate mesh; don't fail geometry gate
    end
end