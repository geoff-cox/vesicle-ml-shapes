% -------------------------------------------------------------------------
% EXTRACTED HELPER for "initQuadtree"
%   - Source: sim_driver_quad_tree_full.m
%   - Extracted: 2025-10-11 11:50:11
%   - Sub-helpers (nested functions) are retained in this file.
% -------------------------------------------------------------------------

function QT = initQuadtree(S)
    root = makeCell(S.grid.A_lim(1), S.grid.A_lim(2), ...
                    S.grid.B_lim(1), S.grid.B_lim(2), 0);
    QT.root  = root;
    QT.queue = {root};
    % empty struct array with same fields as a cell
    QT.cells = repmat(root, 0, 1);
end