% -------------------------------------------------------------------------
% EXTRACTED HELPER for "appendCatalogRow"
%   - Source: sim_driver_quad_tree_full.m
%   - Extracted: 2025-10-11 11:50:12
%   - Sub-helpers (nested functions) are retained in this file.
% -------------------------------------------------------------------------

function appendCatalogRow(S, row)
    L = load(S.paths.catalog_mat,'T'); 
    T = L.T;
    T = [T; struct2table(row)];
    save(S.paths.catalog_mat,'T','-v7.3');
end