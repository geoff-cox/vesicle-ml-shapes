% -------------------------------------------------------------------------
% EXTRACTED HELPER for "saveCache"
%   - Source: sim_driver_quad_tree_full.m
%   - Extracted: 2025-10-11 11:50:12
%   - Sub-helpers (nested functions) are retained in this file.
% -------------------------------------------------------------------------

function saveCache(S, cache)
    f = fullfile(S.paths.run, 'cache.mat');
    save(f,'cache','-v7.3');
end