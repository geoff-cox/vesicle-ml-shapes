% -------------------------------------------------------------------------
% EXTRACTED HELPER for "initSolveCache"
%   - Source: sim_driver_quad_tree_full.m
%   - Extracted: 2025-10-11 11:50:11
%   - Sub-helpers (nested functions) are retained in this file.
% -------------------------------------------------------------------------

function cache = initSolveCache()
    % Simple cache for solved parameter points and their solutions/labels.
    cache = struct();
    cache.keys   = [];        % Nx2 matrix of params
    cache.labels = [];        % Nx1 labels
    cache.E      = [];        % Nx1 energy
    cache.P      = [];        % Nx1 pressure
    cache.sols   = {};        % cell array of bvp solutions (for reuse)
end