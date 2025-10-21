% -------------------------------------------------------------------------
% EXTRACTED HELPER for "resultFileFor"
%   - Source: sim_driver_quad_tree_full.m
%   - Extracted: 2025-10-11 11:50:12
%   - Sub-helpers (nested functions) are retained in this file.
% -------------------------------------------------------------------------

function [path,hash] = resultFileFor(S)
    hash = simpleDataHash(makeSolveKey(S), 'SHA-256');
    path = fullfile(S.paths.results, [hash '.mat']);
end