% -------------------------------------------------------------------------
% EXTRACTED HELPER for "initialGuessFromFile"
%   - Source: sim_driver_quad_tree_full.m
%   - Extracted: 2025-10-11 11:50:12
%   - Sub-helpers (nested functions) are retained in this file.
% -------------------------------------------------------------------------

function initSol = initialGuessFromFile(S, ~)
    p1 = fullfile(S.paths.mats, S.initialGuess);
    if exist(p1,'file')~=2
        p1 = fullfile('InitialShapes', S.initialGuess);
    end
    tmp = load(p1);
    initSol = tmp.Version(1).Solution;
end