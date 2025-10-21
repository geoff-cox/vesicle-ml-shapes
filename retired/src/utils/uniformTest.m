% -------------------------------------------------------------------------
% EXTRACTED HELPER for "uniformTest"
%   - Source: sim_driver_quad_tree_full.m
%   - Extracted: 2025-10-11 11:50:11
%   - Sub-helpers (nested functions) are retained in this file.
% -------------------------------------------------------------------------

function [uniform, mixedEdges] = uniformTest(C, eTol, pTol, tau)
    % Decide whether a cell is "uniform" (same morphology) or "mixed".
    % Heuristics:
    %  1) all corner labels equal --> uniform
    %  2) else, if max energy/pressure difference small AND shape distances small --> uniform
    %  3) track which edges are mixed to seed boundary tracing

    labs = C.cornerLabel;
    if all(labs == labs(1))
        uniform = true; mixedEdges = zeros(0,2); return;
    end

    eSpread = max(C.cornerEnergy)   - min(C.cornerEnergy);
    pSpread = max(C.cornerPressure) - min(C.cornerPressure);

    % Optional: include shape descriptor distances if you store them in meta
    % Here, we only use labels + E/P spreads
    shapeSmall = true; % replace with real test if you export descriptors

    uniform = (eSpread < eTol) && (pSpread < pTol) && shapeSmall;

    % Which edges have differing labels?
    edges = [1 2; 2 3; 3 4; 4 1]; % SW-SE, SE-NE, NE-NW, NW-SW
    mixedEdges = [];
    if ~uniform
        for k = 1:4
            if labs(edges(k,1)) ~= labs(edges(k,2))
                mixedEdges(end+1,:) = edges(k,:); %#ok<AGROW>
            end
        end
    end
end