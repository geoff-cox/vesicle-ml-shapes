function [uniform, mixedEdges] = uniformTest(C, eTol, pTol, tau)
% Decide whether a cell is "uniform" (same morphology) or "mixed".
% Heuristics:
%   0) if any corner is unsolved → not uniform (forces subdivision)
%   1) all corner labels equal → uniform
%   2) else, if energy/pressure spreads small (and optional shape test) → uniform
% Returns:
%   uniform     : logical
%   mixedEdges  : K×2 indices into the corners (SW,SE,NE,NW) that differ in label

    % 0) require solved corners
    if isfield(C,'cornerSolved') && any(~C.cornerSolved)
        uniform = false; mixedEdges = zeros(0,2); return;
    end

    labs = C.cornerLabel;
    % 1) all equal?
    if all(labs == labs(1))
        uniform = true; mixedEdges = zeros(0,2); return;
    end

    % Guard against NaNs
    E = C.cornerEnergy;   E(~isfinite(E)) = inf;
    P = C.cornerPressure; P(~isfinite(P)) = inf;

    eSpread = max(E) - min(E);
    pSpread = max(P) - min(P);

    % Optional shape criterion (placeholder)
    shapeSmall = true; %#ok<NASGU>  % replace if you add descriptors and use 'tau'

    uniform = (eSpread <= eTol) && (pSpread <= pTol);  % keep strict but inclusive

    % Edges in SW,SE,NE,NW order
    edges = [1 2; 2 3; 3 4; 4 1];  % SW-SE, SE-NE, NE-NW, NW-SW
    mixedEdges = [];
    if ~uniform
        for k = 1:4
            if labs(edges(k,1)) ~= labs(edges(k,2))
                mixedEdges(end+1,:) = edges(k,:); %#ok<AGROW>
            end
        end
    end
end
