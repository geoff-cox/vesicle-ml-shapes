function [aS, bS, S_star] = computePhaseScales(A)
    % COMPUTEPHASESCALES  Map area ratio A ∈ [0,1] to phase arc-length scalings.
    % S_star = acos(1 - 2A) ∈ [0, π], aS = S_star/π, bS = (S_star - π)/π.
    arguments
        A (1,1) double {mustBeFinite}
    end
    S_star = acos(1 - 2*A);
    aS = S_star/pi;
    bS = (S_star - pi)/pi;  % note: negative; used with SB = bS*S + π
end