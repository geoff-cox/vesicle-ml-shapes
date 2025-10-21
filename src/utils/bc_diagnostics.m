function [BCmax, report] = bc_diagnostics(sol, bcfun)

    if isfield(sol,'parameters') && ~isempty(sol.parameters)
        lam = sol.parameters(:);
    end
    Ya = sol.y(:,1); Yb = sol.y(:,end);
    res = bcfun(Ya, Yb, lam);
    [BCmax, iMax] = max(abs(res));
    report.idx = iMax;
    report.res = res;

    % --- robust min radius AWAY FROM POLES ---
    s  = sol.x;                    % s in [0, pi]
    rA = abs(sol.y(4,:));          % alpha radius
    rB = abs(sol.y(13,:));         % beta radius

    % choose a small buffer near each pole; purely diagnostic
    hmean = mean(diff(s));
    % heuristics: at least a few mesh spacings, and ~1% of the domain
    buf = max(5*hmean, 0.01*pi);

    mask = (s > buf) & (s < (pi - buf));
    if any(mask)
        rMinAway = min( [ min(rA(mask)), min(rB(mask)) ] );
    else
        % if mesh is too coarse, fall back to whole domain (rare in practice)
        rMinAway = min( [ min(rA), min(rB) ] );
    end

    report.min_r = rMinAway;       % <-- use this in your acceptance gate
end