function sol = continuation_towards_H0(sol0, H0_from, H0_to, sim, stepMax)
% CONTINUATION_TOWARDS_H0  Walks from H0_from -> H0_to in small steps, reusing each solve as init.
% Inputs:
%   sol0     : bvp6c solution struct to start from
%   H0_from  : [H0_1 H0_2] of the warm-start
%   H0_to    : [H0_1 H0_2] target
%   sim      : struct with fields MP (A,V,KA,KB,KG) and TH (delta, opts)
%   stepMax  : maximum step length in H0 per subsolve (e.g., 0.15)
% Output:
%   sol      : bvp6c solution at (approximately) H0_to (last successful step)

    if nargin < 5 || isempty(stepMax), stepMax = 0.15; end

    Hfrom = H0_from(:)'; Hto = H0_to(:)'; 
    d = norm(Hto - Hfrom);
    if ~isfinite(d) || d==0, sol = sol0; return; end

    nSteps = max(1, ceil(d/stepMax));
    path   = Hfrom + (0:nSteps).'/nSteps .* (Hto - Hfrom);
    sol    = sol0;

    % fixed physics for this run
    MP = sim.MP;  [aS,bS] = computePhaseScales(MP.A);
    delta0 = sim.TH.delta;
    opts0  = sim.TH.opts;

    % suppress bvp6c diagnostic warnings
    warnState = warning('off','MATLAB:bvp6c:RelTolNotMet');

    for i = 2:size(path,1)  % path(1,:) is Hfrom
        H0 = path(i,:);
        Par = struct('H0',H0, 'A',MP.A,'V',MP.V,'KA',MP.KA,'KB',MP.KB,'KG',MP.KG, ...
                     'aS',aS,'bS',bS, 'delta',delta0);

        odefun = @(s,y,lam) BendV_Lag_EIGp_DE_impl(s,y,lam,Par);
        bcfun  = @(ya,yb,lam) BendV_Lag_EIGp_BC_impl(ya,yb,lam,Par);

        fprintf('continuation_towards_H0: attempting a bridge solve at %s with δ = %g\n',mat2str(H0),delta0)
        try
           sol = bvp6c(odefun, bcfun, sol, opts0);  % reuse last sol as initSol
           fprintf('continuation_towards_H0: solution accepted\n')
        catch
            fprintf('continuation_towards_H0: solution rejected due to error\n')
            % Bail out early; caller can still try main attempt ladder from the last good sol.
            break
        end
    end
    % restore original warning state after success
    warning(warnState);
end
