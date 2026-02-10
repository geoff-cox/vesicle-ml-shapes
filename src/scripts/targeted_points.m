% 0) Session bootstrap
if exist('bootstrap.m','file'), run('bootstrap.m'); end

% 1) Build sim + paths
sim = makeSimConfig();
sim.debug.short = true;                          % quick meshes
sim = initFoldersAndLogging(sim);
sim = initResultsIO(sim);

% Solver robustness tweaks
sim.solver.delta = max(sim.solver.delta, 0.02);  % wider Taylor buffer near poles
sim.solver.opts  = bvpset(sim.solver.opts,'NMax',4000);  % allow more Newton work

fprintf('[targeted_points] Writing to: %s\n', sim.paths.run);

% 2) Warm the cache by solving a few anchors
cache = bootstrapCache(sim);  % may skip some seeds; that's fine

% Always ensure (0,0) is solved first (best universal seed)
ensureSolved([0 0]);

% Also try a couple of diagonals (often easier than axes)
ensureSolved([1 1]);
ensureSolved([-1 -1]);

% 3) List of "hard" targets we want to reach
targets = [
    1  0;
   -1  0;
    0  1;
    2 -2;
   -2  2
];

% 4) Attempt continuation to each target
for i = 1:size(targets,1)
    pT = targets(i,:);
    fprintf('\n--- Target (%.2f, %.2f) ---\n', pT(1), pT(2));
    ok = climbTo(pT, 0.20);    % start with 0.20 step
    if ~ok, ok = climbTo(pT, 0.10); end  % retry with smaller steps
    if ~ok, ok = climbTo(pT, 0.05); end  % last resort
    if ~ok
        fprintf('FAIL: could not reach (%.2f, %.2f)\n', pT(1), pT(2));
    end
end

% 5) Close log
if isfield(sim,'log') && isfield(sim.log,'fid') && sim.log.fid > 0
    fclose(sim.log.fid);
end

% ---------- helpers (in-file local functions) ----------

function ok = ensureSolved(p)
    % Solve exactly at p (uses cache + attempt ladder inside solveAtParams)
    ok = true;
    try
        assignin('caller','sim',evalin('caller','sim'));  % bring sim/cache scope
        assignin('caller','cache',evalin('caller','cache'));
        sim = evalin('caller','sim'); cache = evalin('caller','cache');
        sim.H0 = p;
        [~, meta] = solveAtParams(sim, cache);
        cache = meta.cache;
        assignin('caller','cache',cache);
        fprintf('Anchor OK @ (%.2f, %.2f) | label=%d | E=%.5g | P=%.5g\n', ...
            p(1), p(2), meta.label, meta.E_total, meta.P_osm);
    catch ME
        ok = false;
        fprintf('Anchor FAIL @ (%.2f, %.2f): %s\n', p(1), p(2), ME.message);
    end
end

function ok = climbTo(pT, ds)
    % Two-leg continuation from the nearest solved cache point to pT,
    % stepping by 'ds' in parameter space and calling solveAtParams at each substep.
    ok = false;
    sim = evalin('caller','sim'); cache = evalin('caller','cache');

    p0 = pickSeed(cache, pT);
    if isempty(p0)
        % fall back to (0,0) if nothing else
        p0 = [0 0];
        if ~ensureSolved(p0), return; end
        sim = evalin('caller','sim'); cache = evalin('caller','cache');
    end

    % Path: (p0_1, p0_2) -> (pT_1, p0_2) -> (pT_1, pT_2)
    leg1 = ramp(p0(1), pT(1), ds); leg1 = [leg1(:), repmat(p0(2), numel(leg1),1)];
    leg2 = ramp(p0(2), pT(2), ds); leg2 = [repmat(pT(1), numel(leg2),1), leg2(:)];
    path = unique([leg1; leg2], 'rows', 'stable');

    try
        for k = 1:size(path,1)
            sim.H0 = path(k,:);
            % slight delta bump as we move can help near poles
            if mod(k,5)==0, sim.solver.delta = min(0.03, sim.solver.delta*1.25); end

            [~, meta] = solveAtParams(sim, cache);
            cache = meta.cache;  % keep cache hot
            fprintf(' step-> (%.3f, %.3f)  | label=%d | E=%.5g | P=%.5g\n', ...
                sim.H0(1), sim.H0(2), meta.label, meta.E_total, meta.P_osm);
        end
        assignin('caller','cache',cache);
        ok = true;
        fprintf('SUCCESS reached (%.2f, %.2f)\n', pT(1), pT(2));
    catch ME
        fprintf('Step FAIL near (%.3f, %.3f): %s\n', sim.H0(1), sim.H0(2), ME.message);
    end
end

function s = ramp(a, b, ds)
    if a == b, s = b; return; end
    n = max(2, ceil(abs(b-a)/ds));
    s = linspace(a,b,n);
end

function p0 = pickSeed(cache, pT)
    p0 = [];
    if ~isfield(cache,'keys') || isempty(cache.keys), return; end
    keys = cache.keys;    % Nx2 solved points
    % Prefer same signs (often same morphology basin), else nearest
    sameSign = all(sign(keys)==sign(pT),2);
    pool = keys(sameSign,:);
    if isempty(pool), pool = keys; end
    [~,j] = min(sum((pool - pT).^2,2));
    p0 = pool(j,:);
end

run('cleanup.m');