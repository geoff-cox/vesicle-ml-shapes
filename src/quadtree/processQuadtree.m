% === FILE: processQuadtree.m ===
% === patch ===
% Only showing the critical initialization lines to change.
% =============

function [task, cache] = processQuadtree(cache, T, MP)
% PROCESSQUADTREE  Adaptive quadtree task scheduler over (H0_1,H0_2).
% Returns:
%   task.params.H0_1, task.params.H0_2   for the next solve
%   cache updated in-place
%
% Assumptions:
%   - T is the unified catalog table with column "entry" (cell of structs).
%   - Driver is responsible for persisting cache via saveCache().

    if nargin < 3, MP = struct(); end
    task = [];
    cache = ensure_defaults(cache);

    if ~isfield(cache,'QT') || isempty(cache.QT) || ~isfield(cache.QT,'queue') || isempty(cache.QT.queue)
        C0 = makeCell(cache.config.bounds, 0);
        cache.QT.queue = {C0};
        cache.QT.cells = C0([]);   % empty struct array of same type
    end

    tol = 1e-12;

    deferredCount = 0;
    queueCount = 0;
    while ~isempty(cache.QT.queue)
        if deferredCount == 0
            queueCount = numel(cache.QT.queue);
        end

        C = cache.QT.queue{1};
        cache.QT.queue(1) = [];

        [C, anyUnknown] = refresh_corners_from_catalog(C, T, MP, tol);

        if anyUnknown
            k = find(~C.cornerSolved,1,'first');

            % If driver has recorded cooldown failures with H0_1/H0_2, respect them:
            if isfield(cache,'failures') && isfield(cache,'iter') && ~isempty(cache.failures)
                k = pick_unblocked_corner(C, cache);
                if isempty(k)
                    % All candidate corners in this cell are currently cooldown-blocked.
                    % Defer this cell for now, but if *every* cell in the queue is deferred,
                    % fall back to scheduling a blocked corner rather than exiting early.
                    deferredCount = deferredCount + 1;
                    if deferredCount >= queueCount
                        % Check if all unsolved corners are permanently
                        % blocked (blockUntil == Inf).  If so, mark them
                        % solved with a failure label and let the cell
                        % proceed to the uniform test instead of returning
                        % a task that the driver will endlessly skip.
                        allPerm = all_unsolved_permanent(C, cache);
                        if allPerm
                            for jj = find(~C.cornerSolved)'
                                C.cornerSolved(jj) = true;
                                C.cornerLabel(jj) = "failed";
                                C.cornerEnergy(jj) = NaN;
                                C.cornerPressure(jj) = NaN;
                            end
                            % All unsolved corners are permanently blocked and have
                            % just been marked as solved with label "failed". Proceed
                            % to the uniform-test / refinement logic in this same
                            % iteration (avoid re-calling refresh_corners_from_catalog
                            % on this cell, which would reset these labels).
                            deferredCount = 0;
                            anyUnknown = false;
                        end

                        % Fallback: choose the first unsolved corner even if blocked, so that
                        % the driver can advance iter and eventually let cooldown expire.
                        k_fallback = find(~C.cornerSolved,1,'first');
                        if isempty(k_fallback)
                            % Defensive: no unsolved corners found; exit loop as before.
                            break
                        end
                        params = struct('H0_1',C.corners(k_fallback,1), ...
                                        'H0_2',C.corners(k_fallback,2));
                        % Re-queue this cell at the front so it will be revisited.
                        cache.QT.queue = [{C}, cache.QT.queue];
                        task = struct('params',params);
                        return
                    else
                        % Defer this cell; move on to other queued cells.
                        cache.QT.queue{end+1} = C;
                        continue
                    end
                end
            end

            params = struct('H0_1',C.corners(k,1),'H0_2',C.corners(k,2));
            cache.QT.queue = [{C}, cache.QT.queue];
            task = struct('params',params);
            return
        else
            [uniform,mixedEdges] = uniformTest(C, cache.config.eTol, cache.config.pTol, cache.config.shapeTau);
            C.isUniform = uniform;
            C.mixedEdges = mixedEdges;

            if uniform || C.depth >= cache.config.maxDepth || numel(cache.QT.cells) >= cache.config.maxCells
                cache.QT.cells = [cache.QT.cells; C];
            else
                [C1,C2,C3,C4] = subdivideCell(C);
                cache.QT.queue(end+1:end+4) = {C1,C2,C3,C4};
            end

            deferredCount = 0;
        end
    end
end

% ---------------- defaults ----------------
function cache = ensure_defaults(cache)
    if ~isfield(cache,'config'), cache.config = struct(); end
    if ~isfield(cache.config,'bounds')
        cache.config.bounds = [-1 1; -1 1];   % [H0_1min H0_1max; H0_2min H0_2max]
    end
    cache.config.maxDepth = defaultArg(cache.config,'maxDepth', 7);
    cache.config.maxCells = defaultArg(cache.config,'maxCells', 4000);
    cache.config.eTol     = defaultArg(cache.config,'eTol',     5e-3);
    cache.config.pTol     = defaultArg(cache.config,'pTol',     5e-3);
    cache.config.shapeTau = defaultArg(cache.config,'shapeTau', 0.08);
end

% ---------------- failure-aware corner selection ----------------
function k = pick_unblocked_corner(C, cache)
    k = [];
    it = cache.iter;
    tol = 1e-12;
    for j = 1:size(C.corners,1)
        if C.cornerSolved(j), continue; end
        h1 = C.corners(j,1); h2 = C.corners(j,2);
        if ~isBlocked(cache.failures, it, h1, h2, tol)
            k = j; return
        end
    end
end

function tf = isBlocked(F, it, h1, h2, tol)
    tf = false;
    if isempty(F), return; end
    if ~isstruct(F), return; end
    if ~all(isfield(F, {'H0_1','H0_2','blockUntil'})), return; end

    for i=1:numel(F)
        if abs(double(F(i).H0_1)-h1)<tol && abs(double(F(i).H0_2)-h2)<tol
            if double(F(i).blockUntil) > it
                tf = true;
            end
            return
        end
    end
end

function tf = all_unsolved_permanent(C, cache)
% True when every unsolved corner is permanently blocked (blockUntil==Inf).
    tf = false;
    unsolved = find(~C.cornerSolved);
    if isempty(unsolved), return; end
    tol = 1e-12;
    F = cache.failures;
    for j = unsolved'
        h1 = C.corners(j,1); h2 = C.corners(j,2);
        perm = false;
        if isstruct(F) && ~isempty(F) && all(isfield(F,{'H0_1','H0_2','blockUntil'}))
            for i=1:numel(F)
                if abs(double(F(i).H0_1)-h1)<tol && abs(double(F(i).H0_2)-h2)<tol
                    perm = isinf(double(F(i).blockUntil));
                    break
                end
            end
        end
        if ~perm, return; end
    end
    tf = true;
end

% ---------------- catalog lookup ----------------
function [C, anyUnknown] = refresh_corners_from_catalog(C, T, MP, tol)
    anyUnknown = false;
    for i=1:4
        H1 = C.corners(i,1); H2 = C.corners(i,2);
        [solved,label,E,P] = lookup_in_catalog(T,H1,H2,MP,tol);
        C.cornerSolved(i)=solved;
        C.cornerLabel(i)=string(label);
        C.cornerEnergy(i)=E;
        C.cornerPressure(i)=P;
        if ~solved, anyUnknown=true; end
    end
end

function [solved,label,E,P] = lookup_in_catalog(T,H1,H2,MP,tol)
    solved=false; label=""; E=NaN; P=NaN;
    if isempty(T) || ~any(T.Properties.VariableNames=="entry"), return; end

    entries = T.entry;
    if isstruct(entries), entries = num2cell(entries);
    elseif ~iscell(entries)
        try, entries = num2cell(entries); catch, return; end
    end

    H1v = cellfun(@(e) g(e,'params','H0_1'), entries);
    H2v = cellfun(@(e) g(e,'params','H0_2'), entries);

    Av  = cellfun(@(e) g(e,'params','A'),  entries);
    Vv  = cellfun(@(e) g(e,'params','V'),  entries);
    KAv = cellfun(@(e) g(e,'params','KA'), entries);
    KBv = cellfun(@(e) g(e,'params','KB'), entries);
    KGv = cellfun(@(e) g(e,'params','KG'), entries);

    tolH0   = tol;
    tolPhys = 1e-12;
    if isstruct(MP) && isfield(MP,'tolPhys') && MP.tolPhys > 0
        tolPhys = MP.tolPhys;
    end

    hit = abs(H1v-H1)<=tolH0 & abs(H2v-H2)<=tolH0 ...
        & abs(Av - MP.A)  <= tolPhys ...
        & abs(Vv - MP.V)  <= tolPhys ...
        & abs(KAv- MP.KA) <= tolPhys ...
        & abs(KBv- MP.KB) <= tolPhys ...
        & abs(KGv- MP.KG) <= tolPhys;

    if any(hit)
        ix = find(hit,1,'last');     % newest
        m  = entries{ix}.meta;
        label = field_or_default(m,'label',"");
        E     = field_or_default(m,'E',NaN);
        P     = field_or_default(m,'P',NaN);
        solved = true;
    end
end

function v = g(e,a,b)
    if isstruct(e) && isfield(e,a) && isstruct(e.(a)) && isfield(e.(a),b)
        v = double(e.(a).(b));
    else
        v = NaN;
    end
end

function v = field_or_default(s,k,def)
    if isstruct(s) && isfield(s,k), v = s.(k); else, v = def; end
end

% ---------------- cell geometry ----------------
function [C1,C2,C3,C4] = subdivideCell(C)
    a1 = C.bounds(1,1); a2 = C.bounds(1,2);
    b1 = C.bounds(2,1); b2 = C.bounds(2,2);
    am = 0.5*(a1 + a2);
    bm = 0.5*(b1 + b2);
    d1 = C.depth + 1;

    C1 = makeCell([a1, am; b1, bm], d1);  % SW
    C2 = makeCell([am, a2; b1, bm], d1);  % SE
    C3 = makeCell([am, a2; bm, b2], d1);  % NE
    C4 = makeCell([a1, am; bm, b2], d1);  % NW
end

function C = makeCell(bounds, depth)
    if nargin < 2, depth = 0; end
    H1 = bounds(1,:); H2 = bounds(2,:);
    SW = [H1(1) H2(1)]; SE = [H1(2) H2(1)];
    NE = [H1(2) H2(2)]; NW = [H1(1) H2(2)];

    C = struct( ...
        'depth',depth, ...
        'bounds',bounds, ...
        'corners', [SW; SE; NE; NW], ...
        'cornerSolved', false(4,1), ...
        'cornerKey', strings(4,1), ...
        'cornerLabel', strings(4,1), ...
        'cornerEnergy', NaN(4,1), ...
        'cornerPressure', NaN(4,1), ...
        'isUniform', false, ...
        'mixedEdges', zeros(0,2) ...
    );
end

% ---------------- uniform test ----------------
function [uniform, mixedEdges] = uniformTest(C, eTol, pTol, tau)
% Uniform if all corners solved and either:
%   (i) labels all equal, OR
%   (ii) energy/pressure spreads small (and optional shape spread small)

    if any(~C.cornerSolved)
        uniform = false; mixedEdges = zeros(0,2); return;
    end

    labs = string(C.cornerLabel);
    if all(labs == labs(1))
        uniform = true; mixedEdges = zeros(0,2); return;
    end

    E = C.cornerEnergy;   E(~isfinite(E)) = inf;
    P = C.cornerPressure; P(~isfinite(P)) = inf;

    eSpread = max(E) - min(E);
    pSpread = max(P) - min(P);

    shapeSmall = true;
    if isfield(C,'cornerRmin') && all(isfinite(C.cornerRmin))
        r = C.cornerRmin;
        shapeSmall = (max(r)-min(r)) <= tau;
    end

    uniform = (eSpread <= eTol) && (pSpread <= pTol) && shapeSmall;

    edges = [1 2; 2 3; 3 4; 4 1];
    mixedEdges = zeros(0,2);
    if ~uniform
        for k = 1:4
            if labs(edges(k,1)) ~= labs(edges(k,2))
                mixedEdges(end+1,:) = edges(k,:); %#ok<AGROW>
            end
        end
    end
end
