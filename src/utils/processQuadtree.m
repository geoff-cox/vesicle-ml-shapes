% === FILE: processQuadtree.m ===
% === patch ===
% Only showing the critical initialization lines to change.
% =============
function [task, cache] = processQuadtree(cache, T, MP)

    if nargin < 3, MP = struct(); end
    task = [];
    cache = ensure_defaults(cache);

    if ~isfield(cache,'QT') || isempty(cache.QT) || isempty(cache.QT.queue)
        if ~isfield(cache,'config') || ~isfield(cache.config,'bounds')
            cache.config.bounds = [-1 1; -1 1];
        end
        C0 = makeCell(cache.config.bounds);
        cache.QT.queue = {C0};
        cache.QT.cells = C0([]);   % <-- FIX: empty struct array, not []
    end

    tol = 1e-12;

    while ~isempty(cache.QT.queue)

        C = cache.QT.queue{1};
        cache.QT.queue(1) = [];

        [C, anyUnknown] = refresh_corners_from_catalog(C, T, MP, tol);

        if anyUnknown
            k = find(~C.cornerSolved,1,'first');

            % Optional: skip blocked corners if cache.failures/blockUntil exists
            if isfield(cache,'failures') && isfield(cache,'iter') && ~isempty(cache.failures)
                k = pick_unblocked_corner(C, cache);
                if isempty(k)
                    cache.QT.queue{end+1} = C; % defer
                    continue
                end
            end

            params = struct('H0_1',C.corners(k,1),'H0_2',C.corners(k,2));
            cache.QT.queue = [{C}, cache.QT.queue];
            task = struct('params',params);
            return
        end

        [uniform,mixedEdges] = uniformTest(C, cache.config.eTol, cache.config.pTol, cache.config.shapeTau);
        C.isUniform = uniform;
        C.mixedEdges = mixedEdges;

        if uniform || C.depth >= cache.config.maxDepth || numel(cache.QT.cells) >= cache.config.maxCells
            cache.QT.cells = [cache.QT.cells; C];
        else
            [C1,C2,C3,C4] = subdivideCell(C);
            cache.QT.queue(end+1:end+4) = {C1,C2,C3,C4};
        end
    end
end

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
    for i=1:numel(F)
        if abs(F(i).H0_1-h1)<tol && abs(F(i).H0_2-h2)<tol
            if isfield(F(i),'blockUntil') && F(i).blockUntil > it
                tf = true;
            end
            return
        end
    end
end

function cache = ensure_defaults(cache)
    if ~isfield(cache,'config'), cache.config = struct(); end
    if ~isfield(cache.config,'bounds')
        % Default box if not provided: [-1, +1] x [-1, +1]
        cache.config.bounds = [-1 1; -1 1];
    end
    cache.config.maxDepth = defaultArg(cache.config,'maxDepth', 6);
    cache.config.maxCells = defaultArg(cache.config,'maxCells', 2000);
    cache.config.eTol     = defaultArg(cache.config,'eTol',     5e-3);
    cache.config.pTol     = defaultArg(cache.config,'pTol',     5e-3);
    cache.config.shapeTau = defaultArg(cache.config,'shapeTau', 0.08);
end

function [C, anyUnknown] = refresh_corners_from_catalog(C, T, MP, tol)
    % For each corner, check if a matching param exists in catalog T.
    % If yes, fill label/energy/pressure and mark solved.
    % A simple exact match on (H0_1,H0_2) is used; you can add tolerances if needed.
    anyUnknown = false;
    for i=1:4
        H1 = C.corners(i,1); H2 = C.corners(i,2);
        [solved,label,E,P] = lookup_in_catalog(T,H1,H2,MP,tol);
        C.cornerSolved(i)=solved; C.cornerLabel(i)=label;
        C.cornerEnergy(i)=E; C.cornerPressure(i)=P;
        if ~solved, anyUnknown=true; end
    end

    function [solved,label,E,P] = lookup_in_catalog(T,H1,H2,MP,tol)
        % LOOKUP_IN_CATALOG  Find a solved corner in catalog by (H0_1,H0_2) and physics.
        % Robust to T.entry being either a cell column OR a struct array.

        solved=false; label=""; E=NaN; P=NaN;
        if isempty(T) || ~any(T.Properties.VariableNames=="entry"), return; end

        entries = T.entry;
        if isstruct(entries)
            entries = num2cell(entries);   % <-- critical: support struct-array catalogs
        elseif ~iscell(entries)
            try
                entries = num2cell(entries);
            catch
                return;
            end
        end

        H1v = cellfun(@(e) g(e,'params','H0_1'), entries);
        H2v = cellfun(@(e) g(e,'params','H0_2'), entries);

        % Physics params
        Av  = cellfun(@(e) g(e,'params','A'),  entries);
        Vv  = cellfun(@(e) g(e,'params','V'),  entries);
        KAv = cellfun(@(e) g(e,'params','KA'), entries);
        KBv = cellfun(@(e) g(e,'params','KB'), entries);
        KGv = cellfun(@(e) g(e,'params','KG'), entries);

        % --- tolerances ---
        tolH0 = tol;
        tolPhys = 1e-12;  % robust default for floating comparisons
        if isstruct(MP) && isfield(MP,'tolPhys') && isnumeric(MP.tolPhys) && MP.tolPhys > 0
            tolPhys = MP.tolPhys;
        end

        hit = abs(H1v-H1)<=tolH0 & abs(H2v-H2)<=tolH0 ...
            & abs(Av - MP.A)  <= tolPhys ...
            & abs(Vv - MP.V)  <= tolPhys ...
            & abs(KAv- MP.KA) <= tolPhys ...
            & abs(KBv- MP.KB) <= tolPhys ...
            & abs(KGv- MP.KG) <= tolPhys;

        if any(hit)
            ix = find(hit,1,'last');   % keep newest (consistent with your original)
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
    
end

function [C1,C2,C3,C4] = subdivideCell(C)
    % Split into 4 children (SW, SE, NE, NW), inherit depth+1.
    a1 = C.corners(1,1);
    a2 = C.corners(2,1);
    b1 = C.corners(2,2);
    b2 = C.corners(3,2);
    am = 0.5*(a1 + a2);
    bm = 0.5*(b1 + b2);
    d1 = C.depth + 1;

    C1 = makeCell([a1, am; b1, bm], d1);  % SW
    C2 = makeCell([am, a2; b1, bm], d1);  % SE
    C3 = makeCell([am, a2; bm, b2], d1);  % NE
    C4 = makeCell([a1, am; bm, b2], d1);  % NW

    function C = makeCell(bounds, depth)
        if nargin < 2
            depth = 0; % Default depth if not provided
        end
        % A quadtree cell covering [a1,a2] x [b1,b2].
        H1 = bounds(1,:); H2 = bounds(2,:);
        SW = [H1(1) H2(1)]; SE = [H1(2) H2(1)];
        NE = [H1(2) H2(2)]; NW = [H1(1) H2(2)];
        C = struct( ...
            'depth',depth, ...
            'corners', [SW; SE; NE; NW], ...
            'cornerSolved', false(4,1), ...
            'cornerKey', strings(4,1), ...
            'cornerLabel', nan(4,1), ...
            'cornerEnergy', nan(4,1), ...
            'cornerPressure', nan(4,1), ...
            'isUniform', false, ...
            'mixedEdges', zeros(0,2) ... % pairs of corner indices with differing labels
        );
    end
end

% === FILE: uniformTest.m (replacement) ===
function [uniform, mixedEdges] = uniformTest(C, eTol, pTol, tau)
% Uniform if:
%   - all corners solved, AND
%   - either all labels equal OR energy/pressure spreads are below thresholds
% mixedEdges returned as index pairs corresponding to physical edges, inferred from coords.

    if isfield(C,'cornerSolved') && any(~C.cornerSolved)
        uniform = false; mixedEdges = zeros(0,2); return;
    end

    labs = string(C.cornerLabel);   % robust to string/cellstr/numeric-ish
    if all(labs == labs(1))
        uniform = true; mixedEdges = zeros(0,2); return;
    end

    E = C.cornerEnergy;   E(~isfinite(E)) = inf;
    P = C.cornerPressure; P(~isfinite(P)) = inf;

    eSpread = max(E) - min(E);
    pSpread = max(P) - min(P);

    % optional descriptor hook (e.g., rMinAway) if you later add it into C:
    shapeSmall = true; %#ok<NASGU>
    if isfield(C,'cornerRmin') && all(isfinite(C.cornerRmin))
        r = C.cornerRmin;
        shapeSmall = (max(r)-min(r)) <= tau;
    end

    uniform = (eSpread <= eTol) && (pSpread <= pTol) && shapeSmall;

    % infer LL/LR/UR/UL indices from coordinates
    idx = inferCornerOrder(C.corners);  % [LL LR UR UL]
    edges = [idx(1) idx(2);  idx(2) idx(3);  idx(3) idx(4);  idx(4) idx(1)];

    mixedEdges = zeros(0,2);
    if ~uniform
        for k = 1:4
            if labs(edges(k,1)) ~= labs(edges(k,2))
                mixedEdges(end+1,:) = edges(k,:); %#ok<AGROW>
            end
        end
    end

    function idx = inferCornerOrder(X)
    % X is 4x2: [H0_1 H0_2]
        x = X(:,1); y = X(:,2);

        % LL: min x and min y
        [~,iLL] = min(x + 10*y);  % weight y to break ties robustly
        % LR: max x and min y
        [~,iLR] = max(x - 10*y);
        % UR: max x and max y
        [~,iUR] = max(x + 10*y);
        % UL: min x and max y
        [~,iUL] = min(x - 10*y);

        idx = [iLL iLR iUR iUL];
        if numel(unique(idx)) < 4
            error('uniformTest: could not infer unique corner ordering from coordinates.');
        end
    end

end


