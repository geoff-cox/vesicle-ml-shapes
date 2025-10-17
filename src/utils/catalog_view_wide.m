% File: src/utils/catalog_view_wide.m
% Create the 15-column table view from the 3-col catalog

function W = catalog_view_wide(simDir)
    T = catalog_load(simDir);
    n = height(T);
    W = table('Size',[n 15], ...
        'VariableTypes', {'datetime','string','double','double','double','double','double','double','double','string','double','double','double','double','double'}, ...
        'VariableNames', {'timestamp','hash','H0_1','H0_2','A','V','KA','KB','KG','label','E','P','BCmax','DEmax','mesh'});

    for i=1:n
        W.timestamp(i) = T.timestamp(i);
        W.hash(i)      = string(T.hash(i));
        e = T.entry{i};
        p = struct(); m = struct();
        if isfield(e,'params'), p = e.params; end
        if isfield(e,'meta'),   m = e.meta;   end

        W.H0_1(i) = getn(p,'H0_1');
        W.H0_2(i) = getn(p,'H0_2');
        W.A(i)    = getn(p,'A');
        W.V(i)    = getn(p,'V');
        W.KA(i)   = getn(p,'KA');
        W.KB(i)   = getn(p,'KB');
        W.KG(i)   = getn(p,'KG');

        W.label(i)  = gets(m,'label');
        W.E(i)      = getn(m,'E');
        W.P(i)      = getn(m,'P');
        W.BCmax(i)  = getn(m,'BCmax');
        W.DEmax(i)  = getn(m,'DEmax');
        W.mesh(i)   = getn(m,'mesh');
    end
end

function v = getn(S, k)
    if isstruct(S) && isfield(S,k) && isnumeric(S.(k)) && ~isempty(S.(k))
        v = double(S.(k));
    else
        v = NaN;
    end
end
function v = gets(S, k)
    if isstruct(S) && isfield(S,k) && ~isempty(S.(k))
        v = string(S.(k));
    else
        v = "";
    end
end
