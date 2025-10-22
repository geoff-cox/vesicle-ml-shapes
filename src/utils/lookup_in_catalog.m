function [solved,label,E,P] = lookup_in_catalog(T,H1,H2,MP,tol)
    solved=false; label=""; E=NaN; P=NaN;
    if isempty(T) || ~any(T.Properties.VariableNames=="entry"), return; end
    H1v = cellfun(@(e) g(e,'params','H0_1'), T.entry);
    H2v = cellfun(@(e) g(e,'params','H0_2'), T.entry);
    Av  = cellfun(@(e) g(e,'params','A'),    T.entry);
    Vv  = cellfun(@(e) g(e,'params','V'),    T.entry);
    KAv = cellfun(@(e) g(e,'params','KA'),   T.entry);
    KBv = cellfun(@(e) g(e,'params','KB'),   T.entry);
    KGv = cellfun(@(e) g(e,'params','KG'),   T.entry);

    hit = abs(H1v-H1)<=tol & abs(H2v-H2)<=tol & Av==MP.A & Vv==MP.V ...
        & KAv==MP.KA & KBv==MP.KB & KGv==MP.KG;

    if any(hit)
        ix = find(hit,1,'last');
        m  = T.entry{ix}.meta;
        label = field_or_default(m,'label',"");
        E     = field_or_default(m,'E',NaN);
        P     = field_or_default(m,'P',NaN);
        solved = true;
    end
end

function v = g(e,a,b)
    if isstruct(e) && isfield(e,a) && isstruct(e.(a)) && isfield(e.(a),b)
        v = double(e.(a).(b)); else, v = NaN; end
end

function v = field_or_default(s,k,def)
    if isstruct(s) && isfield(s,k), v = s.(k); else, v = def; end
end