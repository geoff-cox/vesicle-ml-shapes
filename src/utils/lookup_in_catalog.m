function [solved,label,E,P] = lookup_in_catalog(T, H1, H2, MP, tol)
    solved=false;
    label="";
    E=NaN;
    P=NaN;
    if isempty(T), return; end
    if any(T.Properties.VariableNames=="entry")
        H1v = cellfun(@(e) fget(e,'params','H0_1'), T.entry);
        H2v = cellfun(@(e) fget(e,'params','H0_2'), T.entry);
        hit = abs(H1v-H1)<=tol & abs(H2v-H2)<=tol;
        % optionally also require MP match if provided
        if ~isempty(fieldnames(MP))
            Av  = cellfun(@(e) fget(e,'params','A'),  T.entry);
            Vv  = cellfun(@(e) fget(e,'params','V'),  T.entry);
            KAv = cellfun(@(e) fget(e,'params','KA'), T.entry);
            KBv = cellfun(@(e) fget(e,'params','KB'), T.entry);
            KGv = cellfun(@(e) fget(e,'params','KG'), T.entry);
            hit = hit   & Av  == MP.A ...
                        & Vv  == MP.V ...
                        & KAv == MP.KA ...
                        & KBv == MP.KB ...
                        & KGv == MP.KG;
        end
        if any(hit)
            idx = find(hit,1,'last');
            m = T.entry{idx}.meta;
            [label,E,P] = meta_fields(m); solved=true;
        end
    else
        if all(ismember({'H0_1','H0_2'}, T.Properties.VariableNames))
            hit = abs(T.H0_1-H1)<=tol & abs(T.H0_2-H2)<=tol;
            if any(hit)
                idx = find(hit,1,'last');
                label=string(T.label(idx));
                E=T.E(idx);
                P=T.P(idx);
                solved=true;
            end
        end
    end
end

function v = fget(e,a,b)
    if isstruct(e)  && isfield(e,a) ...
                    && isstruct(e.(a)) ...
                    && isfield(e.(a),b)
        v = double(e.(a).(b));
    else, v = NaN;
    end
end