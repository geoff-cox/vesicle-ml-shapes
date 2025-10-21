% -------------------------------------------------------------------------
% EXTRACTED HELPER for "normalizeForHash"
%   - Source: sim_driver_quad_tree_full.m
%   - Extracted: 2025-10-11 11:50:12
%   - Sub-helpers (nested functions) are retained in this file.
% -------------------------------------------------------------------------

function y = normalizeForHash(x)
    if istable(x)
        x = table2struct(x, 'ToScalar', true);
    end
    if isstruct(x)
        if numel(x) > 1
            y = arrayfun(@normalizeForHash,x,'UniformOutput',false); return;
        end
        f = sort(fieldnames(x));
        s = struct();
        for i=1:numel(f), s.(f{i}) = normalizeForHash(x.(f{i})); end
        y = s; return
    elseif iscell(x)
        y = cellfun(@normalizeForHash,x,'UniformOutput',false); return
    elseif isstring(x)
        y = cellstr(x); return
    elseif isa(x,'categorical')
        y = struct('codes',double(x),'cats',{categories(x)},'ord',isordinal(x)); return
    else
        y = x;
    end
end