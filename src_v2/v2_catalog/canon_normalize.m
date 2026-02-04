function y = canon_normalize(x)
%CANON_NORMALIZE Normalize MATLAB data into JSON-safe, stable ordering.
% - Struct fields sorted recursively
% - Containers.Map -> struct
% - datetime -> ISO-8601 UTC string
% - function_handle -> func2str
% - non-finite numbers -> strings "NaN"/"Inf"/"-Inf"
% - objects -> struct with class + public properties (best-effort)

    if isstruct(x)
        f = sort(fieldnames(x));
        y = struct();
        for i = 1:numel(f)
            y.(f{i}) = canon_normalize(x.(f{i}));
        end
        y = orderfields(y, f);
        return;
    end

    if isa(x, 'containers.Map')
        k = sort(string(x.keys));
        y = struct();
        for i=1:numel(k)
            key = char(k(i));
            y.(matlab.lang.makeValidName(key)) = canon_normalize(x(key));
        end
        y = orderfields(y);
        return;
    end

    if iscell(x)
        y = cellfun(@canon_normalize, x, 'UniformOutput', false);
        return;
    end

    if isstring(x) || ischar(x)
        y = string(x);
        return;
    end

    if isdatetime(x)
        y = iso_utc(x);
        return;
    end

    if isa(x,'function_handle')
        y = "func:" + string(func2str(x));
        return;
    end

    if isnumeric(x) || islogical(x)
        if islogical(x)
            y = x; return;
        end
        if any(~isfinite(x(:)))
            % Replace non-finite elements with strings (JSON-safe)
            y = arrayfun(@nf, double(x), 'UniformOutput', false);
            return;
        end
        y = x;
        return;
    end

    % Best-effort for objects
    try
        cls = class(x);
        props = properties(x);
        s = struct();
        s.__class__ = cls;
        for i=1:numel(props)
            p = props{i};
            try
                s.(p) = canon_normalize(x.(p));
            catch
                s.(p) = "<unreadable>";
            end
        end
        y = orderfields(s);
    catch
        y = "<unsupported:" + string(class(x)) + ">";
    end
end

function out = nf(v)
    if isnan(v), out = "NaN";
    elseif isinf(v) && v>0, out = "Inf";
    elseif isinf(v) && v<0, out = "-Inf";
    else, out = v;
    end
end
