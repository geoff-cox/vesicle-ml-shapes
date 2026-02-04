function s = sql_lit(x)
%SQL_LIT Convert a MATLAB scalar to an SQLite literal.
% Supports: numeric/logical scalar, string/char, empty -> NULL.

    if isempty(x)
        s = "NULL"; return;
    end

    if islogical(x)
        s = string(double(x)); return;
    end

    if isnumeric(x)
        if ~isscalar(x)
            error('sql_lit:nonscalar','Numeric SQL literals must be scalars.');
        end
        if ~isfinite(x)
            s = "NULL";
        else
            s = string(sprintf('%.17g', double(x)));
        end
        return;
    end

    if ischar(x) || isstring(x)
        s = sql_quote(x); return;
    end

    if isdatetime(x)
        s = sql_quote(iso_utc(x)); return;
    end

    error('sql_lit:unsupported','Unsupported type for sql_lit: %s', class(x));
end
