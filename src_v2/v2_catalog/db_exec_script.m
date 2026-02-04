function db_exec_script(db, scriptText)
%DB_EXEC_SCRIPT Execute a SQL script containing multiple statements.

    stmts = split_sql_statements(scriptText);
    for i = 1:numel(stmts)
        s = strtrim(stmts{i});
        if strlength(s)==0, continue; end
        exec(db, string(s));
    end
end

function parts = split_sql_statements(txt)
% Split on semicolons not inside single quotes.
    c = char(string(txt));
    parts = {};
    buf = char.empty(1,0);
    inQuote = false;

    k = 1;
    while k <= numel(c)
        ch = c(k);

        if ch == ''''
            if inQuote
                % escaped quote: ''
                if k < numel(c) && c(k+1) == ''''
                    buf(end+1:end+2) = ''''''; %#ok<AGROW>
                    k = k + 2;
                    continue;
                else
                    inQuote = false;
                end
            else
                inQuote = true;
            end
        end

        if ch == ';' && ~inQuote
            parts{end+1} = string(buf); %#ok<AGROW>
            buf = char.empty(1,0);
        else
            buf(end+1) = ch; %#ok<AGROW>
        end
        k = k + 1;
    end

    if ~isempty(strtrim(buf))
        parts{end+1} = string(buf); %#ok<AGROW>
    end
end