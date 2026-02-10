function s = sql_quote(x)
%SQL_QUOTE Quote a TEXT literal for SQLite (or NULL).
    if isempty(x)
        s = "NULL";
        return;
    end
    t = char(string(x));
    t = strrep(t, '''', ''''''); % escape single quotes
    s = "'" + string(t) + "'";
end
