function db_exec(db, sql)
%DB_EXEC Execute a statement (string).
    exec(db, string(sql));
end
