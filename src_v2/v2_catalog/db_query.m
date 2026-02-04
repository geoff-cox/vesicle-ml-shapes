function rows = db_query(db, sql)
%DB_QUERY Fetch rows from a query. Returns cell array.
    rows = fetch(db, string(sql));
end
