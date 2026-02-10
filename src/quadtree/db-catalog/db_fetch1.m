function v = db_fetch1(db, sql)
%DB_FETCH1 Fetch first cell of first row (or []).
    rows = fetch(db, string(sql));
    if isempty(rows), v = []; return; end
    if iscell(rows)
        v = rows{1,1};
    else
        v = rows(1,1);
    end
end
