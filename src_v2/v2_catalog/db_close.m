function db_close(db)
%DB_CLOSE Close SQLite handle safely.
    if isempty(db), return; end
    try
        close(db);
    catch
    end
end
