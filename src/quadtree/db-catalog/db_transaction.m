function out = db_transaction(db, fcn)
%DB_TRANSACTION Execute fcn(db) inside BEGIN IMMEDIATE / COMMIT.
% Rolls back on error and rethrows.

    exec(db, "BEGIN IMMEDIATE;");
    try
        out = fcn(db);
        exec(db, "COMMIT;");
    catch ME
        try exec(db, "ROLLBACK;"); catch, end %#ok<CTCH>
        rethrow(ME);
    end
end
