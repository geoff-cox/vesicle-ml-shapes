function db = db_open(dbfile)
%DB_OPEN Open SQLite DB with recommended pragmas.
%   db = db_open(fullfile(simDir,'catalog.db'));

    db = sqlite(dbfile, "create");

    % Pragmas (good defaults for append-heavy workloads)
    exec(db, "PRAGMA journal_mode=WAL;");
    exec(db, "PRAGMA synchronous=NORMAL;");
    exec(db, "PRAGMA foreign_keys=ON;");
    exec(db, "PRAGMA temp_store=MEMORY;");
end
