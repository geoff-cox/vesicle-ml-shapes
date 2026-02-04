function id = db_last_insert_rowid(db)
    id = fetch(db, "SELECT last_insert_rowid();");
    id = id{1};
end
