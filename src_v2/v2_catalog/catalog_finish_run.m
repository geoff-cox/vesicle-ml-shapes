function catalog_finish_run(db, run_id, status)
%CATALOG_FINISH_RUN Mark a run finished/aborted with timestamp.
    if nargin < 3, status = "finished"; end
    t = iso_utc(datetime('now','TimeZone','UTC'));
    sql = "UPDATE runs SET finished_at=" + sql_lit(t) + ", status=" + sql_lit(status) + ...
          " WHERE run_id=" + sql_lit(run_id) + ";";
    db_exec(db, sql);
end
