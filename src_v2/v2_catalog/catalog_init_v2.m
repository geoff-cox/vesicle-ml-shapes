function [db, run_id] = catalog_init_v2(resultsRoot, runStruct, ddlPath)
%CATALOG_INIT_V2 Create/open catalog.db, apply DDL if new, insert run row.
%
% resultsRoot: path to SimResults_v2 directory
% runStruct: struct with at least fields:
%   runStruct.model_version
%   runStruct.bounds = [h01min h01max; h02min h02max]
%   runStruct.qt_max_depth, runStruct.qt_max_cells
%   runStruct.MP (optional: A,V,KA,KB,KG)
%   runStruct.SP/TH/MP/etc (anything) -> sim_json
% ddlPath (optional): path to ddl.sql (defaults to resultsRoot/schema/ddl.sql)

    if nargin < 3 || isempty(ddlPath)
        ddlPath = fullfile(resultsRoot, 'schema', 'ddl.sql');
    end

    ensure_dirs(resultsRoot);

    dbfile = fullfile(resultsRoot, 'catalog.db');
    isNew = ~exist(dbfile,'file');

    db = db_open(dbfile);

    if isNew
        if ~exist(ddlPath,'file')
            error('catalog_init_v2:missingDDL','ddl.sql not found at: %s', ddlPath);
        end
        ddl = fileread(ddlPath);
        db_exec_script(db, ddl);
    end

    run_id = new_uuid();
    started_at = iso_utc(datetime('now','TimeZone','UTC'));

    mv = string(runStruct.model_version);
    b = runStruct.bounds;

    qtDepth = runStruct.qt_max_depth;
    qtCells = runStruct.qt_max_cells;

    % Optional physics constants for the run
    A  = getfield_default(runStruct,'A',  []);
    V  = getfield_default(runStruct,'V',  []);
    KA = getfield_default(runStruct,'KA', []);
    KB = getfield_default(runStruct,'KB', []);
    KG = getfield_default(runStruct,'KG', []);

    [sim_json, sim_sha] = canon_json(runStruct);

    sql = "INSERT INTO runs(run_id,started_at,status,model_version," + ...
          "h0_1_min,h0_1_max,h0_2_min,h0_2_max,qt_max_depth,qt_max_cells," + ...
          "A,V,KA,KB,KG,sim_json,sim_json_sha256) VALUES (" + ...
          sql_lit(run_id) + "," + sql_lit(started_at) + "," + sql_lit("running") + "," + sql_lit(mv) + "," + ...
          sql_lit(b(1,1)) + "," + sql_lit(b(1,2)) + "," + sql_lit(b(2,1)) + "," + sql_lit(b(2,2)) + "," + ...
          sql_lit(qtDepth) + "," + sql_lit(qtCells) + "," + ...
          sql_lit_or_null(A) + "," + sql_lit_or_null(V) + "," + ...
          sql_lit_or_null(KA) + "," + sql_lit_or_null(KB) + "," + sql_lit_or_null(KG) + "," + ...
          sql_lit(sim_json) + "," + sql_lit(sim_sha) + ");";

    db_exec(db, sql);
end

function ensure_dirs(root)
    if ~exist(root,'dir'), mkdir(root); end
    d = {fullfile(root,'artifacts','solution'), fullfile(root,'cache'), fullfile(root,'schema')};
    for i=1:numel(d)
        if ~exist(d{i},'dir'), mkdir(d{i}); end
    end
end

function v = getfield_default(S, name, def)
    if isstruct(S) && isfield(S,name), v = S.(name); else, v = def; end
end

function s = sql_lit_or_null(x)
    if isempty(x), s = "NULL"; else, s = sql_lit(x); end
end
