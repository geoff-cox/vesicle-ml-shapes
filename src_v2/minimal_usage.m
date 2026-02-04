% resultsRoot = fullfile(projRoot,'SimResults_v2');
% runStruct should be a plain struct (no giant solver objects if you can avoid it).

runStruct = struct();
runStruct.model_version = sim.SP.ModelVersion;
runStruct.bounds = sim.SP.H0Bounds;
runStruct.qt_max_depth = sim.SP.QTmaxDepth;
runStruct.qt_max_cells = sim.SP.QTmaxCells;
runStruct.SP = sim.SP;
runStruct.TH = sim.TH;
runStruct.MP = sim.MP;

[db, run_id] = catalog_init_v2(resultsRoot, runStruct);

% For each attempted solution:
P = struct('H0_1',h1,'H0_2',h2,'A',sim.MP.A,'V',sim.MP.V,'KA',sim.MP.KA,'KB',sim.MP.KB,'KG',sim.MP.KG);

meta = struct();
meta.status = "ok";
meta.label  = "whatever";
meta.energy_E = E;
meta.pressure_P = P0;
meta.bc_max = BC;
meta.de_max = DE;
meta.mesh_n = meshN;
meta.wall_time_s = tSolve;
meta.seed_hash = "seedHashOptional";
meta.is_best = 1;

solution_id = catalog_insert_solution(db, run_id, sim.SP.ModelVersion, P, meta);

% After writing HDF5 to artifacts/solution/<solution_id>.h5:
relpath = "artifacts/solution/" + solution_id + ".h5";
artifact_id = catalog_add_artifact(db, solution_id, "solution_h5", "h5", relpath, "", []);

% Finish run:
catalog_finish_run(db, run_id, "finished");
db_close(db);
