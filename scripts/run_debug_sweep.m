% --- run_debug_sweep.m
run('bootstrap.m');

% Open the driver and force debug mode (short run) at runtime:
sim_driver_quad_tree;   % uses sim.debug.short = true inside the driver

% This should:
% - create SimResults/<DATE>/run_<timestamp>/
% - emit OPfile.txt (log), catalog.mat (table) and a solutions/ dir
