% --- verify_cache_hits.m
run('bootstrap.m');

% Run the driver twice (short mode) and confirm disk-cache hits on second pass
sim_driver_quad_tree;  % pass 1
sim_driver_quad_tree;  % pass 2 (should log "... HIT disk cache ..." for many points)

% Grep OPfile.txt (or just open) to verify cache paths and acceptance lines.
type(fullfile('SimResults', dir('SimResults').name, 'OPfile.txt'));
