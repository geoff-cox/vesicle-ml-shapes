% plot_phase_map.m
% Visualization: scatter plot of solved points colored by morphology label.
%
% Usage: run from the src/ directory (or add src/ to path first).

here = fileparts(mfilename('fullpath'));
srcRoot = fileparts(here);
addpath(genpath(srcRoot));

projRoot = fileparts(srcRoot);
simDir   = fullfile(projRoot, 'sim-results');
T = catalog_load(simDir);

assert(height(T) > 0, 'Catalog is empty — run simulations first.');

% Extract H0_1, H0_2, label from the entry structs
H1  = cellfun(@(e) defaultArg(e.params,'H0_1',NaN), T.entry);
H2  = cellfun(@(e) defaultArg(e.params,'H0_2',NaN), T.entry);
lab = cellfun(@(e) defaultArg(e.meta,'label',0),     T.entry);

valid = isfinite(H1) & isfinite(H2);

figure('Color','w'); hold on; box on; grid on;
gscatter(H1(valid), H2(valid), lab(valid), [], '.', 14);
xlabel('H_0^{(1)}'); ylabel('H_0^{(2)}');
title('Morphology map (quadtree samples)'); axis tight;