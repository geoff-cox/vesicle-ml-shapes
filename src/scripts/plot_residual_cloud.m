% plot_residual_cloud.m
% Visualization: log-log scatter of BC vs DE residuals for accepted solutions.
%
% Usage: run from the src/ directory (or add src/ to path first).

here = fileparts(mfilename('fullpath'));
srcRoot = fileparts(here);
addpath(genpath(srcRoot));

projRoot = fileparts(srcRoot);
simDir   = fullfile(projRoot, 'sim-results');
T = catalog_load(simDir);

assert(height(T) > 0, 'Catalog is empty — run simulations first.');

BCmax = cellfun(@(e) defaultArg(e.meta,'BCmax',NaN), T.entry);
DEmax = cellfun(@(e) defaultArg(e.meta,'DEmax',NaN), T.entry);

valid = isfinite(BCmax) & isfinite(DEmax);

figure('Color','w');
loglog(BCmax(valid), DEmax(valid), '.'); grid on;
xlabel('BC residual (max)'); ylabel('DE residual (max)');
title('Acceptance cloud'); xline(1e-6,'--'); yline(2e-1,'--');