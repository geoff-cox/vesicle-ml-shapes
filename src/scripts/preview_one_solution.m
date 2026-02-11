% preview_one_solution.m
% Visualization: 2D profile of alpha and beta phase shapes.
%
% Usage: run from the src/ directory (or add src/ to path first).

here = fileparts(mfilename('fullpath'));
srcRoot = fileparts(here);
addpath(genpath(srcRoot));

projRoot = fileparts(srcRoot);
resDir   = fullfile(projRoot, 'SimResults', 'hashed_results');
d = dir(fullfile(resDir, '*.mat'));
assert(~isempty(d), 'No hashed results found yet.');

picks = 1:2:length(d);
for pick = picks
    data = fullfile(d(pick).folder, d(pick).name);
    S = load(data, 'result', 'meta');

    sol = S.result.sol;
    r_A = sol.y(4,:);  z_A = sol.y(5,:);
    r_B = sol.y(13,:); z_B = sol.y(14,:);

    E_total = defaultArg(S.meta, 'E', NaN);
    P_osm   = defaultArg(S.meta, 'P', NaN);

    figure('Color','w'); hold on; axis equal; grid on;
    plot(r_A, z_A, '-', 'LineWidth',1.5);
    plot(r_B, z_B, '-', 'LineWidth',1.5);
    line([0 0], ylim, 'Color',[0.5 0.5 0.5], 'LineStyle','--');
    xlabel('r'); ylabel('z');
    title(sprintf('Shape preview | E=%.4g | P=%.4g', E_total, P_osm));
    legend({'Phase (1)','Phase (2)'},'Location','best');
end