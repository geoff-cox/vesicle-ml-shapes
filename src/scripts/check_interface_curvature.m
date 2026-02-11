% check_interface_curvature.m
% Diagnostic: compute mean-curvature jump at the interface from a solution.
%
% Usage: run from the src/ directory (or add src/ to path first).

here = fileparts(mfilename('fullpath'));
srcRoot = fileparts(here);
addpath(genpath(srcRoot));

% load a recent accepted solution
projRoot = fileparts(srcRoot);
resDir   = fullfile(projRoot, 'sim-results', 'hashed_results');
d = dir(fullfile(resDir, '*.mat'));
assert(~isempty(d), 'No hashed results found yet.');

S = load(fullfile(d(end).folder, d(end).name), 'result');
sol = S.result.sol;

% pull out neck indices
HA = sol.y(2,:);   HB = sol.y(11,:);

% mean curvature 2H = psi' + sin(psi)/r (here, H stored directly in y(2,:), y(11,:))
H_jump = 2*HB(end) - 2*HA(end);
fprintf('Interface mean-curvature jump (2H_B - 2H_A) ~ %.4e\n', H_jump);