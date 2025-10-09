% --- verify_residuals.m
run('bootstrap.m');

% find most recent run:
d = dir(fullfile('SimResults','run_*')); if isempty(d), d = dir('SimResults'); end
[~,i] = max([d.datenum]); runDir = fullfile(d(i).folder, d(i).name);

L = load(fullfile(runDir,'catalog.mat')); T = L.T;

% Basic acceptance gates (tweak as needed)
BC_OK = T.BCmax <= 1e-6;
DE_OK = T.DEmax <= 2e-1;

fprintf('BC pass: %d/%d\n', sum(BC_OK), height(T));
fprintf('DE pass: %d/%d\n', sum(DE_OK), height(T));

% Quick scatter: DE vs BC colored by label
figure('Color','w'); 
loglog(T.BCmax, T.DEmax, '.'); grid on;
xlabel('BC residual (max)'); ylabel('DE residual (max)');
title('Gate diagnostics (all accepted points)');
