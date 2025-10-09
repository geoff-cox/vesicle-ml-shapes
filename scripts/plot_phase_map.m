% --- plot_phase_map.m
run('bootstrap.m');

% Load newest catalog.mat (table T with H0_1,H0_2,label)
d = dir(fullfile('SimResults','**','catalog.mat')); [~,i]=max([d.datenum]); L = load(fullfile(d(i).folder,d(i).name));
T = L.T;

figure('Color','w'); hold on; box on; grid on;
gscatter(T.H0_1, T.H0_2, T.label, [], '.', 14);
xlabel('H_0^{(1)}'); ylabel('H_0^{(2)}');
title('Morphology map (quadtree samples)'); axis tight;
