% --- plot_residual_cloud.m
run('bootstrap.m');

d = dir(fullfile('SimResults','**','catalog.mat')); [~,i]=max([d.datenum]); L = load(fullfile(d(i).folder,d(i).name));
T = L.T;

figure('Color','w'); loglog(T.BCmax, T.DEmax, '.'); grid on;
xlabel('BC residual (max)'); ylabel('DE residual (max)');
title('Acceptance cloud'); xline(1e-6,'--'); yline(2e-1,'--');
