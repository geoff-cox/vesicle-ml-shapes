run('bootstrap.m');

% Pick one solved hash file in SimResults/<DATE>/results/*.mat
resDir = fullfile('SimResults'); d = dir(fullfile(resDir,'**','results','*.mat'));
assert(~isempty(d), 'No hashed results found yet.');
picks = 1:2:length(d); %[3 5 7 10 14 22];
for pick = picks
    data = fullfile(d(pick).folder, d(pick).name);   % newest
    S = load(data,'sol','E_total','P_osm');
    
    % Variables (by convention in your code)
    s   = S.sol.x;
    y   = S.sol.y;
    r_A = y(4,:); z_A = y(5,:);
    r_B = y(13,:); z_B = y(14,:);
    
    % Plot the 2D profile (both phases)
    figure('Color','w'); hold on; axis equal; grid on;
    plot(r_A, z_A, '-', 'LineWidth',1.5);
    plot(r_B, z_B, '-', 'LineWidth',1.5);
    line([0 0], ylim, 'Color',[0.5 0.5 0.5], 'LineStyle','--'); % axis of rotation
    xlabel('r'); ylabel('z'); title(sprintf('Shape preview | E=%.4g | P=%.4g',S.E_total,S.P_osm));
    legend({'Phase (1)','Phase (2)'},'Location','best');
end

run('cleanup.m');