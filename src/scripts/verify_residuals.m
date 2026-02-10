disp('lick it!')
% runDir = 'SimResults';
% d = dir(runDir);
% latest_date = datetime("01-Jan-0000");
% for k = 1:length(d)
%     try
%         run_date = datetime(d(k).name);
%         if run_date > latest_date
%             latest_date = run_date; % Update latest_date if a more recent date is found
%             runDir = fullfile(runDir, string(latest_date));
%         end
%     catch
%         continue; % Skip to the next iteration if there's an error in date conversion
%     end
% end
% 
% runDir = fullfile(runDir, "run_*");
% d = dir(runDir); 
% if isempty(d)
%     d = dir(runDir);
% end
% 
% % find most recent run:
% [~,i] = max([d.datenum]);
% runDir = fullfile(d(i).folder, d(i).name);
% 
% L = load(fullfile(runDir,'catalog.mat')); T = L.T;
% 
% % Basic acceptance gates (tweak as needed)
% BC_OK = T.BCmax <= 1e-6;
% DE_OK = T.DEmax <= 2e-1;
% 
% fprintf('BC pass: %d/%d\n', sum(BC_OK), height(T));
% fprintf('DE pass: %d/%d\n', sum(DE_OK), height(T));
% 
% % Quick scatter: DE vs BC colored by label
% figure('Color','w'); 
% loglog(T.BCmax, T.DEmax, '.'); grid on;
% xlabel('BC residual (max)'); ylabel('DE residual (max)');
% title('Gate diagnostics (all accepted points)');