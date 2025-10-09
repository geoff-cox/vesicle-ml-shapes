% --- check_interface_curvature.m
run('bootstrap.m');

% load a recent accepted solution
resDir = fullfile('SimResults'); d = dir(fullfile(resDir,'**','results','*.mat'));
S = load(fullfile(d(end).folder,d(end).name),'sol');

% pull out neck indices (last node of alpha equals first of beta in your split setup)
psiA = S.sol.y(3,:);  rA = S.sol.y(4,:);  HA = S.sol.y(2,:);
psiB = S.sol.y(12,:); rB = S.sol.y(13,:); HB = S.sol.y(11,:);

% mean curvature 2H = psi' + sin(psi)/r (here, H stored directly in y(2,:), y(11,:))
H_jump = 2*HB(end) - 2*HA(end);
fprintf('Interface mean-curvature jump (2H_B - 2H_A) ~ %.4e\n', H_jump);
