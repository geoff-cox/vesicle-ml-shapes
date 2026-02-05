% INSPECT_SOLUTION_POLES  Diagnostic tool to examine pole behavior in solved vesicles
%
% This script loads solved vesicle shapes and checks:
%   1. Whether north pole r returns to axis (r→0)
%   2. Pole regularity conditions
%   3. Boundary condition residuals at poles
%   4. Energy and volume integration accuracy
%
% Usage:
%   inspect_solution_poles()  % Uses catalog to load recent solutions
%   inspect_solution_poles(hash)  % Inspect specific solution by hash
%   inspect_solution_poles('random', N)  % Inspect N random solutions

function inspect_solution_poles(varargin)
    
    % Parse inputs
    if nargin == 0
        mode = 'recent';
        n_samples = 5;
    elseif nargin == 1
        mode = 'specific';
        hash_id = varargin{1};
    elseif nargin == 2
        mode = varargin{1};
        n_samples = varargin{2};
    end
    
    fprintf('\n=== Solution Pole Inspection ===\n\n');
    
    % Load catalog
    catalog_path = 'sim-results';
    if ~exist(fullfile(catalog_path, 'catalog.mat'), 'file')
        fprintf('Error: catalog.mat not found. Run simulation first.\n');
        return;
    end
    
    T = catalog_load(catalog_path);
    fprintf('Loaded catalog with %d entries\n', height(T));
    
    % Select solutions to inspect
    if strcmp(mode, 'recent')
        fprintf('Inspecting %d most recent solutions...\n\n', n_samples);
        [~, sort_idx] = sort(T.timestamp, 'descend');
        inspect_idx = sort_idx(1:min(n_samples, height(T)));
    elseif strcmp(mode, 'random')
        fprintf('Inspecting %d random solutions...\n\n', n_samples);
        inspect_idx = randperm(height(T), min(n_samples, height(T)));
    elseif strcmp(mode, 'specific')
        fprintf('Inspecting specific solution: %s\n\n', hash_id);
        inspect_idx = find(strcmp(T.hash, hash_id));
        if isempty(inspect_idx)
            fprintf('Error: Hash not found in catalog\n');
            return;
        end
    end
    
    % Inspect each solution
    for i = 1:length(inspect_idx)
        idx = inspect_idx(i);
        entry = T.entry{idx};
        hash = T.hash{idx};
        
        fprintf('--- Solution %d/%d: %s ---\n', i, length(inspect_idx), hash(1:12));
        fprintf('  Parameters: H0_1=%.3f, H0_2=%.3f\n', ...
            entry.params.H0_1, entry.params.H0_2);
        
        % Load solution
        sol_file = fullfile(catalog_path, 'hashed_results', [hash '.mat']);
        if ~exist(sol_file, 'file')
            fprintf('  ⚠️  Solution file not found: %s\n\n', sol_file);
            continue;
        end
        
        data = load(sol_file);
        sol = data.result.sol;
        
        % Extract state at poles
        y_south = sol.y(:, 1);    % s = 0 (south pole)
        y_north = sol.y(:, end);  % s = π (north pole)
        
        % Parse state vectors
        [QAs, HAs, PAs, rAs, zAs, ~, sAs, VAs, EAs] = parse_state(y_south(1:9));
        [QBn, HBn, PBn, rBn, zBn, ~, sBn, VBn, EBn] = parse_state(y_north(10:18));
        
        % Report south pole
        fprintf('\n  South Pole (α-phase, s=0):\n');
        fprintf('    Q = %.6f  (should be ~0)\n', QAs);
        fprintf('    H = %.6f\n', HAs);
        fprintf('    P = %.6f  (should be ~0)\n', PAs);
        fprintf('    r = %.6f  (should be ~0)\n', rAs);
        fprintf('    z = %.6f  (should be ~0)\n', zAs);
        fprintf('    s = %.6f  (should be ~0)\n', sAs);
        fprintf('    V = %.6f  (should be ~0)\n', VAs);
        fprintf('    E = %.6f  (should be ~0)\n', EAs);
        
        % Report north pole
        fprintf('\n  North Pole (β-phase, s=π):\n');
        fprintf('    Q = %.6f  (should be ~0)\n', QBn);
        fprintf('    H = %.6f\n', HBn);
        fprintf('    P = %.6f  (should be ~π = %.6f)\n', PBn, pi);
        fprintf('    r = %.6f  ❓ IS THIS ZERO?\n', rBn);
        fprintf('    z = %.6f\n', zBn);
        fprintf('    s = %.6f  (should be ~0)\n', sBn);
        fprintf('    V = %.6f  (should be ~0)\n', VBn);
        fprintf('    E = %.6f  (should be ~0)\n', EBn);
        
        % Check critical values
        south_ok = (abs(QAs) < 1e-5) && (abs(PAs) < 1e-5) && ...
                   (abs(rAs) < 1e-5) && (abs(zAs) < 1e-5);
        north_ok = (abs(QBn) < 1e-5) && (abs(PBn - pi) < 1e-5);
        north_axis = (abs(rBn) < 1e-3);  % Generous threshold
        
        fprintf('\n  Validation:\n');
        if south_ok
            fprintf('    ✅ South pole: All BC satisfied\n');
        else
            fprintf('    ❌ South pole: BC violations detected\n');
        end
        
        if north_ok
            fprintf('    ✅ North pole: Q and P conditions satisfied\n');
        else
            fprintf('    ❌ North pole: Q or P conditions violated\n');
        end
        
        if north_axis
            fprintf('    ✅ North pole: r ≈ 0 (returns to axis)\n');
        else
            fprintf('    ⚠️  North pole: r = %.6f (NOT on axis)\n', rBn);
            fprintf('        This suggests neck/tubule geometry\n');
        end
        
        % Energy and volume totals
        fprintf('\n  Integrated Quantities:\n');
        fprintf('    Total Volume (α): %.6f\n', y_south(8));
        fprintf('    Total Volume (β): %.6f\n', y_north(17));
        fprintf('    Total Energy (α): %.6f\n', y_south(9));
        fprintf('    Total Energy (β): %.6f\n', y_north(18));
        
        fprintf('\n');
    end
    
    fprintf('=== Inspection Complete ===\n\n');
end

function [Q, H, P, r, z, L, s, V, E] = parse_state(y)
    Q = y(1);
    H = y(2);
    P = y(3);
    r = y(4);
    z = y(5);
    L = y(6);
    s = y(7);
    V = y(8);
    E = y(9);
end
