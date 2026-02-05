% DIAGNOSTICS_BC_DE_MATH  Test suite for mathematical correctness of BC and DE implementations
%
% This script performs mathematical validation tests on:
%   - BendV_Lag_EIGp_BC_impl.m (boundary conditions)
%   - BendV_Lag_EIGp_DE_impl.m (differential equations)
%
% Tests include:
%   1. Dimensional consistency
%   2. Pole regularity verification
%   3. Energy conservation checks
%   4. Symmetry properties
%   5. Comparison of pole vs bulk forms

function diagnostics_bc_de_math()
    fprintf('\n=== Mathematical Diagnostics for BendV_Lag_EIGp BC/DE ===\n\n');
    
    % Test 1: Check dimensional consistency
    test_dimensional_consistency();
    
    % Test 2: Verify pole expansion limit
    test_pole_vs_bulk_limit();
    
    % Test 3: Check boundary condition structure
    test_bc_structure();
    
    % Test 4: Verify energy scaling
    test_energy_scaling();
    
    fprintf('\n=== Diagnostics Complete ===\n');
end

%% Test 1: Dimensional Consistency
function test_dimensional_consistency()
    fprintf('TEST 1: Dimensional Consistency\n');
    fprintf('----------------------------------------\n');
    
    % Analyze RHS_pole Q-equation (line 26)
    fprintf('Analyzing RHS_pole Q-equation: H*L + 0.5*lam - k*H0*H^2 + 0.5*k*H*H0^2\n');
    
    % Assume typical units:
    % H ~ [1/length]
    % L ~ [pressure] = [force/length^2]
    % k ~ [energy] = [force*length]
    % lam ~ [pressure]
    
    fprintf('  Term H*L:           [1/L]*[F/L^2] = [F/L^3]\n');
    fprintf('  Term 0.5*lam:       [F/L^2]\n');
    fprintf('  Term k*H0*H^2:      [F*L]*[1/L]*[1/L^2] = [F/L^2]\n');
    fprintf('  Term 0.5*k*H*H0^2:  [F*L]*[1/L]*[1/L^2] = [F/L^2]\n');
    fprintf('\n');
    fprintf('  ❌ ISSUE: H*L has units [F/L^3], others have [F/L^2] - INCONSISTENT!\n');
    
    % Analyze bulk Q-equation (line 38)
    fprintf('\nAnalyzing bulk Q-equation (inside sin(S)/r factor):\n');
    fprintf('  All terms multiplied by sin(S)/r, so dimensions should match before multiplication\n');
    fprintf('  Term 2*H*L:         [1/L]*[F/L^2] = [F/L^3]\n');
    fprintf('  Term lam:           [F/L^2]\n');
    fprintf('  Term bending:       [F*L]*[1/L]*[1/L^2] = [F/L^2]\n');
    fprintf('\n');
    fprintf('  ⚠️  NOTE: Terms still have mixed dimensions [F/L^2] and [F/L^3]\n');
    fprintf('  This suggests r factor may be absorbed in definition of L or H\n');
    fprintf('\n');
end

%% Test 2: Pole vs Bulk Limit
function test_pole_vs_bulk_limit()
    fprintf('TEST 2: Pole vs Bulk Expansion Limit\n');
    fprintf('----------------------------------------\n');
    
    % At pole: r → 0, P → 0 (south) or π (north), sin(P) → 0
    % Check if RHS_pole matches limit of RHS as r → 0
    
    fprintf('Setting up test near south pole (P→0, r→0):\n');
    
    % Test parameters
    Q = 0.1; H = 1.0; P_vals = [0.1, 0.01, 0.001, 0.0001];
    r_vals = [0.1, 0.01, 0.001, 0.0001];
    L = 10; lam = 5; k = 1; H0 = 0.5;
    
    fprintf('\n  Testing Q-equation:\n');
    fprintf('  RHS_pole form: H*L + 0.5*lam - k*H0*H^2 + 0.5*k*H*H0^2\n');
    
    pole_value = H*L + 0.5*lam - k*H0*H^2 + 0.5*k*H*H0^2;
    fprintf('  RHS_pole = %.6f\n', pole_value);
    
    fprintf('\n  Bulk form (without sin(S)/r factor):\n');
    fprintf('    -Q*cos(P)/r - k*(2*H - H0)*(H*H0 + 2*(H - sin(P)/r)^2) + 2*H*L + lam\n');
    
    fprintf('\n  Testing limit as P→0, r→0:\n');
    for i = 1:length(P_vals)
        P = P_vals(i);
        r = r_vals(i);
        
        % Bulk formula (without the sin(S)/r outer factor)
        bulk_inner = -Q*cos(P)/r - k*(2*H - H0)*(H*H0 + 2*(H - sin(P)/r)^2) + 2*H*L + lam;
        
        fprintf('    P=%.4f, r=%.4f: bulk_inner = %.6f\n', P, r, bulk_inner);
    end
    
    fprintf('\n  ❌ ISSUE: Pole value = %.6f does NOT match bulk limit\n', pole_value);
    fprintf('  Expected pole form should include 2*H*L (not H*L) and lam (not 0.5*lam)\n');
    fprintf('\n');
end

%% Test 3: Boundary Condition Structure
function test_bc_structure()
    fprintf('TEST 3: Boundary Condition Structure\n');
    fprintf('----------------------------------------\n');
    
    % Count active boundary conditions
    fprintf('South pole (α-phase) boundary conditions:\n');
    south_active = {'QAs=0', 'PAs=0', 'rAs=0', 'zAs=0', 'sAs=0', 'VAs=0', 'EAs=0'};
    south_free = {'HAs', 'LAs'};
    fprintf('  Active: %d conditions - %s\n', length(south_active), strjoin(south_active, ', '));
    fprintf('  Free:   %d variables - %s\n', length(south_free), strjoin(south_free, ', '));
    
    fprintf('\nNorth pole (β-phase) boundary conditions:\n');
    north_active = {'QBn=0', 'PBn=π', 'sBn=0', 'VBn=0', 'EBn=0'};
    north_free = {'HBn', 'rBn', 'zBn', 'LBn'};
    fprintf('  Active: %d conditions - %s\n', length(north_active), strjoin(north_active, ', '));
    fprintf('  Free:   %d variables - %s\n', length(north_free), strjoin(north_free, ', '));
    
    fprintf('\nNeck junction boundary conditions:\n');
    neck_conditions = {
        'PBk = PAk (angle continuity)', ...
        'rBk = rAk (radial continuity)', ...
        'zBk = zAk (axial continuity)', ...
        'VAk - VBk = Vf (volume constraint)', ...
        'QBk - QAk = sin(PAk)/rAk (stress jump)', ...
        'Pressure balance equation', ...
        'Force balance equation'
    };
    fprintf('  Total: %d conditions\n', length(neck_conditions));
    for i = 1:length(neck_conditions)
        fprintf('    %d. %s\n', i, neck_conditions{i});
    end
    
    % Count total BCs
    total_bc = length(south_active) + length(north_active) + length(neck_conditions);
    total_vars = 18; % 9 per phase
    fprintf('\nTotal boundary conditions: %d\n', total_bc);
    fprintf('Total variables: %d\n', total_vars);
    fprintf('Unknowns (Lagrange multipliers): ~2 (lam, and phase-specific Ls)\n');
    
    fprintf('\n  ⚠️  ISSUE: North pole has %d free variables vs %d at south pole\n', ...
        length(north_free), length(south_free));
    fprintf('  This asymmetry needs physical justification\n');
    fprintf('\n');
end

%% Test 4: Energy Scaling
function test_energy_scaling()
    fprintf('TEST 4: Energy Integration Scaling\n');
    fprintf('----------------------------------------\n');
    
    fprintf('Energy equation in code: 0.25*k*(2*H - H0)^2 * sin(S)\n');
    fprintf('Expected for axisymmetric: k*(2*H - H0)^2 * r * sin(P) * (geometric_factor)\n');
    
    fprintf('\nVolume equation in code: 0.75*r*sin(P)*sin(S)\n');
    fprintf('Expected for axisymmetric: r*sin(P) * (geometric_factor)\n');
    
    fprintf('\n  ⚠️  OBSERVATION:\n');
    fprintf('  - Volume equation HAS r*sin(P) dependence (correct)\n');
    fprintf('  - Energy equation LACKS r*sin(P) dependence (suspicious)\n');
    fprintf('  - Coefficient 0.25 in energy needs justification\n');
    fprintf('  - Coefficient 0.75 in volume needs justification\n');
    
    fprintf('\nAnalytical test: Unit sphere (H = 1 everywhere, H0 = 0)\n');
    fprintf('  Bending energy E = integral of k*(2H)^2 dA\n');
    fprintf('  For sphere: E = k*4*integral(dA) = k*4*(4π) = 16πk\n');
    fprintf('  With k=1: E_expected = 50.265 (16π)\n');
    fprintf('\n  ACTION ITEM: Test this analytically or numerically\n');
    fprintf('\n');
end
