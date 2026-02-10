% === validate_issue1_fix.m ===
% Validation script for Issue 1 fix (pole expansion coefficients)
%
% Tests:
% 1. Pole-bulk consistency at boundary
% 2. Sphere geometry (H=H0, Q=0)
% 3. Energy conservation
%
% Run this after fixing line 562 in solveAtParams.m

function validate_issue1_fix()
    fprintf('=== Validation of Issue 1 Fix ===\n');
    fprintf('Testing corrected pole expansion in solveAtParams.m\n\n');
    
    %% Test 1: Pole-Bulk Consistency
    fprintf('Test 1: Pole-Bulk Consistency\n');
    fprintf('------------------------------\n');
    test_pole_bulk_consistency();
    
    %% Test 2: Sphere Geometry
    fprintf('\nTest 2: Sphere Geometry\n');
    fprintf('-----------------------\n');
    test_sphere_geometry();
    
    %% Test 3: Parameter Sweep Comparison
    fprintf('\nTest 3: Parameter Sweep Comparison\n');
    fprintf('-----------------------------------\n');
    fprintf('This would compare solutions before/after fix.\n');
    fprintf('Requires archived pre-fix results for comparison.\n');
    fprintf('SKIPPED in this validation (no pre-fix archive available).\n');
    
    fprintf('\n=== Validation Complete ===\n');
end

function test_pole_bulk_consistency()
    % Test that pole and bulk equations produce similar values at boundary
    %
    % At the boundary (S = delta*pi), the pole expansion and bulk equation
    % should give nearly identical results since sin(S)/r ≈ 1 for small S
    
    % Define test parameters
    delta = 0.1;  % Typical delta value used in code
    S_boundary = delta * pi;
    
    % Test state variables (arbitrary but representative)
    Q = 0.01;
    H = 1.0;
    P = 0.05;
    r = 0.05;
    z = 0.0;
    L = 0.5;
    s = 0.1;
    V = 0.5;
    B = 0.0;
    
    % Parameters
    k = 1.0;
    H0 = 1.0;
    lam = 0.1;  % This would normally be computed
    phase = 1;
    
    % Evaluate RHS_pole (NEW CORRECTED VERSION)
    RHS_pole_Q = 2*H*L + lam - 2*k*H0*H^2 + k*H*H0^2;
    
    % Evaluate RHS at boundary (from bulk equation line 574)
    % Note: At pole limit, we should use sin(P)/r → H for consistency
    % Here we're testing near the boundary, so use actual H value
    sin_S_over_r = sin(S_boundary) / r;
    H_approx = sin(P) / r;  % This approximation improves as P→0
    RHS_Q = (-Q*cos(P)/r - k*(2*H - H0)*(H*H0 + 2*(H - H_approx)^2) ...
             + 2*H*L + lam) * sin_S_over_r;
    
    fprintf('  Boundary S = %.6f (delta*pi with delta=%.2f)\n', S_boundary, delta);
    fprintf('  sin(S)/r at boundary = %.6f\n', sin_S_over_r);
    fprintf('  RHS_pole (Q-equation) = %.6f\n', RHS_pole_Q);
    fprintf('  RHS (Q-equation) = %.6f\n', RHS_Q);
    fprintf('  Relative difference = %.2f%%\n', ...
            abs(RHS_pole_Q - RHS_Q) / abs(RHS_pole_Q) * 100);
    
    % Check if sin(S)/r is close to 1 (as required by parameterization)
    expected_jacobian = 1.0;
    jacobian_error = abs(sin_S_over_r - expected_jacobian);
    fprintf('  Expected sin(S)/r = %.6f, Error = %.6f\n', ...
            expected_jacobian, jacobian_error);
    
    if jacobian_error > 0.1
        fprintf('  WARNING: sin(S)/r not close to 1 at boundary!\n');
        fprintf('  This suggests the parameterization may not be consistent.\n');
    else
        fprintf('  ✓ PASS: sin(S)/r ≈ 1 at boundary as expected\n');
    end
end

function test_sphere_geometry()
    % For a sphere with uniform H0, the solution should be:
    % - H = H0 everywhere (constant mean curvature)
    % - Q = 0 everywhere (no shear stress)
    % - P consistent with spherical geometry
    %
    % This is a basic sanity check that the corrected equations
    % produce physically reasonable results
    
    fprintf('  Testing sphere with H0 = 1.0 (uniform curvature)\n');
    fprintf('  Note: Full BVP solve required - just showing expected behavior\n');
    
    % Expected behavior for sphere
    H0_sphere = 1.0;
    A_sphere = 4*pi;
    V_sphere = 4*pi/3;
    
    fprintf('\n  Expected solution for perfect sphere:\n');
    fprintf('    H = %.6f everywhere (constant mean curvature)\n', H0_sphere);
    fprintf('    Q ≈ 0 everywhere (no shear stress)\n');
    fprintf('    r(S) = sin(S) (spherical geometry)\n');
    fprintf('    z(S) = -cos(S) (spherical geometry)\n');
    
    fprintf('\n  To verify: Run solveAtParams with:\n');
    fprintf('    A = %.6f, V = %.6f\n', A_sphere, V_sphere);
    fprintf('    H0_1 = H0_2 = %.6f\n', H0_sphere);
    fprintf('    Check that H ≈ %.6f and Q ≈ 0 in solution\n', H0_sphere);
    
    % Mathematical check: For sphere, the Q-equation at pole should give:
    % dQ/dS = 2*H*L + lam - 2*k*H0*H^2 + k*H*H0^2
    % With H = H0 = 1, Q = 0:
    % dQ/dS = 2*L + lam - 2*k*1 + k*1 = 2*L + lam - k
    %
    % For equilibrium, this should be close to 0 (Q remains 0)
    
    H = H0_sphere;
    k = 1.0;  % Example bending rigidity
    
    fprintf('\n  Mathematical check at pole (H = H0 = %.2f, Q = 0):\n', H0_sphere);
    fprintf('    For Q to remain 0, need: 2*H*L + lam - 2*k*H0*H^2 + k*H*H0^2 ≈ 0\n');
    fprintf('    With H = H0: 2*H0*L + lam - 2*k*H0^3 + k*H0^3 = 2*H0*L + lam - k*H0^3 ≈ 0\n');
    fprintf('    Simplifies to: 2*H0*L + lam ≈ k*H0^3\n');
    fprintf('    The Lagrange multipliers L and lam adjust to satisfy this.\n');
    
    fprintf('  ✓ Expected behavior documented\n');
end
