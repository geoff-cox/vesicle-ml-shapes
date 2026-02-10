# Issue 1: Corrected Fix and Validation Plan

**Date:** February 10, 2026  
**Status:** CRITICAL BUG CONFIRMED - FIX REQUIRED

---

## Summary

Issue 1 from `copilot-audit/MATHEMATICAL_ANALYSIS.md` is **VALID**. The pole expansion in `src/utils/solveAtParams.m` has incorrect coefficients in the Q-equation.

**All four coefficients are wrong by a factor of 2.**

---

## The Error

### Location
- **File**: `src/utils/solveAtParams.m`
- **Line**: 562 (first line of RHS_pole anonymous function)

### Current Code (INCORRECT)
```matlab
RHS_pole = @(Q, H, P, r, z, L, s, V, B, S, k, H0, phase) [ ...
    H*L + 0.5*lam - k*H0*H^2 + 0.5*k*H*H0^2;  % ← LINE 562: ALL WRONG
    0;
    H;
    phase;
    0;
    0;
    1;
    0.75*r*sin(P)*sin(S);
    0.25*k*(2*H - H0)^2 * sin(S);
];
```

### Corrected Code (CORRECT)
```matlab
RHS_pole = @(Q, H, P, r, z, L, s, V, B, S, k, H0, phase) [ ...
    2*H*L + lam - 2*k*H0*H^2 + k*H*H0^2;  % ← CORRECTED
    0;
    H;
    phase;
    0;
    0;
    1;
    0.75*r*sin(P)*sin(S);
    0.25*k*(2*H - H0)^2 * sin(S);
];
```

### Specific Changes
| Coefficient | Current (Wrong) | Should Be (Correct) | Change |
|-------------|-----------------|---------------------|--------|
| H*L term | `H*L` | `2*H*L` | Multiply by 2 |
| λ term | `0.5*lam` | `lam` | Multiply by 2 |
| -k*H₀*H² term | `-k*H0*H^2` | `-2*k*H0*H^2` | Multiply by 2 |
| k*H*H₀² term | `0.5*k*H*H0^2` | `k*H*H0^2` | Multiply by 2 |

---

## Mathematical Justification

### The Parameterization

The code uses two different parameterizations:

**Pole Expansion (near S=0 or S=π):**
```matlab
% Line 568 of solveAtParams.m
ds/dS = 1
```

**Bulk Equation (S away from poles):**
```matlab
% Line 580 of solveAtParams.m  
ds/dS = sin(S)/r
```

### The Critical Insight

For these to be consistent at the pole boundary, we must have:
```
sin(S)/r = 1  at the pole
```

This means **there is NO factor of 1/2 Jacobian transformation** at the pole. The pole expansion should exactly match the bulk equation evaluated at the pole limit.

### Derivation of Correct Formula

Starting from bulk Q-equation (line 574):
```matlab
dQ/dS = (-Q*cos(P)/r - k*(2*H - H0)*(H*H0 + 2*(H - sin(P)/r)^2) 
         + 2*H*L + lam)*sin(S)/r
```

Taking the pole limit (S→0):
1. Q → 0 (no shear at pole)
2. P → 0 (tangent angle)
3. sin(P)/r → H (mean curvature)
4. **sin(S)/r → 1** (from parameterization consistency)

Substituting:
```
dQ/dS = (0 - k*(2*H - H0)*(H*H0 + 2*(H - H)^2) + 2*H*L + lam) × 1
      = -k*(2*H - H0)*H*H0 + 2*H*L + lam
      = -2*k*H^2*H0 + k*H*H0^2 + 2*H*L + lam
```

This is exactly what the corrected code has.

---

## Why the Error Occurred

Someone implementing the code likely:
1. Incorrectly assumed `sin(S)/r = 1/2` at the pole
2. Multiplied the bulk equation by 1/2 to get the pole expansion
3. Failed to verify against the parameterization

The evidence: all coefficients are exactly half of what they should be.

---

## Impact Analysis

### Severity
**CRITICAL** - This is a fundamental error in the physics equations.

### What's Affected
- ✗ All solutions in `SimResults/` catalog
- ✗ Stress distributions (Q values throughout domain)
- ✗ Mean curvature fields (H values)
- ✗ Energy calculations
- ✗ Lagrange multipliers (L, λ)

### Why Results May Still Look Reasonable
1. **Relative scaling preserved**: All terms wrong by same factor
2. **Lagrange multipliers adapt**: L and λ are computed to satisfy constraints
3. **BVP solver converges**: Finds solutions to the (incorrect) equations
4. **Visual appearance**: Shapes may look qualitatively correct

### Why Results Are Still Wrong
- Absolute stress values incorrect
- Energy values incorrect  
- Cannot compare to theoretical predictions
- Physical interpretation compromised

---

## Fix Implementation

### Required Code Changes

**File**: `src/utils/solveAtParams.m`

**Change at line 562:**
```matlab
% OLD (line 562)
    H*L + 0.5*lam - k*H0*H^2 + 0.5*k*H*H0^2;

% NEW (line 562)  
    2*H*L + lam - 2*k*H0*H^2 + k*H*H0^2;
```

### Implementation Notes
1. This is a single-line change
2. No other files need modification
3. The fix is mathematically straightforward
4. All existing results should be invalidated

---

## Validation Plan

### 1. Unit Test: Sphere Geometry
A sphere with uniform H₀ should give:
- H = H₀ everywhere
- Q = 0 everywhere (no shear)
- P = arcsin(H*r) (consistent with geometry)

**Test:**
```matlab
% Solve for sphere: H0_1 = H0_2 = H0_sphere
% Verify H ≈ H0_sphere everywhere
% Verify Q ≈ 0 everywhere
```

### 2. Consistency Check: Pole vs Bulk
At the pole boundary (S = delta*π), the pole and bulk equations should produce nearly identical values.

**Test:**
```matlab
% Evaluate RHS_pole at S = delta*π
% Evaluate RHS at S = delta*π
% Verify |RHS_pole - RHS| < tolerance
```

### 3. Energy Conservation
Total energy should be constant along solution paths.

**Test:**
```matlab
% Integrate energy along S
% Verify E_total matches boundary values
```

### 4. Regression Test: Compare Solutions
Run a small parameter sweep before and after the fix.

**Test:**
```matlab
% Generate solutions with old code
% Apply fix
% Generate solutions with new code
% Compare shapes, energies, stress distributions
% Document differences
```

### 5. Published Results Check
If this code has been used in publications, compare to any published values.

---

## Post-Fix Actions Required

### Immediate
1. ✅ Apply the fix to `solveAtParams.m` line 562
2. ✅ Run validation tests (sphere, consistency, energy)
3. ✅ Document fix in commit message and code comments

### Short-term
1. ⚠️ Invalidate all results in `SimResults/` 
2. ⚠️ Clear catalog.mat and cache.mat (or mark as pre-fix)
3. ⚠️ Re-run representative parameter sweeps
4. ⚠️ Compare old vs new results quantitatively

### Long-term
1. ⚠️ Review any publications using old results
2. ⚠️ Consider publishing erratum if needed
3. ⚠️ Add regression tests to prevent future errors
4. ⚠️ Document the parameterization more clearly

---

## Additional Fixes Needed

While fixing line 562, also check if similar errors exist elsewhere:

### Other Files to Check
1. `src/utils/solveAtParams_v2.m` (if it exists)
2. `src_v2/` directory (any alternative implementations)
3. Any other files that implement RHS_pole

### Pattern to Look For
- Factors of 0.5 or 2 that might indicate Jacobian confusion
- Comments mentioning "sin(S)/r = 1/2" (incorrect)
- Any hardcoded Jacobian factors

---

## Documentation Updates Needed

### In Code Comments
Add to `solveAtParams.m` near line 560:
```matlab
% CRITICAL: Pole expansion parameterization
% -----------------------------------------
% The parameterization requires sin(S)/r = 1 at poles (S=0 or S=π)
% This is enforced by setting ds/dS = 1 in RHS_pole (line 568)
% vs ds/dS = sin(S)/r in bulk RHS (line 580)
%
% Therefore: NO factor of 1/2 Jacobian transformation at pole
% Pole coefficients match bulk equation evaluated at pole limit
%
% Previous versions incorrectly assumed sin(S)/r = 1/2 at pole
% This was fixed on [DATE] - see Issue 1 assessment docs
```

### In Project Documentation
Update `PROJECT_AUDIT.md` to mention:
- The parameterization scheme (sin(S)/r = 1 at poles)
- The fix applied to line 562
- That pre-fix results are invalid

---

## Testing Script Template

```matlab
%% Validation Script for Issue 1 Fix
% Tests the corrected pole expansion

%% Test 1: Sphere geometry
fprintf('Test 1: Sphere geometry\n');
params = struct('A', 4*pi, 'V', 4*pi/3, 'KA', 1, 'KB', 1, 'KG', 0, ...
                'H0_1', 1, 'H0_2', 1);
% ... solve and verify H ≈ 1, Q ≈ 0

%% Test 2: Pole-bulk consistency  
fprintf('Test 2: Pole-bulk consistency\n');
% ... evaluate both RHS at boundary

%% Test 3: Energy conservation
fprintf('Test 3: Energy conservation\n');
% ... integrate and verify

%% Test 4: Regression comparison
fprintf('Test 4: Compare to old results\n');
% ... load old solution, compare shapes
```

---

## Summary Checklist

- [x] Error identified and understood
- [x] Mathematical derivation of correct formula
- [x] Impact analysis completed
- [ ] Fix applied to code
- [ ] Validation tests run
- [ ] Results compared (old vs new)
- [ ] Documentation updated
- [ ] Existing results invalidated/marked
- [ ] Publications reviewed (if applicable)

---

## References

- **Error location**: `src/utils/solveAtParams.m`, line 562
- **Original audit**: `copilot-audit/MATHEMATICAL_ANALYSIS.md`, lines 55-105
- **Corrected assessment**: `ASSESSMENT_ISSUE_1.md`
- **Parameterization proof**: Lines 568, 580 of solveAtParams.m

---

**Next Step**: Apply the fix and run validation tests.
