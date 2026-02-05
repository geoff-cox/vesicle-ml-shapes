# Annotated Source Code: Issues Identified

This document shows the exact problematic lines with detailed annotations.

---

## File: BendV_Lag_EIGp_DE_impl.m

### Lines 25-35: RHS_pole Function - CRITICAL ISSUES

```matlab
24:     % RHS_pole handles singular pole expansion; RHS is bulk form.
25:     RHS_pole = @(Q, H, P, r, z, L, s, V, B, S, k, H0, phase) [ ...
26:         H*L + 0.5*lam - k*H0*H^2 + 0.5*k*H*H0^2;    % ❌ CRITICAL ERROR
27:         0;
28:         H;
29:         phase;                                       % ⚠️ UNUSED PARAMETER
30:         0;
31:         0;
32:         1;
33:         0.75*r*sin(P)*sin(S);
34:         0.25*k*(2*H - H0)^2 * sin(S);               % ⚠️ MISSING r*sin(P)?
35:     ];
```

#### Line 26 - CRITICAL: Multiple Errors in Q-Equation

**Current code:**
```matlab
H*L + 0.5*lam - k*H0*H^2 + 0.5*k*H*H0^2;
```

**Problems:**

1. **❌ Wrong coefficient on H*L:**
   - Current: `H*L`
   - Expected: `2*H*L` (factor of 2 missing)
   - Compare with bulk (line 38): uses `2*H*L`

2. **❌ Wrong coefficient on lam:**
   - Current: `0.5*lam`
   - Expected: `lam` (factor of 2 error)
   - Compare with bulk (line 38): uses `lam`

3. **❌ Wrong bending force form:**
   - Current: `-k*H0*H^2 + 0.5*k*H*H0^2`
   - This does NOT match d/dH[κ(2H - H₀)²]
   - Derivative: d/dH[κ(2H - H₀)²] = 2κ(2H - H₀)·2 = 4κ(2H - H₀)
   - Expanding: 4κ(2H - H₀) = 8κH - 4κH₀
   - Current form: -κH₀H² + 0.5κHH₀² cannot simplify to this!

4. **❌ DIMENSIONAL INCONSISTENCY:**
   ```
   H*L:            [1/L]·[F/L²] = [F/L³]
   0.5*lam:        [F/L²]
   k*H0*H^2:       [F·L]·[1/L]·[1/L²] = [F/L²]
   0.5*k*H*H0^2:   [F·L]·[1/L]·[1/L²] = [F/L²]
   ```
   Terms have incompatible units!

**Proposed fix (pending verification):**
```matlab
2*H*L + lam - k*(2*H - H0)*(H*H0 + 2*H^2);
```
This assumes pole limit with sin(P)/r → 1 at the pole.

---

#### Line 29 - MINOR: Unused Parameter

```matlab
phase;    % ⚠️ Parameter passed but never used
```

**Issue:** The `phase` variable (+1 for α, -1 for β) is evaluated but doesn't affect any calculation.

**Action:** Either remove it or implement its intended purpose (e.g., sign of curvature normal).

---

#### Line 34 - MODERATE: Missing Geometric Factors?

```matlab
0.25*k*(2*H - H0)^2 * sin(S);    % ⚠️ Where is r*sin(P)?
```

**Issue:** Energy integration for axisymmetric surface should include r*sin(P):
```
dE/ds = κ(2H - H₀)² × dA/ds = κ(2H - H₀)² × 2πr sin(P)
```

**Compare with Volume (line 33):**
```matlab
0.75*r*sin(P)*sin(S);    % ✓ Includes r*sin(P)
```

**Questions:**
1. Why does volume include r*sin(P) but energy doesn't?
2. What is the 0.25 coefficient? (Expected: 2π for surface area)
3. Is there a coordinate transformation that absorbs r?

**Action:** Verify against analytical test case (sphere).

---

### Lines 37-47: RHS (Bulk) Function

```matlab
37:     RHS = @(Q, H, P, r, z, L, s, V, B, S, k, H0) [ ...
38:         (-Q*cos(P)/r - k*(2*H - H0)*(H*H0 + 2*(H - sin(P)/r)^2) + 2*H*L + lam)*sin(S)/r;
39:         Q/(2*k)*sin(S)/r;
40:         (2*H - sin(P)/r)*sin(S)/r;
41:         cos(P)*sin(S)/r;
42:         sin(P)*sin(S)/r;
43:         0;
44:         sin(S)/r;
45:         0.75*r*sin(P)*sin(S);
46:         0.25*k*(2*H - H0)^2 * sin(S);
47:     ];
```

#### Line 38 - Reference for Comparison

This is the correct bulk form. Note:
- Uses `2*H*L` (not H*L)
- Uses `lam` (not 0.5*lam)
- Includes full bending force term

#### Line 46 - Same Energy Issue

```matlab
0.25*k*(2*H - H0)^2 * sin(S);    % ⚠️ Same as pole (line 34)
```

Both pole and bulk use identical energy formula. This is suspicious because:
- At pole, geometric factors should simplify differently
- Missing r*sin(P) factors

---

## File: BendV_Lag_EIGp_BC_impl.m

### Lines 14-27: South Pole Boundary Conditions

```matlab
14:     % α-phase (s)outh pole conditions
15:     south_pole = num2cell(y_poles(1:9));
16:     [QAs, HAs, PAs, rAs, zAs, LAs, sAs, VAs, EAs] = deal(south_pole{:});
17:     res_south = [
18:         QAs        % ✓ Q = 0 at pole
19:         % HAs     % ⚠️ FREE (no constraint)
20:         PAs        % ✓ P = 0 at south pole
21:         rAs        % ✓ r = 0 (on axis)
22:         zAs        % ✓ z = 0 (origin)
23:         % LAs     % ⚠️ FREE (Lagrange multiplier)
24:         sAs        % ✓ s = 0
25:         VAs        % ✓ V = 0 (volume starts at 0)
26:         EAs        % ✓ E = 0 (energy starts at 0)
27:     ];
```

**Status:** ✓ LOOKS CORRECT

**Constrained:** 7 variables (Q, P, r, z, s, V, E)  
**Free:** 2 variables (H, L)

**Rationale:**
- H free: Required for regularity (limiting curvature value)
- L free: Lagrange multiplier determined by global constraint

---

### Lines 29-42: North Pole Boundary Conditions - MODERATE CONCERN

```matlab
29:     % β-phase (n)orth pole conditions
30:     north_pole = num2cell(y_poles(10:18));
31:     [QBn, HBn, PBn, rBn, zBn, LBn, sBn, VBn, EBn] = deal(north_pole{:});
32:     res_north = [
33:         QBn        % ✓ Q = 0 at pole
34:         % HBn     % ⚠️ FREE (no constraint)
35:         PBn - pi   % ✓ P = π at north pole
36:         % rBn     % ⚠️ FREE - Why?
37:         % zBn     % ⚠️ FREE - Why?
38:         % LBn     % ⚠️ FREE (Lagrange multiplier)
39:         sBn        % ✓ s = 0
40:         VBn        % ✓ V = 0
41:         EBn        % ✓ E = 0
42:     ];
```

**Status:** ⚠️ ASYMMETRIC vs South Pole

**Constrained:** 5 variables (Q, P, s, V, E)  
**Free:** 4 variables (H, r, z, L)

#### Line 36-37 - Why are r and z FREE?

**Question:** Should north pole be on axis (r=0)?

**Possible scenarios:**
1. **On-axis vesicle:** North pole returns to axis → r should be 0
2. **Tubular neck:** North pole is at end of tube → r > 0

**Current implementation:** Leaves r free, letting solver decide.

**Action needed:**
1. Check solved solutions: does r(north) → 0 naturally?
2. Document design choice: why r, z are free
3. Consider adding constraint if physically required

**Tool:** Run `inspect_solution_poles.m` to check actual solutions.

---

### Lines 52-64: Neck Junction Conditions

```matlab
52:     res_neck = [
53:         PBk - PAk                          % ✓ Angle continuity
54:         rBk - rAk                          % ✓ Radial continuity
55:         zBk - zAk                          % ✓ Axial continuity
56:         (VAk - VBk) - Vf;                  % ✓ Volume constraint
57:         (QBk - QAk) - sin(PAk)/rAk;        % ✓ Stress jump
58:         kB*(2*HBk - H0(2)) ...             % Pressure balance
59:             - kA*(2*HAk - H0(1)) ...
60:             + kG*(sin(PAk)/rAk);
61:         kA*(2*HAk - H0(1))*(HAk - sin(PAk)/rAk + H0(1)/2) ...  % Force balance
62:         - kB*(2*HBk - H0(2))*(HBk - sin(PAk)/rAk + H0(2)/2) ...  % ⚠️ Uses PAk!
63:         - (cos(PAk)/rAk) + LBk - LAk
64:     ];
```

#### Line 62 - MODERATE: Asymmetric Phase Coupling

```matlab
- kB*(2*HBk - H0(2))*(HBk - sin(PAk)/rAk + H0(2)/2) ...  % ⚠️ Uses PAk, rAk
```

**Issue:** β-phase term uses α-phase geometry (`sin(PAk)/rAk`).

**Analysis:**
- Lines 53-54 enforce `PBk = PAk` and `rBk = rAk`
- Therefore: `sin(PBk)/rBk = sin(PAk)/rAk`
- Mathematically equivalent!

**But:**
- Not manifestly symmetric under phase exchange
- Harder to verify independently
- If BCs change, could break

**Better form (equivalent):**
```matlab
- kB*(2*HBk - H0(2))*(HBk - sin(PBk)/rBk + H0(2)/2) ...
```

**Status:** Mathematically correct but poor code style.

---

## Summary of Annotated Issues

### CRITICAL (Must Fix)

**BendV_Lag_EIGp_DE_impl.m, Line 26:**
- Wrong coefficients (H*L → 2*H*L, 0.5*lam → lam)
- Wrong bending force form
- Dimensional inconsistency

### MODERATE (Needs Verification)

**BendV_Lag_EIGp_DE_impl.m, Lines 34, 46:**
- Missing r*sin(P) in energy integration
- Unexplained coefficient 0.25

**BendV_Lag_EIGp_BC_impl.m, Lines 36-37:**
- North pole r, z left free (asymmetric with south pole)
- Needs physical justification or constraint

**BendV_Lag_EIGp_BC_impl.m, Line 62:**
- Uses α-phase geometry for both phases
- Mathematically OK but confusing

### MINOR (Clean Up)

**BendV_Lag_EIGp_DE_impl.m, Line 29:**
- Unused `phase` parameter

---

## Testing Checklist

- [ ] Run `diagnostics_bc_de_math.m` to verify dimensional analysis
- [ ] Run `inspect_solution_poles.m` to check if r(north) → 0
- [ ] Test energy conservation with analytical sphere case
- [ ] Verify pole-bulk limit numerically
- [ ] Consult original derivation papers
- [ ] Implement fixes for critical issues
- [ ] Re-run all simulations to validate

---

**End of Annotated Source Code Analysis**
