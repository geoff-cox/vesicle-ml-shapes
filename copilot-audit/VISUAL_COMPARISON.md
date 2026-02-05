# Visual Comparison: Pole vs Bulk Equations

This document provides side-by-side comparison of the pole and bulk implementations to highlight mathematical inconsistencies.

---

## Issue 1: Q-Equation Inconsistency

### Bulk Form (BendV_Lag_EIGp_DE_impl.m, line 38)

```matlab
RHS = @(Q, H, P, r, z, L, s, V, B, S, k, H0) [ ...
    (-Q*cos(P)/r 
     - k*(2*H - H0)*(H*H0 + 2*(H - sin(P)/r)^2) 
     + 2*H*L 
     + lam
    )*sin(S)/r;
    ...
];
```

**Terms in bulk Q-equation (inside sin(S)/r factor):**
1. `-Q*cos(P)/r` - Geometric stress term
2. `-k*(2*H - H0)*(H*H0 + 2*(H - sin(P)/r)^2)` - Bending force
3. `+ 2*H*L` - Pressure-curvature coupling  ⬅️ NOTE: **2×H×L**
4. `+ lam` - Lagrange multiplier  ⬅️ NOTE: **lam** (coefficient 1)

---

### Pole Form (BendV_Lag_EIGp_DE_impl.m, line 26)

```matlab
RHS_pole = @(Q, H, P, r, z, L, s, V, B, S, k, H0, phase) [ ...
    H*L 
    + 0.5*lam 
    - k*H0*H^2 
    + 0.5*k*H*H0^2;
    ...
];
```

**Terms in pole Q-equation:**
1. `H*L` - Pressure-curvature coupling  ⬅️ **WRONG: Should be 2×H×L**
2. `+ 0.5*lam` - Lagrange multiplier  ⬅️ **WRONG: Should be lam**
3. `- k*H0*H^2` - Partial bending term?
4. `+ 0.5*k*H*H0^2` - Partial bending term?

---

### Side-by-Side Comparison

| Term | Bulk | Pole | Match? |
|------|------|------|--------|
| Pressure coupling | `2*H*L` | `H*L` | ❌ Factor of 2 error |
| Lagrange mult. | `lam` | `0.5*lam` | ❌ Factor of 2 error |
| Bending force | `k*(2*H - H0)*(...)` | `-k*H0*H^2 + 0.5*k*H*H0^2` | ❌ Different form |
| Geometric term | `-Q*cos(P)/r` | (absent) | ⚠️ Should vanish at pole |

---

### Mathematical Expectation

At the pole (r → 0, P → 0), the bulk equation should simplify to a regular form. 

**Expected pole limit:**
```matlab
% As r→0, P→0:
% - Q*cos(P)/r → finite (Q also→0)  
% - sin(P)/r → 1 (using L'Hôpital)
% - H - sin(P)/r → H - 1

% Pole limit should be approximately:
2*H*L + lam - k*(2*H - H0)*(H*H0 + 2*(H-1)^2)
```

This does NOT match the current pole implementation!

---

## Issue 2: Dimensional Analysis

### Units of Each Term

Assume standard units:
- H: [1/Length] (curvature)
- L: [Pressure] = [Force/Length²]
- k: [Energy] = [Force × Length]
- lam: [Pressure] = [Force/Length²]
- Q: [dimensionless]

### Bulk Q-Equation (before sin(S)/r factor)

```matlab
-Q*cos(P)/r                                      → [1/L]
-k*(2*H - H0)*(H*H0 + 2*(H - sin(P)/r)^2)       → [F·L]·[1/L]·[1/L²] = [F/L²]
+2*H*L                                           → [1/L]·[F/L²] = [F/L³]
+lam                                             → [F/L²]
```

**Problem:** Terms have mixed units [F/L²] and [F/L³] !

**Possible resolution:** 
- Maybe L actually has units [F/L³] not [F/L²]?
- Or there's a hidden r factor in the definition?

### Pole Q-Equation

```matlab
H*L              → [1/L]·[F/L²] = [F/L³]  (if L is pressure)
0.5*lam          → [F/L²]
-k*H0*H^2        → [F·L]·[1/L]·[1/L²] = [F/L²]
+0.5*k*H*H0^2    → [F·L]·[1/L]·[1/L²] = [F/L²]
```

**Same problem:** H*L has units [F/L³] while other terms have [F/L²]

**This is mathematically inconsistent!**

---

## Issue 3: Energy Equation Comparison

### Bulk Energy (line 46)

```matlab
0.25*k*(2*H - H0)^2 * sin(S);
```

### Pole Energy (line 34)

```matlab
0.25*k*(2*H - H0)^2 * sin(S);
```

**They're identical!** This is suspicious because:

1. **Expected for axisymmetric:**
   ```
   dE/ds = energy_density × dA/ds = κ(2H - H₀)² × 2πr sin(P)
   ```

2. **Compare with Volume equation (lines 33, 45):**
   ```matlab
   0.75*r*sin(P)*sin(S);  % ← HAS r and sin(P) factors!
   ```

3. **Question:** Why does Volume integration include `r*sin(P)` but Energy doesn't?

### Coefficient Mystery

**Energy:** `0.25*k*(2*H - H0)^2`  
**Volume:** `0.75*r*sin(P)`

Where do 0.25 and 0.75 come from?

**Standard axisymmetric volume element:**
```
dV = π*r² * dz = π*r² * sin(P) * ds
```

For a shell: dV ≈ 2π*r * (thickness) * sin(P) * ds

**Standard energy:**
```
dE = κ*(2H - H₀)² * dA = κ*(2H - H₀)² * 2π*r*sin(P) * ds
```

The factors 0.25 and 0.75 suggest some normalization or scaling, but it's not documented.

---

## Issue 4: Boundary Condition Asymmetry

### South Pole (α-phase)

```matlab
res_south = [
    QAs        % ✓ = 0
    % HAs     % ✗ FREE
    PAs        % ✓ = 0
    rAs        % ✓ = 0
    zAs        % ✓ = 0
    % LAs     % ✗ FREE
    sAs        % ✓ = 0
    VAs        % ✓ = 0
    EAs        % ✓ = 0
];
```

**Constrained:** 7 variables  
**Free:** 2 variables (H, L)

---

### North Pole (β-phase)

```matlab
res_north = [
    QBn        % ✓ = 0
    % HBn     % ✗ FREE
    PBn - pi   % ✓ = π
    % rBn     % ✗ FREE ← Why free?
    % zBn     % ✗ FREE ← Why free?
    % LBn     % ✗ FREE
    sBn        % ✓ = 0
    VBn        % ✓ = 0
    EBn        % ✓ = 0
];
```

**Constrained:** 5 variables  
**Free:** 4 variables (H, r, z, L)

---

### Why the Asymmetry?

**Hypotheses:**

1. **z-offset:** South pole is coordinate origin (z=0), but north pole z is determined by integration ✓ Makes sense

2. **r-coordinate:** 
   - At south pole: r=0 (on axis) ✓
   - At north pole: r=? 
     - If vesicle returns to axis: should be r=0
     - If vesicle has tubular neck: r>0 at north pole
   - Current implementation leaves it FREE

3. **This needs validation:** Check actual solved solutions to see if r(north) → 0 or not

**Action:** Run `inspect_solution_poles.m` to check!

---

## Issue 5: Neck Force Balance

### Current Implementation (line 61-64)

```matlab
kA*(2*HAk - H0(1))*(HAk - sin(PAk)/rAk + H0(1)/2) ...
- kB*(2*HBk - H0(2))*(HBk - sin(PAk)/rAk + H0(2)/2) ...  % ← Uses PAk!
- (cos(PAk)/rAk) + LBk - LAk
```

**Notice:** Both terms use `sin(PAk)/rAk` (α-phase geometry)

### Boundary Conditions at Neck

From lines 52-57:
```matlab
res_neck = [
    PBk - PAk          % ← Enforces PBk = PAk
    rBk - rAk          % ← Enforces rBk = rAk
    zBk - zAk          % ← Enforces zBk = zAk
    ...
];
```

### Analysis

Since the BCs enforce `PBk = PAk` and `rBk = rAk`, we have:
```
sin(PBk)/rBk = sin(PAk)/rAk
```

So using `PAk` for both phases is **mathematically equivalent** to using phase-specific values.

**However:**
- The formulation is **not manifestly symmetric**
- It creates an artificial coupling
- Makes independent verification harder
- If BCs are modified, this could break

**Verdict:** Mathematically correct but poor code design

---

## Recommended Fixes

### Fix 1: Correct Q-Equation in Pole (CRITICAL)

**Current (line 26):**
```matlab
H*L + 0.5*lam - k*H0*H^2 + 0.5*k*H*H0^2;
```

**Proposed fix (pending derivation verification):**
```matlab
2*H*L + lam - k*(2*H - H0)*(H*H0 + 2*H^2);
```

**Rationale:**
- Match coefficient factors (2×H×L, lam)
- Use consistent bending force form
- Ensure dimensional consistency

**IMPORTANT:** This fix is tentative and must be verified against the original mathematical derivation!

---

### Fix 2: Add Documentation Comments

Add detailed comments explaining:

```matlab
% Boundary conditions at poles:
% - H free: Required for regularity (H = limiting value of mean curvature)
% - L free: Lagrange multiplier determined by global volume constraint
% - At south pole (s=0): r=0, z=0 (coordinate origin on axis)
% - At north pole (s=π): r, z free (determined by integration)
%   Note: For on-axis vesicles, r naturally → 0. For necked/tubular shapes, r > 0.
```

---

### Fix 3: Validate Energy Integration

Test with analytical case:
```matlab
% Test: Perfect sphere with R=1, H=1/R=1, H0=0
% Expected energy: E = κ * ∫(2H)² dA = κ * 4 * 4πR² = 16πκ
% With κ=1: E_expected = 50.265
```

Run this test to verify the 0.25 coefficient and missing r-dependence.

---

## Summary Table of Issues

| # | Issue | File | Line | Severity | Status |
|---|-------|------|------|----------|--------|
| 1 | Wrong coefficients in Q-pole | DE_impl.m | 26 | CRITICAL | Needs fix |
| 2 | Dimensional inconsistency | DE_impl.m | 26 | CRITICAL | Needs fix |
| 3 | Asymmetric phase coupling | BC_impl.m | 62 | MODERATE | Need clarification |
| 4 | Underconstrained north pole | BC_impl.m | 36-38 | MODERATE | Need validation |
| 5 | Missing r in energy | DE_impl.m | 34,46 | MODERATE | Need derivation |
| 6 | Poor documentation | BC_impl.m | various | MODERATE | Add comments |
| 7 | Unused phase parameter | DE_impl.m | 29 | MINOR | Clean up |

---

**Next Steps:**
1. Consult original derivation to verify correct pole expansion
2. Run `diagnostics_bc_de_math.m` tests
3. Run `inspect_solution_poles.m` on existing solutions
4. Implement fixes
5. Re-validate all solutions
