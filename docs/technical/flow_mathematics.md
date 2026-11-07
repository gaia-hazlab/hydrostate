# Flow equation, finite-volume residual and linearization

## Defined physical and mathematical scope

`hydrostate.flow.RichardsModel` operates on a uniform Cartesian 3-D grid with
an explicit active-cell mask. Arrays have order `[z, y, x]`; coordinates ascend
and `z` is elevation in meters. Inactive interfaces and unspecified exterior
faces are impermeable. The present conductivity is scalar/isotropic. This is
not a rotated terrain-coordinate solver, general unstructured discretization,
full anisotropic tensor implementation or a MODFLOW backend.

For constant-density water with atmospheric pore air, the implemented mixed
balance is

```text
H = z + psi                         total hydraulic head [m]
q = -K(psi) grad(H)                  Darcy volume flux [m/s]
partial W(psi)/partial t + div(q) = f
W(psi) = theta(psi) + Ss max(psi, 0)
```

`f` has units 1/s and is positive into the cell. The public kernel instead
accepts integrated `sources_m3_s = f V`. Root extraction is a negative source
if included by a caller. Pressure head is distinct from total head, and
`u_w - u_a = rho_w g psi`. The saturated `Ss` term is an explicit piecewise
linear elastic-storage approximation, not an additional specific-yield term.
It is counted in stored water and never substituted for retention drainage.

The implementation also exposes an explicitly mathematical linear mode:
constant conductivity and `W = C psi`, with `C > 0`. This storage is **not**
physical volumetric soil water content. Linear verification deliberately uses
the same topology, boundary preparation, face fluxes, residual, assembly and
advance methods as nonlinear valley runs.

## Constitutive closure

For `psi < 0`, the van Genuchten–Mualem laws use

```text
m = 1 - 1/n; a = -alpha psi; Se = (1+a^n)^(-m)
theta = theta_r + (theta_s-theta_r) Se
K = Ks sqrt(Se) [1 - (1-Se^(1/m))^m]^2
```

For `psi >= 0`, `Se=1`, `theta=theta_s`, and `K=Ks`.
Require `n>1`, `alpha>0`, `Ks>0`, `Ss>=0` and admissible water contents.
Effective saturation `Se` is not liquid saturation `theta/porosity`.
[van Genuchten, 1980](https://doi.org/10.2136/sssaj1980.03615995004400050002x),
[Mualem, 1976](https://doi.org/10.1029/WR012i003p00513).

On the smooth unsaturated branch,

```text
Se' = m n alpha a^(n-1) (1+a^n)^(-m-1)
theta' = (theta_s-theta_r) Se'
B = 1-(1-Se^(1/m))^m
dB/dSe = (1-Se^(1/m))^(m-1) Se^(1/m-1)
K' = Ks [0.5 Se^(-0.5) B^2 + 2 sqrt(Se) B dB/dSe] Se'
```

The independent verification oracle uses this chain rule, and compares it
with numerical conductivity derivatives of the actual kernel. Storage and
capacity are checked on both sides of saturation, away from the junction.
No differentiability claim is made at `psi=0`: the elastic capacity is
piecewise, and the chosen conductivity derivative can become singular as the
unsaturated branch approaches saturation. There is no hysteresis, ice,
preferential-flow domain, vapor or variable-density closure in this kernel.

## Discrete flux and backward Euler residual

Every internal face has one stored oriented volume flux. For adjacent active
cells `i` and `j`, center distance `d`, and face area `A`,

```text
Kij = 2 Ki Kj / (Ki+Kj)
Tij = (A/d) Kij                      [m2/s]
Qij = Tij (Hi-Hj)                    [m3/s, positive i -> j]
```

The same `Qij` enters cell `i` with positive and cell `j` with negative sign.
For a prescribed total head `Hb` at an exterior face,
`Qib = Ki A/(d/2) (Hi-Hb)`. This is a half-cell connection using current cell
conductivity. A prescribed flux boundary supplies an **outward integrated
volume rate**, not a depth rate; callers must multiply a depth flux by area.
Distinct faces of a corner cell may have distinct boundary conditions, while
duplicate conditions on the same physical face are rejected.

The solved residual is

```text
Ri = Vi [W(psi_i_new)-W(psi_i_old)]/dt
     + sum_internal Qij + sum_boundary Qib - si = 0.
```

The actual nonlinear storage difference is evaluated at every residual call.
Replacing it by `C_old (psi_new-psi_old)` would solve a different equation and
can destroy water conservation. The mixed formulation is motivated by the
storage/conservation analysis of
[Celia et al., 1990](https://doi.org/10.1029/WR026i007p01483).

Summing the **actual** residual over active cells gives

```text
sum Ri = [sum Vi Delta W]/dt + sum exterior Qib - sum si.
```

Internal exchanges cancel algebraically. This identity explains closure,
but does not prove spatial accuracy: an inaccurate face approximation can
still cancel exactly. The runner therefore also compares analytical fields,
individual analytical Darcy fluxes and refinement orders.

## Modified Picard matrix versus true Newton Jacobian

The solver matrix freezes conductances at the current trial state while
retaining current storage capacity:

```text
Aii = Vi W'(psi_i)/dt + sum_j Tij + sum_head_boundaries Tib
Aij = -Tij
```

A small capacity floor is used in this *iteration matrix*; it does not replace
storage in the residual or add water to the accepted state. State-dependent
conductivity is recalculated on subsequent iterations. Damping limits trial
pressure changes, and accepted states must meet both local volume-residual and
pending Picard-correction tolerances. The pending correction is computed even
on the first trial: a small absolute volume residual does not imply a small
head error when capacity/conductance is tiny. This is an approximate
linearized correction diagnostic, not a bound from the true Newton Jacobian.
There is no accepted-state clipping.

This is not the full nonlinear Newton Jacobian. For an internal face, the true
derivatives also contain

```text
dQij/dpsi_i = Tij + (dTij/dpsi_i)(Hi-Hj)
dQij/dpsi_j = -Tij + (dTij/dpsi_j)(Hi-Hj),
```

with analogous conductivity derivatives at head boundaries. Omitting these
terms is a stated Picard linearization, not a coding mistake if the exact
nonlinear residual converges. In linear mode those derivative terms vanish:
the actual solver matrix must equal the independently differentiated residual.
The verification suite checks that equality, symmetry, positive eigenvalues,
nonpositive neighbor entries and storage row sums for a small no-flow system.

With positive linear storage, the matrix is
`diag(CV/dt) + B^T diag(T) B + boundary_diagonal`, hence positive definite.
That argument does not prove nonlinear Newton symmetry or universal Picard
convergence. A purely saturated incompressible sealed problem can require a
head reference and lies outside that positive-storage argument.

## Geometry and floating-point limits

The two-point flux approximation suits the stated orthogonal scalar-K mesh.
Mask steps are actual no-flow/exposed Cartesian faces, not a smooth terrain
surface. Hydrostatic tests on masks check topology and gravity cancellation;
they do not quantify stair-step geometry error. Harmonic means across sharp
fronts/layers and cell-K exterior connections require solution refinement.
Denominator guards and floating-point underflow at extremely small
conductivities also delimit dry-state behavior; the smooth verification range
does not establish reliable behavior at arbitrarily large suction.
