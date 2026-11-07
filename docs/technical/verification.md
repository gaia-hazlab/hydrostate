# Analytical verification experiments and measured scope

Run `python -m hydrostate.verification --output runs/verification/report.json`
with the project environment. The committed
[evidence](verification_summary.json) identifies the tested `flow.py` by SHA-256.
The runner uses that actual kernel; it does not solve an independent simplified
stencil and transfer the result to the valley model by assumption.

## Method and norms

A manufactured solution prescribes a smooth field, differentiates the
**continuum** equation analytically to obtain its source, and applies the
prescribed field at initial/boundary locations. Numerical error is then measured
against that field. Generating the source by applying the numerical residual
to the desired answer would make a tautological test and is not done here.
Order verification is established practice; it does not validate the physical
model against observations.
[Salari and Knupp, Sandia SAND2000-1444](https://doi.org/10.2172/759450),
[NASA verification report](https://ntrs.nasa.gov/citations/20080010138).

Equal-volume Cartesian cells permit volume-normalized
`L2 = sqrt(sum(error_i^2 Vi)/sum Vi)` and `Linf = max |error_i|`.
Inactive cells are excluded. `h=max(dx,dy,dz)` and observed pairwise order is
`p=log(E_coarse/E_fine)/log(h_coarse/h_fine)`; time orders use the time-step
ratio. The spatial grids use 5, 9 and 17 cells per axis and unequal physical
lengths `(1,1.2,0.8) m`. Odd grids preserve a cell at the sinusoid's exact
central maximum; this removes changing sampling of the maximum from the
linear maximum-norm order. Ratios need not be exactly two.

Required gates are specified in the runner and evaluated separately from
nonlinear diagnostics. Required spatial orders lie between 1.8 and 2.3 and
time orders between 0.85 and 1.15. Local solver tolerances are set well below
the discretization errors. The report records errors and orders, rather than
only a Boolean result; failure exits nonzero.

## Required analytical checks

1. **Constitutive and storage derivatives.** Independent analytical VG
   retention/capacity/conductivity formulas and the conductivity chain rule are
   compared against the kernel, including a centered difference of its actual
   conductivity. The derivative of `storage()` is checked against `capacity()`
   at smooth negative- and positive-pressure points. The branch junction is
   excluded explicitly.
2. **Affine 3-D Darcy flow.** `H=4+0.3x-0.2y+0.4z` has known constant vector
   flux and zero divergence for constant K. Check exact state preservation,
   every shared internal flux and every exterior head-boundary flux on an
   unequal-spacing grid. This exposes axis, gravity, sign, area and half-cell
   distance mistakes that a global water balance may miss.
3. **Prescribed source/flux.** A one-cell linear-storage problem has
   `Delta psi=(source-outward_flux) dt/(CV)`. This checks integrated source and
   outward flux conventions independently of head boundaries.
4. **Heterogeneous masked hydrostatics.** Constant total head with pressures
   on both sides of saturation must remain at rest through four orders of K
   heterogeneity and a stair-step active mask. This checks gravity cancellation
   and active topology, not the geometric accuracy of real topography.
5. **Nonlinear closed redistribution.** A nonhydrostatic unsaturated field must
   change while preserving total mixed storage with sealed faces. This is a
   conservation check, not a quantitative infiltration-front benchmark.
6. **Linear residual and actual solver matrix.** Independently perturb the
   actual residual, compare its derivatives with the matrix used to solve it,
   and verify diffusion/storage structure. A nonlinear frozen-K Picard matrix
   is deliberately not claimed to equal a true Newton Jacobian.

## Space and time order experiments

For the linear-storage constant-K branch,

```text
C partial H/partial t - K laplacian(H) = f
P = sin(pi*x/Lx) sin(pi*y/Ly) sin(pi*z/Lz)
H_exact = 3 m + A(1+t/tau)P
f = A [C/tau + K(kx²+ky²+kz²)(1+t/tau)] P
A=0.2 m; C=0.25 1/m; K=0.05 m/s; tau=1 s; final t=0.3 s.
```

The pressure passed to the model is `psi=H-z`. All six face groups receive
analytical **total** head evaluated at face centers; the source is evaluated
at cell centers and multiplied by volume. Linear time dependence makes the
backward Euler temporal derivative exact, isolating spatial error without
using an artificially discretized source. These are mathematical test
coefficients, not Washington soil estimates.

The time experiment uses `H=3+0.3 exp(-2t)`, constant in space, with
`f=-C*2*0.3 exp(-2t)`. Its continuum/FV spatial Laplacian is exactly zero.
Using 4, 8 and 16 steps to `t=0.4 s` isolates first-order backward Euler error
without contaminating the rate by an unresolved spatial mode.

The committed run measures spatial L2/Linf orders approximately 2.010 and
2.003, and temporal L2 orders 0.940 and 0.968 (Linf 0.931 and 0.963).
The fine spatial RMS head error is approximately `1.72e-4 m` (Linf `4.88e-4 m`); exact measured
values, local residuals and source hash are in the JSON. These rates concern
the stated smooth problems and linear mode, not storm peak convergence.

## Nonlinear MMS: retained diagnostic and unresolved order

A stationary smooth unsaturated target is

```text
psi = -0.8 + 0.12 sin(pi*x/Lx) sin(pi*y/Ly) sin(pi*z/Lz)
H = psi + z
f = -K'(psi) [|grad psi|² + partial psi/partial z]
    -K(psi) laplacian(psi).
```

The `partial psi/partial z` term comes from gravity in
`div[-K(psi)(grad psi+e_z)]`; omitting it would change the manufactured PDE.
`K'` is obtained by the analytical VG–Mualem chain rule, independently of the
kernel's spatial operator. Since the continuum target is stationary, its time
derivative is zero. Starting from it, one actual implicit step measures
spatial error in maintaining that balance.

Both a 1 s trial and a 0.1 s trial are retained. On 5/9/17 grids the 1 s trial
has L2 orders about 1.81/1.93 and Linf 1.37/1.66. The shorter trial has L2
1.25/1.57 and Linf 0.75/1.18. Errors decrease, but asymptotic second-order
Linf behavior is **not established**. These results are diagnostic, have no
formal `passed` value, and are excluded from the required acceptance flag.

A smooth orthogonal scalar-K solution is expected to approach second-order
head accuracy under adequate refinement and boundary treatment. Current
cell-conductivity half-cell boundary fluxes can produce localized errors, and
a short numerical response can be underresolved relative to its diffusion
length. Those are hypotheses explaining the measured behavior, not a proof
of its cause. Directional refinement, finer asymptotic meshes and boundary
flux comparisons are needed before promoting this nonlinear order claim.
The report preserves the unfavorable short-duration results rather than
changing a tolerance to classify them as accepted.

## Still required for scientific use

Independent Richards-code or analytical infiltration/capillary-rise benchmarks
must cover dry fronts, saturation crossings, perched layers, strong contrasts
and moving water tables. Surface ponding, bank/stream exchange and runoff need
their own coupled benchmarks. Test masks and vertical/horizontal refinement
for a particular valley and test its hydrograph, profile pressures and volumes;
there is no inference from a smooth MMS to an unresolved event peak.

Verification also does not establish observation-operator correctness,
assimilation skill, ecological calibration, mechanical triggering/deformation,
or suitability for operational warnings. Those claims require their own
tests, parameter evidence, observations and uncertainty evaluation.
