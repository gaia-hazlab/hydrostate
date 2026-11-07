# Technical reports and verification evidence

This directory records the equations, numerical implementation and measured
verification scope of `hydrostate`. Passing a code-verification test does not
establish field accuracy, mesh convergence of the valley example, or skill at
predicting hazards.

- [Flow equations, discrete residual and solver matrix](flow_mathematics.md)
  derives the implemented operator and distinguishes modified Picard from a
  full Newton Jacobian.
- [Manufactured solutions and verification interpretation](verification.md)
  gives analytical sources, boundary values, norms, gates and remaining gaps.
- [Machine-readable measured verification evidence](verification_summary.json)
  contains errors/orders, required-check results, nonlinear diagnostics, source
  hash, dependency version and execution time.
- [3D valley resolution sensitivity](valley_refinement.json) records a separate
  24-hour time-step and staircase-grid probe. It is application sensitivity
  evidence, not another accepted manufactured-solution order test. Interpretation
  and reproduction commands are in the [use-case guide](../use_cases.md).
- [Conservative material changes and hydraulic restarts](material_updates.md)
  derives fixed-volume storage reconciliation and records rejection branches,
  complete hydraulic state and restart-equivalence checks.
- [Gaia model interchange](../datahub_contract.md) documents typed static and
  evolving properties, exact events, interval ledgers and full coupled-memory
  payloads verified through Gaia's actual storage APIs.

Regenerate the evidence against the current source:

```sh
PYTHONPATH=src python -m hydrostate.verification --output docs/technical/verification_summary.json
PYTHONPATH=src python -m unittest discover -s tests -p test_verification.py -v
```

Use the project Python environment with NumPy/SciPy. The runner exits nonzero
when a **required** acceptance check fails, and rejects nonfinite JSON. Its
overall `passed` flag covers the explicitly stated required scope. Nonlinear
MMS diagnostics have `passed: null`, `acceptance_required: false`, and report
that asymptotic maximum-norm accuracy is not established; they are not quietly
counted as successful verification.

Current evidence supports second-order spatial head convergence for smooth
constant-conductivity/linear-storage 3-D manufactured flow and first-order
backward Euler time convergence, plus constitutive, flux/sign, hydrostatic,
storage and linear-assembly checks. Smooth nonlinear MMS errors decrease, but
their maximum-norm convergence remains below the limiting second-order target
on the measured grids, especially at short response times. Sharp fronts,
contrasting layers, water-table crossings, terrain stair errors, surface
ponding and coupled DHSVM/MODFLOW behavior still require separate benchmarks.

The numerical valley and older hillslope toy are demonstrations. Neither their
budget closure nor these smooth-solution tests demonstrate convergence of a
storm outlet hydrograph. Scientific applications require application-specific
mesh/time refinement and independent observations/reference solutions.
