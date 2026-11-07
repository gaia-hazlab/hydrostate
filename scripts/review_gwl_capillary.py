"""Reproduce a closed-transfer counterexample in a pinned external GWL checkout.

This is an audit diagnostic, not a test claiming the external model passes.
No upstream files or observational data are changed. Run with Hydrostate's
numerics extra and supply a trusted checkout of the reviewed public source.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import numpy as np

REVIEWED_COMMIT = '435445f17b90e1cdc4ee4e12a9e1cb17338f7e3f'


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checkout', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--initial-dtw-m', type=float, default=.5)
    args = parser.parse_args()
    root = args.checkout.resolve()
    commit = subprocess.check_output(['git', '-C', str(root), 'rev-parse', 'HEAD'], text=True).strip()
    if commit != REVIEWED_COMMIT:
        raise ValueError(f'Requires reviewed commit {REVIEWED_COMMIT}; got {commit}')
    source = root / 'src/models/water_budget.py'
    # HEAD alone is insufficient when the working source has local edits.
    expected = subprocess.check_output(['git', '-C', str(root), 'show', 'HEAD:src/models/water_budget.py'])
    if source.read_bytes() != expected:
        raise ValueError('Reviewed water-budget source has local edits')
    spec = importlib.util.spec_from_file_location('reviewed_gwl_water_budget', source)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    wp, fc, sat, depth, sy, dtw = .10, .28, .45, 1., .30, args.initial_dtw_m
    if not np.isfinite(dtw) or dtw <= 0:
        raise ValueError('initial DTW must be finite positive')
    result = module.coupled_water_budget(np.zeros((1, 1)), np.zeros((1, 1)), wp, fc, sat,
        root_depth_m=depth, specific_yield=sy, wt_depth0_m=dtw, init='wp', dt_days=1.)
    soil_delta = float((float(result.theta[0, 0]) - wp) * depth * 1000)
    gw_delta = float((dtw - float(result.wt_depth_m[0, 0])) * sy * 1000)
    external_out = float(result.runoff_mm[0, 0] + result.interflow_mm[0, 0]
                         + result.baseflow_mm[0, 0] + result.aet_mm[0, 0])
    residual = soil_delta + gw_delta + external_out
    record = {'scope': 'One closed root-zone/aquifer transfer; source audit, not field validation',
        'repository': 'https://github.com/gaia-hazlab/gwl-space-time-smooth', 'commit': commit,
        'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
        'inputs': {'precip_mm': 0., 'pet_mm': 0., 'theta_wp': wp, 'theta_fc': fc,
                   'theta_sat': sat, 'root_depth_m': depth, 'specific_yield': sy,
                   'water_table_depth_initial_m': dtw, 'init': 'wp', 'dt_days': 1.},
        'soil_storage_change_mm': soil_delta, 'aquifer_storage_change_mm': gw_delta,
        'external_input_mm': 0., 'external_output_mm': external_out,
        'recorded_internal_capillary_mm': float(result.cap_rise_mm[0, 0]),
        'combined_water_residual_mm': residual,
        'expected_combined_water_residual_mm': 0.,
        'absolute_tolerance_mm': 1e-4, 'combined_budget_passed': abs(residual) <= 1e-4,
        'interpretation': 'Capillary water is credited to soil but not debited from aquifer in this reviewed source. Float32 output explains small rounding differences. This counterexample does not assess other regimes or reported field skill.'}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2, allow_nan=False) + '\n')
    print(json.dumps(record, indent=2))
    return 0  # successful audit execution; model acceptance is explicitly recorded above


if __name__ == '__main__':
    raise SystemExit(main())
