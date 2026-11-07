#!/usr/bin/env python3
"""Ingest CSV sensors and compare them to precomputed hydrostate model outputs."""

from __future__ import annotations

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from hydrostate.observations import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())
