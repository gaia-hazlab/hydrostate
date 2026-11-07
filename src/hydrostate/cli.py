"""Lazy command dispatch for the Hydrostate research workflows."""

from __future__ import annotations

import argparse
import sys


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    commands = {
        "valley": "run the synthetic 3D valley experiment",
        "verify": "run the mathematical verification experiments",
        "sensors": "ingest and validate sensor observations",
    }
    if args and args[0] in commands:
        command, rest = args[0], args[1:]
        if command == "valley":
            from .valley import main as entrypoint
        elif command == "verify":
            from .verification import main as entrypoint
        else:
            from .observations import main as entrypoint
        return int(entrypoint(rest) or 0)
    parser = argparse.ArgumentParser(
        prog="hydrostate",
        description="Watershed vadose/groundwater research, synthetic experiments and verification.",
    )
    parser.add_argument("command", nargs="?", choices=commands, help="; ".join(f"{k}: {v}" for k, v in commands.items()))
    parser.parse_args(args)
    parser.print_help()
    return 0
