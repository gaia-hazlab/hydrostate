from __future__ import annotations

from pathlib import Path
import argparse
import json
import sys

from puget_modflow.config import PipelineConfig, write_default_config
from puget_modflow.datasources import (
    download_prism_monthly,
    download_sciencebase_attached_files,
    fetch_json,
    sciencebase_children_json_url,
    sciencebase_item_json_url,
    source_catalog,
    write_download_manifest,
)
from puget_modflow.model import build_mf6_template, run_modflow_workspace
from puget_modflow.pipeline import write_pipeline_plan
from puget_modflow.precip import aggregate_prism_to_domain_csv
from puget_modflow.recharge import (
    read_monthly_precip_csv,
    recharge_series_from_monthly_precip,
    write_recharge_timeseries,
)


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="puget-modflow",
        description="Puget Sound MODFLOW pipeline utilities.",
    )
    subparsers = parser.add_subparsers(required=True)

    init = subparsers.add_parser("init-config", help="write the default TOML config")
    init.add_argument("--output", default="config/puget_sound.toml")
    init.add_argument("--overwrite", action="store_true")
    init.set_defaults(func=cmd_init_config)

    plan = subparsers.add_parser("plan", help="write a JSON pipeline execution plan")
    add_config_arg(plan)
    plan.add_argument("--output", default="runs/puget_sound/pipeline_plan.json")
    plan.set_defaults(func=cmd_plan)

    sources = subparsers.add_parser("sources", help="print configured source URLs")
    add_config_arg(sources)
    sources.add_argument("--output")
    sources.set_defaults(func=cmd_sources)

    fetch_ref = subparsers.add_parser(
        "fetch-reference",
        help="fetch ScienceBase metadata for configured USGS reference datasets",
    )
    add_config_arg(fetch_ref)
    fetch_ref.add_argument("--output", default="data/raw/usgs_reference")
    fetch_ref.add_argument("--children", action="store_true", help="also fetch child item metadata")
    fetch_ref.add_argument(
        "--download-files",
        action="store_true",
        help="download files attached directly to configured ScienceBase items",
    )
    fetch_ref.set_defaults(func=cmd_fetch_reference)

    fetch_prism = subparsers.add_parser("fetch-prism", help="download monthly PRISM grids")
    add_config_arg(fetch_prism)
    fetch_prism.add_argument("--start", required=True, help="YYYY-MM")
    fetch_prism.add_argument("--end", required=True, help="YYYY-MM")
    fetch_prism.add_argument("--output", default="data/raw/prism/monthly_ppt")
    fetch_prism.add_argument("--format", choices=["nc", "asc", "bil"], default=None)
    fetch_prism.set_defaults(func=cmd_fetch_prism)

    aggregate = subparsers.add_parser(
        "aggregate-prism",
        help="aggregate downloaded PRISM rasters to a domain-mean monthly CSV",
    )
    add_config_arg(aggregate)
    aggregate.add_argument("--input", required=True)
    aggregate.add_argument("--output", default="data/processed/monthly_precip.csv")
    aggregate.set_defaults(func=cmd_aggregate_prism)

    recharge = subparsers.add_parser(
        "recharge-from-csv",
        help="convert monthly precipitation CSV to recharge time series",
    )
    add_config_arg(recharge)
    recharge.add_argument("--precip-csv", required=True)
    recharge.add_argument("--output", default="data/processed/recharge_timeseries.csv")
    recharge.set_defaults(func=cmd_recharge_from_csv)

    build = subparsers.add_parser(
        "build-template",
        help="build a small MODFLOW 6 template model for smoke testing",
    )
    add_config_arg(build)
    build.add_argument("--recharge")
    build.add_argument("--output", default="runs/template_mf6")
    build.add_argument("--exe", default="mf6")
    build.set_defaults(func=cmd_build_template)

    run = subparsers.add_parser("run", help="run MODFLOW in an existing workspace")
    run.add_argument("--workspace", required=True)
    run.add_argument("--exe", default="mf6")
    run.add_argument("--timeout", type=int)
    run.set_defaults(func=cmd_run)

    return parser


def add_config_arg(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--config", default="config/puget_sound.toml")


def cmd_init_config(args: argparse.Namespace) -> int:
    path = write_default_config(args.output, overwrite=args.overwrite)
    print(path)
    return 0


def cmd_plan(args: argparse.Namespace) -> int:
    config = PipelineConfig.from_toml(args.config)
    path = write_pipeline_plan(config, args.output)
    print(path)
    return 0


def cmd_sources(args: argparse.Namespace) -> int:
    config = PipelineConfig.from_toml(args.config)
    payload = source_catalog(config)
    text = json.dumps(payload, indent=2)
    if args.output:
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(text, encoding="utf-8")
        print(output)
    else:
        print(text)
    return 0


def cmd_fetch_reference(args: argparse.Namespace) -> int:
    config = PipelineConfig.from_toml(args.config)
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    written = []
    downloads = []
    for dataset in config.reference_datasets:
        if not dataset.sciencebase_item_id:
            continue
        item_id = dataset.sciencebase_item_id
        item = fetch_json(sciencebase_item_json_url(item_id))
        item_path = output / f"{item_id}.json"
        item_path.write_text(json.dumps(item, indent=2), encoding="utf-8")
        written.append(item_path)
        if args.download_files:
            downloads.extend(download_sciencebase_attached_files(item, output / item_id))
        if args.children:
            children = fetch_json(sciencebase_children_json_url(item_id))
            child_path = output / f"{item_id}.children.json"
            child_path.write_text(json.dumps(children, indent=2), encoding="utf-8")
            written.append(child_path)
    if downloads:
        manifest = write_download_manifest(downloads, output / "sciencebase_download_manifest.csv")
        written.append(manifest)
    for path in written:
        print(path)
    if not written:
        print("No ScienceBase item ids configured; see docs/data_sources.md")
    return 0


def cmd_fetch_prism(args: argparse.Namespace) -> int:
    config = PipelineConfig.from_toml(args.config)
    records = download_prism_monthly(
        element=config.precipitation.element,
        start=args.start,
        end=args.end,
        output_dir=args.output,
        region=config.precipitation.region,
        resolution=config.precipitation.resolution,
        fmt=args.format,
    )
    manifest = write_download_manifest(records, Path(args.output) / "download_manifest.csv")
    print(manifest)
    return 0


def cmd_aggregate_prism(args: argparse.Namespace) -> int:
    config = PipelineConfig.from_toml(args.config)
    output = aggregate_prism_to_domain_csv(
        input_dir=args.input,
        output_csv=args.output,
        bounds_wgs84=config.domain.bounds_wgs84,
    )
    print(output)
    return 0


def cmd_recharge_from_csv(args: argparse.Namespace) -> int:
    config = PipelineConfig.from_toml(args.config)
    precip = read_monthly_precip_csv(args.precip_csv)
    recharge = recharge_series_from_monthly_precip(precip, config.recharge)
    output = write_recharge_timeseries(recharge, args.output)
    print(output)
    return 0


def cmd_build_template(args: argparse.Namespace) -> int:
    config = PipelineConfig.from_toml(args.config)
    workspace = build_mf6_template(
        config,
        args.output,
        recharge_csv=args.recharge,
        exe_name=args.exe,
    )
    print(workspace)
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    result = run_modflow_workspace(
        args.workspace,
        exe_name=args.exe,
        timeout=args.timeout,
    )
    if result.stdout:
        print(result.stdout)
    if result.stderr:
        print(result.stderr, file=sys.stderr)
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
