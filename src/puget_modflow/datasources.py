from __future__ import annotations

from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Iterable
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import csv
import json
import re
import shutil

from puget_modflow.config import PipelineConfig


PRISM_BASE_URL = "https://services.nacse.org/prism/data/get"
SCIENCEBASE_ITEM_URL = "https://www.sciencebase.gov/catalog/item"
NWIS_GW_LEVELS_URL = "https://waterservices.usgs.gov/nwis/gwlevels/"

VALID_PRISM_ELEMENTS = {"ppt", "tmin", "tmax", "tmean", "tdmean", "vpdmin", "vpdmax"}
VALID_PRISM_RESOLUTIONS = {"400m", "800m", "4km"}
VALID_PRISM_REGIONS = {"us", "ak", "hi", "pr"}


@dataclass(frozen=True)
class DownloadRecord:
    source: str
    url: str
    path: str


def prism_data_url(
    element: str,
    date_token: str | date,
    *,
    region: str = "us",
    resolution: str = "4km",
    fmt: str | None = None,
    long_term: bool = False,
) -> str:
    """Build a PRISM v2 web-service URL for one grid."""
    element = element.lower()
    region = region.lower()
    resolution = resolution.lower()
    token = _date_token(date_token)

    if element not in VALID_PRISM_ELEMENTS:
        raise ValueError(f"Unsupported PRISM element: {element}")
    if region not in VALID_PRISM_REGIONS:
        raise ValueError(f"Unsupported PRISM region: {region}")
    if resolution not in VALID_PRISM_RESOLUTIONS:
        raise ValueError(f"Unsupported PRISM resolution: {resolution}")

    parts = [PRISM_BASE_URL, region, resolution, element, token]
    if long_term:
        parts.append("lt")
    url = "/".join(parts)
    if fmt:
        url = f"{url}?{urlencode({'format': fmt})}"
    return url


def prism_release_date_url(
    element: str,
    date_token: str | date,
    *,
    region: str = "us",
    resolution: str = "4km",
    json_output: bool = True,
) -> str:
    token = _date_token(date_token)
    url = "/".join(
        [
            PRISM_BASE_URL,
            "releaseDate",
            region.lower(),
            resolution.lower(),
            element.lower(),
            token,
        ]
    )
    if json_output:
        url = f"{url}?{urlencode({'json': 'true'})}"
    return url


def month_tokens(start: str | date, end: str | date) -> list[str]:
    start_date = _parse_month(start)
    end_date = _parse_month(end)
    if start_date > end_date:
        raise ValueError("start must be before or equal to end")

    tokens: list[str] = []
    year = start_date.year
    month = start_date.month
    while (year, month) <= (end_date.year, end_date.month):
        tokens.append(f"{year:04d}{month:02d}")
        month += 1
        if month == 13:
            month = 1
            year += 1
    return tokens


def sciencebase_item_json_url(item_id: str) -> str:
    if not item_id:
        raise ValueError("ScienceBase item id is required")
    return f"{SCIENCEBASE_ITEM_URL}/{item_id}?format=json"


def sciencebase_children_json_url(item_id: str) -> str:
    if not item_id:
        raise ValueError("ScienceBase item id is required")
    return f"{SCIENCEBASE_ITEM_URL}/{item_id}/children?format=json"


def nwis_groundwater_levels_url(
    *,
    bbox: tuple[float, float, float, float] | None = None,
    start: str | date | None = None,
    end: str | date | None = None,
    site_status: str = "all",
    agency: str = "USGS",
    parameter_cd: str = "72019",
    output_format: str = "json",
) -> str:
    """Build an NWIS groundwater-level query URL.

    The default parameter code 72019 is depth to water level below land surface.
    """
    params: dict[str, str] = {
        "format": output_format,
        "siteType": "GW",
        "siteStatus": site_status,
        "agencyCd": agency,
        "parameterCd": parameter_cd,
    }
    if bbox:
        params["bBox"] = ",".join(f"{value:.8g}" for value in bbox)
    if start:
        params["startDT"] = _iso_date(start)
    if end:
        params["endDT"] = _iso_date(end)
    return f"{NWIS_GW_LEVELS_URL}?{urlencode(params)}"


def source_catalog(config: PipelineConfig) -> dict[str, object]:
    return {
        "references": [item.__dict__ for item in config.reference_datasets],
        "observations": [item.__dict__ for item in config.observation_sources],
        "precipitation": {
            "source": config.precipitation.source,
            "sample_url": prism_data_url(
                config.precipitation.element,
                config.model.start_date.strftime("%Y%m"),
                region=config.precipitation.region,
                resolution=config.precipitation.resolution,
            ),
        },
        "nwis_groundwater_levels_url": nwis_groundwater_levels_url(
            bbox=config.domain.bounds_wgs84,
            start=config.model.start_date,
            end=config.model.end_date,
        ),
    }


def fetch_json(url: str) -> dict[str, object]:
    request = Request(url, headers={"User-Agent": "puget-modflow/0.1"})
    with urlopen(request, timeout=60) as response:
        return json.loads(response.read().decode("utf-8"))


def download_url(url: str, output_dir: str | Path, filename: str | None = None) -> Path:
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    request = Request(url, headers={"User-Agent": "puget-modflow/0.1"})
    with urlopen(request, timeout=120) as response:
        resolved_name = filename or _filename_from_headers(response.headers) or _filename_from_url(url)
        destination = output_path / resolved_name
        with destination.open("wb") as handle:
            shutil.copyfileobj(response, handle)
    return destination


def download_prism_monthly(
    *,
    element: str,
    start: str | date,
    end: str | date,
    output_dir: str | Path,
    region: str = "us",
    resolution: str = "4km",
    fmt: str | None = None,
) -> list[DownloadRecord]:
    records: list[DownloadRecord] = []
    for token in month_tokens(start, end):
        url = prism_data_url(
            element,
            token,
            region=region,
            resolution=resolution,
            fmt=fmt,
        )
        path = download_url(url, output_dir)
        records.append(DownloadRecord(source="PRISM", url=url, path=str(path)))
    return records


def sciencebase_attached_files(item: dict[str, object]) -> list[dict[str, str]]:
    files = item.get("files", [])
    if not isinstance(files, list):
        return []
    records: list[dict[str, str]] = []
    for file_info in files:
        if not isinstance(file_info, dict):
            continue
        name = str(file_info.get("name") or "")
        url = str(file_info.get("downloadUri") or file_info.get("url") or "")
        if name and url:
            records.append({"name": name, "url": url})
    return records


def download_sciencebase_attached_files(
    item: dict[str, object],
    output_dir: str | Path,
) -> list[DownloadRecord]:
    records: list[DownloadRecord] = []
    for attached_file in sciencebase_attached_files(item):
        path = download_url(
            attached_file["url"],
            output_dir,
            filename=attached_file["name"],
        )
        records.append(
            DownloadRecord(
                source="ScienceBase",
                url=attached_file["url"],
                path=str(path),
            )
        )
    return records


def write_download_manifest(records: Iterable[DownloadRecord], path: str | Path) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["source", "url", "path"])
        writer.writeheader()
        for record in records:
            writer.writerow(record.__dict__)
    return output


def _date_token(value: str | date) -> str:
    if isinstance(value, date):
        return value.strftime("%Y%m%d")
    token = value.replace("-", "")
    if not re.fullmatch(r"\d{4}|\d{6}|\d{8}", token):
        raise ValueError(f"Date token must be YYYY, YYYYMM, or YYYYMMDD: {value}")
    return token


def _parse_month(value: str | date) -> date:
    if isinstance(value, date):
        return date(value.year, value.month, 1)
    token = value.strip()
    if re.fullmatch(r"\d{4}-\d{2}", token):
        year, month = token.split("-")
        return date(int(year), int(month), 1)
    if re.fullmatch(r"\d{6}", token):
        return date(int(token[:4]), int(token[4:]), 1)
    if re.fullmatch(r"\d{4}-\d{2}-\d{2}", token):
        parsed = date.fromisoformat(token)
        return date(parsed.year, parsed.month, 1)
    raise ValueError(f"Month must be YYYY-MM, YYYYMM, or ISO date: {value}")


def _iso_date(value: str | date) -> str:
    if isinstance(value, date):
        return value.isoformat()
    if re.fullmatch(r"\d{4}-\d{2}-\d{2}", value):
        return value
    if re.fullmatch(r"\d{4}-\d{2}", value):
        return f"{value}-01"
    raise ValueError(f"Date must be YYYY-MM-DD or YYYY-MM: {value}")


def _filename_from_headers(headers: object) -> str | None:
    content_disposition = headers.get("Content-Disposition")  # type: ignore[attr-defined]
    if not content_disposition:
        return None
    match = re.search(r'filename="?([^";]+)"?', content_disposition)
    return match.group(1) if match else None


def _filename_from_url(url: str) -> str:
    token = url.rstrip("/").split("/")[-1].split("?")[0]
    return f"{token}.download"
