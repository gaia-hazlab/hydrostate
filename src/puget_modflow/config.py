from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import date, datetime
from pathlib import Path
from typing import Any
import tomllib


DEFAULT_CONFIG_TOML = """[project]
name = "puget_sound_groundwater"
workspace = "runs/puget_sound"

[domain]
name = "southeast_puget_sound"
description = "USGS active model area in Pierce and southern King Counties; expand only after building a new calibrated regional grid."
bounds_wgs84 = [-123.30, 46.75, -121.70, 47.55]
model_crs = "EPSG:2927"

[model]
engine = "modflow6_template"
reference_model = "usgs_southeast_puget_sound_modflow_nwt"
length_units = "feet"
time_units = "days"
start_date = "2005-01-01"
end_date = "2015-12-31"
nlay = 13
nrow = 32
ncol = 32
delr = 500.0
delc = 500.0
top = 300.0
bottom = -700.0

[precipitation]
source = "PRISM"
region = "us"
resolution = "4km"
element = "ppt"
time_step = "monthly"
units = "millimeters"

[recharge]
method = "swb_preferred_precip_fraction_fallback"
default_fraction = 0.18
min_fraction = 0.02
max_fraction = 0.45
units = "feet_per_day"

[recharge.seasonal_factor]
jan = 1.20
feb = 1.15
mar = 1.05
apr = 0.95
may = 0.75
jun = 0.55
jul = 0.35
aug = 0.30
sep = 0.45
oct = 0.85
nov = 1.10
dec = 1.25

[[reference_datasets]]
name = "USGS 2024 spatial data for southeast Puget Sound groundwater characterization"
doi = "10.5066/P9JFKLMG"
sciencebase_item_id = "62fa7dfed34e25bd092879e6"
role = "GIS framework, hydrogeologic surfaces, observations, and model support data"

[[reference_datasets]]
name = "USGS 2023 MODFLOW-NWT groundwater-flow model near Puget Sound"
doi = "10.5066/P9LU1PMQ"
sciencebase_item_id = "62d823d6d34e83b67d134f25"
role = "Published calibrated steady-state and transient model inputs and outputs"

[[reference_datasets]]
name = "USGS 2023 Soil-Water Balance model for Puyallup and Chambers-Clover basins"
doi = "10.5066/P925FHGQ"
sciencebase_item_id = "5e010753e4b0b207aa033db2"
role = "Preferred recharge estimation workflow for the active model area"

[[observation_sources]]
name = "USGS NWIS groundwater levels"
url = "https://waterservices.usgs.gov/nwis/gwlevels/"
role = "Historical groundwater-level observations for calibration and validation"

[[observation_sources]]
name = "Washington Ecology EIM and well-log/GIS data"
url = "https://ecology.wa.gov/water-shorelines/water-quality/groundwater/groundwater-resources"
role = "State groundwater studies, well logs, and GIS well datasets"
"""


MONTH_KEYS = (
    "jan",
    "feb",
    "mar",
    "apr",
    "may",
    "jun",
    "jul",
    "aug",
    "sep",
    "oct",
    "nov",
    "dec",
)


@dataclass(frozen=True)
class ProjectConfig:
    name: str = "puget_sound_groundwater"
    workspace: str = "runs/puget_sound"


@dataclass(frozen=True)
class DomainConfig:
    name: str
    description: str
    bounds_wgs84: tuple[float, float, float, float]
    model_crs: str


@dataclass(frozen=True)
class ModelConfig:
    engine: str
    reference_model: str
    length_units: str
    time_units: str
    start_date: date
    end_date: date
    nlay: int
    nrow: int
    ncol: int
    delr: float
    delc: float
    top: float
    bottom: float


@dataclass(frozen=True)
class PrecipitationConfig:
    source: str
    region: str
    resolution: str
    element: str
    time_step: str
    units: str


@dataclass(frozen=True)
class RechargeConfig:
    method: str
    default_fraction: float
    min_fraction: float
    max_fraction: float
    units: str
    seasonal_factor: dict[str, float] = field(default_factory=dict)


@dataclass(frozen=True)
class ReferenceDataset:
    name: str
    doi: str
    sciencebase_item_id: str
    role: str


@dataclass(frozen=True)
class ObservationSource:
    name: str
    url: str
    role: str


@dataclass(frozen=True)
class PipelineConfig:
    project: ProjectConfig
    domain: DomainConfig
    model: ModelConfig
    precipitation: PrecipitationConfig
    recharge: RechargeConfig
    reference_datasets: tuple[ReferenceDataset, ...] = ()
    observation_sources: tuple[ObservationSource, ...] = ()

    @classmethod
    def from_toml(cls, path: str | Path) -> "PipelineConfig":
        with Path(path).open("rb") as handle:
            data = tomllib.load(handle)
        return cls.from_mapping(data)

    @classmethod
    def from_mapping(cls, data: dict[str, Any]) -> "PipelineConfig":
        project = ProjectConfig(**data.get("project", {}))

        domain_data = data["domain"]
        domain = DomainConfig(
            name=domain_data["name"],
            description=domain_data.get("description", ""),
            bounds_wgs84=tuple(float(v) for v in domain_data["bounds_wgs84"]),
            model_crs=domain_data["model_crs"],
        )

        model_data = data["model"]
        model = ModelConfig(
            engine=model_data["engine"],
            reference_model=model_data["reference_model"],
            length_units=model_data["length_units"],
            time_units=model_data["time_units"],
            start_date=_parse_date(model_data["start_date"]),
            end_date=_parse_date(model_data["end_date"]),
            nlay=int(model_data["nlay"]),
            nrow=int(model_data["nrow"]),
            ncol=int(model_data["ncol"]),
            delr=float(model_data["delr"]),
            delc=float(model_data["delc"]),
            top=float(model_data["top"]),
            bottom=float(model_data["bottom"]),
        )

        precip = PrecipitationConfig(**data["precipitation"])

        recharge_data = dict(data["recharge"])
        seasonal_factor = {
            month: float(recharge_data.get("seasonal_factor", {}).get(month, 1.0))
            for month in MONTH_KEYS
        }
        recharge_data.pop("seasonal_factor", None)
        recharge = RechargeConfig(
            method=recharge_data["method"],
            default_fraction=float(recharge_data["default_fraction"]),
            min_fraction=float(recharge_data["min_fraction"]),
            max_fraction=float(recharge_data["max_fraction"]),
            units=recharge_data["units"],
            seasonal_factor=seasonal_factor,
        )

        references = tuple(
            ReferenceDataset(**item) for item in data.get("reference_datasets", [])
        )
        observation_sources = tuple(
            ObservationSource(**item) for item in data.get("observation_sources", [])
        )
        return cls(
            project=project,
            domain=domain,
            model=model,
            precipitation=precip,
            recharge=recharge,
            reference_datasets=references,
            observation_sources=observation_sources,
        )

    def to_dict(self) -> dict[str, Any]:
        result = asdict(self)
        result["model"]["start_date"] = self.model.start_date.isoformat()
        result["model"]["end_date"] = self.model.end_date.isoformat()
        return result


def write_default_config(path: str | Path, overwrite: bool = False) -> Path:
    output = Path(path)
    if output.exists() and not overwrite:
        raise FileExistsError(f"{output} already exists; pass overwrite=True to replace it")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(DEFAULT_CONFIG_TOML, encoding="utf-8")
    return output


def _parse_date(value: Any) -> date:
    if isinstance(value, date):
        return value
    if isinstance(value, datetime):
        return value.date()
    if isinstance(value, str):
        return date.fromisoformat(value)
    raise TypeError(f"Unsupported date value: {value!r}")
