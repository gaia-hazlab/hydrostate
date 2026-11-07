from __future__ import annotations

from dataclasses import dataclass
from datetime import date
from pathlib import Path
import calendar
import csv

from puget_modflow.config import MONTH_KEYS, RechargeConfig


@dataclass(frozen=True)
class MonthlyPrecipRecord:
    date: date
    precip_in: float


@dataclass(frozen=True)
class RechargeRecord:
    stress_period: int
    date: date
    recharge_ft_per_day: float
    precip_in: float
    fraction: float
    seasonal_factor: float


def monthly_precip_to_recharge_rate_ft_per_day(
    *,
    precip_in: float,
    year: int,
    month: int,
    fraction: float,
    seasonal_factor: float = 1.0,
) -> float:
    """Convert monthly precipitation depth to a mean daily recharge flux."""
    if precip_in < 0:
        raise ValueError("precipitation cannot be negative")
    if fraction < 0:
        raise ValueError("recharge fraction cannot be negative")
    days = calendar.monthrange(year, month)[1]
    recharge_ft = (precip_in / 12.0) * fraction * seasonal_factor
    return recharge_ft / days


def recharge_series_from_monthly_precip(
    records: list[MonthlyPrecipRecord],
    config: RechargeConfig,
) -> list[RechargeRecord]:
    fraction = _bounded_fraction(
        config.default_fraction,
        min_fraction=config.min_fraction,
        max_fraction=config.max_fraction,
    )
    output: list[RechargeRecord] = []
    for index, record in enumerate(records):
        month_key = MONTH_KEYS[record.date.month - 1]
        seasonal_factor = float(config.seasonal_factor.get(month_key, 1.0))
        recharge = monthly_precip_to_recharge_rate_ft_per_day(
            precip_in=record.precip_in,
            year=record.date.year,
            month=record.date.month,
            fraction=fraction,
            seasonal_factor=seasonal_factor,
        )
        output.append(
            RechargeRecord(
                stress_period=index,
                date=record.date,
                recharge_ft_per_day=recharge,
                precip_in=record.precip_in,
                fraction=fraction,
                seasonal_factor=seasonal_factor,
            )
        )
    return output


def read_monthly_precip_csv(path: str | Path) -> list[MonthlyPrecipRecord]:
    records: list[MonthlyPrecipRecord] = []
    with Path(path).open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError("precipitation CSV has no header")
        precip_field = _find_field(
            reader.fieldnames,
            ("precip_in", "precip_inches", "ppt_in", "precip_mm", "ppt_mm"),
        )
        date_field = _find_field(reader.fieldnames, ("date", "month"))
        for row in reader:
            raw_precip = float(row[precip_field])
            precip_in = raw_precip / 25.4 if precip_field.lower().endswith("_mm") else raw_precip
            records.append(
                MonthlyPrecipRecord(
                    date=_parse_record_date(row[date_field]),
                    precip_in=precip_in,
                )
            )
    return records


def write_recharge_timeseries(records: list[RechargeRecord], path: str | Path) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "stress_period",
        "date",
        "recharge_ft_per_day",
        "precip_in",
        "fraction",
        "seasonal_factor",
    ]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for record in records:
            writer.writerow(
                {
                    "stress_period": record.stress_period,
                    "date": record.date.isoformat(),
                    "recharge_ft_per_day": f"{record.recharge_ft_per_day:.12g}",
                    "precip_in": f"{record.precip_in:.12g}",
                    "fraction": f"{record.fraction:.12g}",
                    "seasonal_factor": f"{record.seasonal_factor:.12g}",
                }
            )
    return output


def read_recharge_timeseries(path: str | Path) -> list[RechargeRecord]:
    records: list[RechargeRecord] = []
    with Path(path).open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            records.append(
                RechargeRecord(
                    stress_period=int(row["stress_period"]),
                    date=date.fromisoformat(row["date"]),
                    recharge_ft_per_day=float(row["recharge_ft_per_day"]),
                    precip_in=float(row["precip_in"]),
                    fraction=float(row["fraction"]),
                    seasonal_factor=float(row["seasonal_factor"]),
                )
            )
    return records


def constant_monthly_precip_records(
    *,
    start: date,
    end: date,
    precip_in: float,
) -> list[MonthlyPrecipRecord]:
    if start > end:
        raise ValueError("start must be before or equal to end")
    records: list[MonthlyPrecipRecord] = []
    year = start.year
    month = start.month
    while (year, month) <= (end.year, end.month):
        records.append(MonthlyPrecipRecord(date=date(year, month, 1), precip_in=precip_in))
        month += 1
        if month == 13:
            month = 1
            year += 1
    return records


def _bounded_fraction(value: float, *, min_fraction: float, max_fraction: float) -> float:
    return max(min_fraction, min(max_fraction, value))


def _find_field(fieldnames: list[str], candidates: tuple[str, ...]) -> str:
    normalized = {name.lower(): name for name in fieldnames}
    for candidate in candidates:
        if candidate in normalized:
            return normalized[candidate]
    raise ValueError(f"Missing one of fields: {', '.join(candidates)}")


def _parse_record_date(value: str) -> date:
    value = value.strip()
    if len(value) == 7:
        return date.fromisoformat(f"{value}-01")
    parsed = date.fromisoformat(value)
    return date(parsed.year, parsed.month, 1)
