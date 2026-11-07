# Data Sources

This pipeline is designed around public, authoritative datasets that can be
downloaded and versioned outside the repository.

## Hydrogeologic Framework and MODFLOW Reference

The strongest public basis for a calibrated Puget Sound groundwater model is
the USGS southeast Puget Sound work in King and Pierce Counties.

- USGS conceptual framework and groundwater budget, 2024:
  `https://doi.org/10.3133/sir20245026v1`
- USGS numerical groundwater-flow model report, 2024:
  `https://doi.org/10.3133/sir20245026v2`
- USGS spatial data release, 2024:
  `https://doi.org/10.5066/P9JFKLMG`
- USGS MODFLOW-NWT model data release, 2023:
  `https://doi.org/10.5066/P9LU1PMQ`
- USGS Soil-Water Balance model data release, 2023:
  `https://doi.org/10.5066/P925FHGQ`

Key model facts from the 2024 USGS numerical report:

- Active model area: about 887 square miles in Pierce and King Counties.
- Time period: steady-state average for 2005-01-01 through 2015-12-31, plus
  monthly transient simulation over the same period.
- Grid: 13 layers, 416 rows, 433 columns, 500 ft by 500 ft cells.
- Calibration targets: measured heads in wells and lakes, plus estimated
  baseflow at selected streamgages.
- Important uncertainty sources: recharge, baseflow, hydraulic properties,
  water use, and the three-dimensional hydrogeologic structure.

For broader Puget Sound exploratory work, use the regional USGS framework:

- USGS Professional Paper 1424-D, Hydrogeologic framework of the Puget Sound
  aquifer system:
  `https://doi.org/10.3133/pp1424D`

That regional framework groups the Puget Sound aquifer system into regional
hydrogeologic units such as alluvial valley aquifers, surficial semiconfining
units, Fraser aquifer units, confining units, and Puget aquifer units. It is a
coarser conceptual framework than the 2024 southeast Puget Sound model and must
be recalibrated if used for prediction.

## Precipitation Forcing

PRISM Climate Group provides gridded precipitation fields suitable for monthly
and daily recharge forcing:

- Dataset documentation:
  `https://data.prism.oregonstate.edu/PRISM_datasets.pdf`
- Web service documentation:
  `https://prism.oregonstate.edu/documents/PRISM_downloads_web_service.pdf`
- Service pattern:
  `https://services.nacse.org/prism/data/get/us/4km/ppt/YYYYMM`

The PRISM service returns one grid per request. Current documentation states
that monthly data are available from 1895 to present and daily data from 1981 to
present for time-series variables including precipitation.

## Groundwater-Level Observations

USGS NWIS groundwater level service:

- Documentation:
  `https://nwis.waterservices.usgs.gov/docs/groundwater-levels/groundwater-levels-details/`
- Service endpoint:
  `https://waterservices.usgs.gov/nwis/gwlevels/`

Washington Department of Ecology groundwater resources:

- Groundwater resources:
  `https://ecology.wa.gov/water-shorelines/water-quality/groundwater/groundwater-resources`

The Ecology page links to the Environmental Information Management groundwater
search, well-log search, and GIS well datasets. These are useful for building
local observation sets, checking well construction, and separating confined from
unconfined observation targets.

## Practical Notes

- Keep raw ScienceBase, PRISM, and NWIS downloads under `data/raw/`.
- Keep processed time series and gridded recharge products under
  `data/processed/`.
- Record all source URLs and download dates in run manifests.
- Do not treat precipitation as recharge without calibration. In the Puget
  Lowland, evapotranspiration, land cover, soil storage, stormwater routing,
  glacial till, outwash, and perched systems all affect recharge timing and
  magnitude.
