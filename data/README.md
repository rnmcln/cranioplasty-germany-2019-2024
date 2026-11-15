# Input data

The analysis uses aggregate counts and summary statistics retrieved from the InEK
DatenBrowser (https://datenbrowser.inek.org; free registration required). The
retrieved values are not redistributed in this repository. To reproduce the analysis,
retrieve them as described below and save them as CSV files in `data/input/`, using
the column layouts in `data/templates/`.

## Data source

For each data year 2019–2024, the annual DRG data delivery under §21 KHEntgG grouped
to the following year's aG-DRG classification was used (for example, "Datenlieferung
DRG 2023 gruppiert nach 2024"). Values were retrieved on 2–3 July 2026. Counts below the
minimum cell size of the DatenBrowser are suppressed at source and were not imputed.

## Selections

Procedure codes were selected in any procedure position. A selection with several
codes returns episodes carrying at least one of them.

| Label in input files | OPS codes selected |
|---|---|
| `Autologous` | 5-020.2 |
| `Broad alloplastic` (alloplastic total) | 5-020.70, 5-020.71, 5-020.72, 5-020.74, 5-020.75 |
| `CAD` | 5-020.71, 5-020.72, 5-020.74, 5-020.75 |

The simple alloplastic group is computed in `code/01_analysis.py` as alloplastic total
minus CAD. Subcodes 5-020.73, .76 and .77 were not populated in any data year.

Filters added to each procedure selection:

- in-hospital death: discharge reason "07 Tod";
- coded infection: secondary diagnosis T81.4 or T85.72 (ICD-10-GM); T85.81 was also
  queried and was suppressed in every group and year.

## Files

| File | Content | Source of values |
|---|---|---|
| `mortality.csv` | year, group, OPS codes, discharge episodes, in-hospital deaths, mortality (%) | summary block; discharge-reason filter 07 |
| `complications.csv` | year, group, discharge episodes, episodes with T85.72, episodes with T81.4 | case count with the secondary-diagnosis filter |
| `demographics.csv` | year, group, episodes, sex (%), mean and SD length of stay, short/normal/long-stay (%), PCCL 0 (%), adults (%) | summary block; adults = sum of the age bands 18–29 years to 80+ years |
| `secondary_diagnoses_autologous_2023.csv` | eight most frequent secondary diagnoses, autologous group, 2023 | secondary-diagnosis output |
| `drg_mix_autologous_2023.csv` | autologous episodes in DRG families B20 and B02, 2023 | case count with DRG filter |
| `drg_cost_b20.csv` | mean per-case cost, cases and mean length of stay for DRGs B20A–B20E, main departments, 2022–2024; case-weighted family mean | cost report browser (whole-DRG values) |

The case-weighted B20 family mean weights each DRG by its total number of cases. Cost
values are reported for whole DRGs and are not specific to cranioplasty.
