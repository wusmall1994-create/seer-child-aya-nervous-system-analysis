# SEER childhood and AYA subsequent nervous system tumor analysis

This repository contains the R analysis code for a population-based SEER
multiple-primary standardized incidence ratio (MP-SIR) study of malignant brain
and other nervous system tumors occurring as the first subsequent malignancy
after cancer diagnosed before age 40 years.

## Repository scope

Only analysis code and documentation are included. This repository does **not**
contain:

- SEER data or SEER*Stat exports;
- individual-level or case-listing records;
- derived tables, figures, spreadsheets, or logs;
- manuscript or submission files; or
- credentials, local paths, or personal information.

SEER Research Data are subject to the SEER Research Data Agreement and must be
obtained directly from the National Cancer Institute. Do not redistribute SEER
data through this repository.

## Requirements

- R 4.6.1 or a compatible recent R release
- R packages: `ggplot2`, `patchwork`, `svglite`, `data.table`, and `survival`
- SEER*Stat 9.0.43 for producing the required input exports

Install the R dependencies with:

```r
install.packages(c("ggplot2", "patchwork", "svglite", "data.table", "survival"))
```

## Local directory structure

Create the following directories next to the analysis script:

```text
project/
├── 02_raw_exports/   # local SEER*Stat exports; never commit
├── 03_analysis/
│   └── run_all.R
├── 04_results/       # generated locally; never commit
└── 05_logs/          # generated locally; never commit
```

The expected input filenames are declared in the `input_files` vector near the
top of `run_all.R`. The code expects the column names produced by the study's
SEER*Stat MP-SIR and case-listing sessions.

## Run

From a terminal:

```sh
Rscript 03_analysis/run_all.R
```

By default, the project root is inferred from the location of `run_all.R`. For
testing an unchanged code checkout against an existing local project, set the
`SEER_PROJECT_DIR` environment variable to that project root. Never commit the
local value.

The script creates aggregate tables, figure files, source-data tables, and QA
logs under `04_results/` and `05_logs/`.

## September 2026 revision analyses

Run the original pipeline first to produce the aggregate input tables, followed by:

```sh
Rscript 03_analysis/revision_analysis.R
Rscript 03_analysis/cif_competing_risk.R
```

Both scripts support `SEER_PROJECT_DIR`. Optional output overrides are
`SEER_REVISION_OUTPUT_DIR` and `SEER_CIF_OUTPUT_DIR`.

`revision_analysis.R` implements the shared 15-category first-cancer model,
1,999-replicate parametric bootstrap (seed 20260921), 200 rounding scenarios
(seed 20260922), and the shared-category forest figure. It reads the original
aggregate tables and the local 05 first-cancer export.

`cif_competing_risk.R` reads the local full-cohort export
`15_full_cohort_index_competing_risk.txt`. Required columns include Patient ID,
Event Number, Reason for Exit, age recode, Survival months, the event-site
variable, Months Since Index (Calculated), and Person Time Years (Calculated).
The session uses a 2-month exclusion and exits at the next malignant tumor.
The analysis uses Aalen–Johansen estimates with death and other first subsequent
malignancies as competing events; study end and loss to follow-up are censoring.
There are 275,732 included persons and 345 eligible target events: one at exactly
2 months and 344 later. Eligible zero-time exits are retained with a numerical
epsilon of 1e-8 years. These counts do not establish the absence of ineligible
tumors during the initial exclusion window. A separate direct recursion checks
the cumulative-incidence point estimates against `survival::survfit`.

Only aggregate results are written by the release scripts. The local author
version's patient-level export has been removed from the release version.
The published scripts retain the same statistical calculations.

## Privacy safeguards

The event case listing is read locally only. Patient identifiers are used in
memory for record linkage and validation. The export functions stop if a
`Patient ID` column is present, and the pipeline scans generated CSV headers
before completion. No individual-level data should be committed or shared.

The `.gitignore` blocks common SEER, data, results, manuscript, credential, and
case-listing file types. Review `git status` and the staged diff before every
commit; `.gitignore` is an additional safeguard, not a substitute for review.

## Reproducibility note

The repository does not include SEER*Stat session files because they may contain
local database paths and other environment-specific settings. Users must
configure equivalent MP-SIR sessions under their own SEER data agreement.

## License

No license has been assigned. All rights are reserved unless the repository
owner adds a license.
