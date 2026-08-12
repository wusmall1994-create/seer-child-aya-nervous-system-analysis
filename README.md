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
- R packages: `ggplot2`, `patchwork`, and `svglite`
- SEER*Stat 9.0.43 for producing the required input exports

Install the R dependencies with:

```r
install.packages(c("ggplot2", "patchwork", "svglite"))
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
