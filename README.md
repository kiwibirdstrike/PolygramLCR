# PolygramLCR

This repository contains the source code and replication materials for logistic
polygram regression with three-part compositional predictors.

## Repository Layout

- `replication/polygramLCR/`: local R package source and numbered replication
  scripts.
- `replication/polygramLCR/scripts/`: scripts for simulation generation, model
  fitting, result summaries, and figure generation.
- `replication/polygramLCR/data/`: processed ATUS analysis data used by the
  real-data scripts.
- `replication/ATUS/`: ATUS public-use source files and documentation.
- `replication/outputs/`: generated outputs. This directory is ignored by Git
  because the full fitted outputs are large and can be regenerated from the
  scripts below.

## Requirements

The current replication environment used R 4.4.3. The package uses compiled C++
code through Rcpp, so macOS users should have Xcode Command Line Tools installed
and Windows users should have Rtools installed.

Install the required R packages:

```bash
Rscript -e 'install.packages(c("Rcpp", "RTriangle", "mgcv", "e1071"), repos = "https://cloud.r-project.org")'
```

## Install the Polygram Package

From the repository root, install the local package before running any
replication script:

```bash
R CMD INSTALL replication/polygramLCR
```

If reinstalling from a working tree that already has local compiled artifacts,
use `R CMD INSTALL --preclean replication/polygramLCR` to force a clean C++
rebuild.

Optional smoke test after installation:

```bash
Rscript replication/polygramLCR/tests/polygram-regression.R
```

## Reproduce the Analyses

Run the scripts from the package root, not from inside the `scripts/` directory.
The scripts use paths relative to `replication/polygramLCR`.

```bash
cd replication/polygramLCR

Rscript scripts/01_generate_simulation_data.R

Rscript scripts/02_1_fit_simulation_logistic.R
Rscript scripts/02_2_fit_simulation_gam.R
Rscript scripts/02_3_fit_simulation_svm.R
Rscript scripts/02_4_fit_simulation_polygram.R
Rscript scripts/03_make_simulation_results.R

Rscript scripts/04_1_fit_real_data_logistic.R
Rscript scripts/04_2_fit_real_data_gam.R
Rscript scripts/04_3_fit_real_data_svm.R
Rscript scripts/04_4_fit_real_data_polygram.R
Rscript scripts/05_make_real_data_results.R
```

Generated files are written under `replication/outputs/`:

- `replication/outputs/simulation/data/`: generated simulation train/test data.
- `replication/outputs/simulation/fits/`: fitted simulation models and
  prediction tables.
- `replication/outputs/simulation/results/`: simulation summaries and active
  edge tables.
- `replication/outputs/real_data/fits/`: fitted ATUS models and prediction
  tables.
- `replication/outputs/real_data/results/`: ATUS summary tables.
- `replication/outputs/figures/`: generated replication figures.

The simulation uses 30 replicates. SIM00, SIM01, and SIM02 use 5,000 training
observations and 5,000 test observations per replicate. SIM03 uses 10,000
training observations and 5,000 test observations per replicate.

## Notes for GitHub

The complete generated output directory is roughly 2 GB, and several generated
CSV files exceed GitHub's normal 100 MB file limit. For that reason,
`replication/outputs/` is excluded from Git and should be regenerated locally
with the commands above.

Manuscript files under the local `paper/` directory are also excluded from this
repository.
