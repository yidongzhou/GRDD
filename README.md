# Geodesic Regression Discontinuity Design (GRDD) — Replication Materials

Replication code for **Kurisu, Zhou, Otsu, and Müller (2025)**: simulations, bootstrap inference, a UK compositional election application, and a supplementary Taipei air-quality analysis. All scripts assume the working directory is this repository folder (`GRDD/`) unless `GRDD_ROOT` is set.

## Repository structure

```
GRDD/
├── README.md                 # This file
├── REPLICATION.md            # Map from paper figures/tables to commands
├── install_packages.R        # Install R package dependencies (CRAN)
├── figure1.R                 # Figure 1 — simulation plots (uses precomputed RData)
├── figure2.R                 # Figure 2 — UK election
├── figure3.R                 # Figure 3 — Taipei hourly CO (lines)
├── figure4.R                 # Figure 4 — Taipei 3D surfaces (HTML)
├── R/                        # GRDD implementation and local Fréchet helpers
│   ├── paths.R               # Project root (GRDD_ROOT or getwd())
│   ├── grdd.R
│   ├── grdd_inference.R
│   └── …                     # lfr_*.R, kerFctn.R, local_linear.R, lcm.R
├── scripts/
│   ├── simulation_estimation.R    # Long-run: network bias MC → output/rdata/estimation.RData
│   ├── simulation_inference.R     # Long-run: size/power MC → output/rdata/inference.RData
│   ├── application_uk_election.R  # UK analysis; loads uk.RData to skip grdd() when cached
│   └── supplement_taipei_metro.R # Taipei smoothing + GRDD + surfaces → output/rdata/taipei_metro.RData
├── data/                     # Raw inputs (CSVs, .dta)
│   ├── incumbency_advantage_*.csv
│   └── all_full_report_stations_94_07_hourly.dta  # Large; required for Taipei supplement
└── output/
    ├── figures/              # Generated PDFs and HTML
    └── rdata/                # estimation.RData, inference.RData, uk.RData, taipei_metro.RData
```

## Software requirements

- **R** (4.0 or later recommended)
- **R packages** (CRAN): see `install_packages.R` for the full list. Core usage draws on `ggplot2`, `rdrobust`, `truncnorm`, `foreach`, `doSNOW`, `haven`, `plotly`, `fdapace`, `pracma`, `manifold`, `Matrix`, `osqp`, `trust`, and dependencies.
- **Parallel simulations**: `scripts/simulation_*.R` use `parallel` + `doSNOW`; worker count defaults to `min(50, detectCores() - 1)`.

### Install packages

From the `GRDD/` directory:

```bash
Rscript install_packages.R
```

Record versions for the journal archive:

```r
sessionInfo()
```

## Reproducing results (quick reference)

| Goal | Command |
|------|---------|
| Figure 1 (a)–(b) | `Rscript figure1.R` |
| Figure 2 | `Rscript figure2.R` |
| Figure 3 | `Rscript figure3.R` (runs `supplement_taipei_metro.R` if needed) |
| Figure 4 | `Rscript figure4.R` |
| Re-run network bias simulation | `Rscript scripts/simulation_estimation.R` |
| Re-run inference simulation | `Rscript scripts/simulation_inference.R` |
| UK analysis only | `Rscript -e 'source("scripts/application_uk_election.R")'` |
| Taipei analysis only | `Rscript scripts/supplement_taipei_metro.R` |

**Outputs**

- Figures: `output/figures/` (PDF for Figures 1–3; HTML widgets for Figure 4).
- UK cache: `output/rdata/uk.RData` (written by `scripts/application_uk_election.R`).
- Taipei supplement: `output/rdata/taipei_metro.RData` (written by `scripts/supplement_taipei_metro.R`).

## Runtime expectations

- **Figure 1**: seconds (loads precomputed `output/rdata/estimation.RData` and `inference.RData`).
- **Figure 2**: with `uk.RData` present, skips the slow `grdd()` fit; inference, `rdrobust`, and local fits still run.
- **Figures 3–4**: long on first run (smoothing, functional GRDD, bootstrap, plotly); faster once `taipei_metro.RData` exists.
- **`simulation_estimation.R` / `simulation_inference.R`**: many hours (Monte Carlo × bootstrap; parallelization helps).

## Environment variables (optional)

| Variable | Effect |
|----------|--------|
| `GRDD_ROOT` | Absolute path to this repo if you cannot `setwd()` to `GRDD/` |
| `GRDD_FORCE_REFIT_UK` | Set to `1` to ignore cached objects in `output/rdata/uk.RData` and re-fit UK GRDD (when sourcing `application_uk_election.R`) |
| `GRDD_FORCE_REFIT_TAIPEI` | Set to `1` to rebuild `output/rdata/taipei_metro.RData` when running `figure3.R` / `figure4.R` |

## Data sources

Replication inputs follow the original public releases; place files in `data/` as in this repository.

- **United Kingdom (Figure 2).** Constituency-level vote shares and incumbency-advantage measures from Eggers & Spirling (2017): *Incumbency Effects and the Strength of Party Preferences: Evidence from Multiparty Elections in the United Kingdom*, *The Journal of Politics* 79, 903–920. CSV filenames match their replication materials.

- **Taipei air quality (Figures 3–4).** Hourly station CO (`all_full_report_stations_94_07_hourly.dta`, large file) from Chen & Whalley (2012): *Green Infrastructure: The Effects of Urban Rail Transit on Air Quality*, *American Economic Journal: Economic Policy* 4, 58–97. Data are used here for the supplementary functional-data illustration.

## Citation

If you use this code or methods, please cite the paper:

> Kurisu, D., Zhou, Y., Otsu, T., and Müller, H. G. (2025). Regression discontinuity designs for functional data and random objects in geodesic spaces. *arXiv* preprint arXiv:2506.18136. https://arxiv.org/abs/2506.18136

BibTeX (abbreviated):

```bibtex
@article{kurisu2025grdd,
  title   = {Regression discontinuity designs for functional data and random objects in geodesic spaces},
  author  = {Kurisu, D. and Zhou, Y. and Otsu, T. and M{\"u}ller, H.-G.},
  journal = {arXiv preprint arXiv:2506.18136},
  year    = {2025}
}
```

## License

This repository is released under the **MIT License**; see the [`LICENSE`](LICENSE) file in the root of the repository.
