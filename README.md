# Block-Calibrated Shape-Aware Screening of Cryptocurrency Volatility Episodes

R implementation of the curvature-informed risk diagnostic (CIRD) for Bitcoin, Ethereum, XRP and Litecoin. Authors: Mervenur Sözen and Çağlar Sözen.

The analysis uses cumulative forward realized-variance curves and is **ex post**: a full curve is observable after its forward horizon has elapsed. The baseline uses a 30-day horizon, seven-day origin blocks, and chronological training/calibration/test splits with 30-origin embargoes.

## Reproduce the analysis

From the repository root, with R installed:

```r
source("install_dependencies.R")
source("cird_analysis.R")
```

Or from a terminal:

```sh
Rscript install_dependencies.R
Rscript cird_analysis.R
Rscript verify_results.R results reported
```

The default analysis checks the four frozen Yahoo Finance CSV inputs against `data/checksums_md5.csv`. It uses those inputs without a market-data download and writes to `results/`. For a separate output directory, set the environment variable `CIRD_OUTPUT_DIR` before running. The script records the R session and checksums for each run. Package installation requires internet access; the analysis itself uses the supplied frozen inputs.

## Contents

- `cird_analysis.R`: complete empirical analysis, benchmarks, component ablations, weight and design sensitivity, Monte Carlo experiments, random-score diagnostic, dated reference mapping, tables and figures.
- `data/`: four frozen price files covering 2019–2024, with MD5 checksums.
- `install_dependencies.R`: explicit installation of required packages.
- `verify_results.R`: numerical output comparison and checks of all eight plot-data sets and PNG dimensions.
- `reported/tables/`: numerical reference outputs; the calendar mapping uses full March and November calendar months.
- `reported/audit/`: reporting-session information and documented reproduction checks.

The script generates four main tables, fourteen Appendix tables, and eight figures. It also exports full-precision numerical results, fitted objects and block calendars. Figure bytes can differ across R graphics devices and fonts; the verification script checks their source data and dimensions rather than requiring identical PNG hashes.

## Interpretation

Split-conformal calibration bounds block exceedance under exchangeability of block scores. Chronological cryptocurrency blocks do not establish that assumption. Observed calibration-gap and weight rankings are descriptive. A small gap can occur for an uninformative score and does not establish screening skill. The one-sided rank benchmark uses the same upper-tail direction as the range reference. Level-only CIRD and terminal FRV have identical rank-based alarm sets. The high–low range reference is a descriptive evaluation aid with overlapping windows and repeated peaks; it does not enter fitting or calibration.

Monte Carlo experiments evaluate calibration, including a volatility-level break. They do not supply detection labels or establish detection superiority. Dated episode mappings describe forward-window overlap and reference peak dates, with no claim of advance warning or causation.

## License and data source

The R code is covered by the repository's MIT license. The frozen data originate from Yahoo Finance; that code license does not grant additional rights in source-provider data. Cite the associated article when using the analysis.
