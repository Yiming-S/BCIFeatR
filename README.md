# BCIFeatR

Feature Extraction Toolkit for EEG-based Brain–Computer Interfaces in R.

## Overview

BCIFeatR provides a unified train/test interface for extracting features from
multi-channel EEG trial data. It is designed for offline BCI decoding pipelines
where a consistent API across diverse feature families simplifies
experimentation.

### Supported feature methods

| Method | Description |
|--------|-------------|
| `logvar` | Log-variance of each channel |
| `logvar_pca` | Log-variance after PCA projection |
| `CSP` | Common Spatial Patterns |
| `FBCSP` | Filter-Bank CSP |
| `FBCSSP` | Filter-Bank CSSP (with time-delay embedding) |
| `TS` | Tangent-space projection from covariance matrices |
| `ACM_TS` | Augmented Covariance Matrix + tangent space |
| `Riemannian` | Riemannian mean + log-map with optional geodesic/FGDA filtering |
| `ATM` | Avalanche Transition Matrix |

### Additional utilities

- **Covariance estimation** — Ledoit–Wolf and OAS shrinkage estimators
- **Riemannian geometry** — SPD manifold operations (Riemannian mean, log/exp maps, geodesic filtering, FGDA)
- **Bandpass filtering** — Butterworth filter bank and data-driven frequency-band selection
- **Feature selection** — Fisher score, PCA, and MIBIF
- **Classifiers** — Fisher LDA and elastic-net multiclass classifier

## Installation

```r
# Install dependencies first (especially when installing from local/source tarballs)
install.packages(c("geigen", "gsignal"))

# Install from GitHub (requires devtools or remotes)
remotes::install_github("Yiming-S/BCIFeatR", dependencies = TRUE)
```

## Quick start

```r
library(BCIFeatR)

# x: list of trial matrices (samples × channels)
# y: factor of class labels (length == number of trials)

# Train — returns extracted features + fitted object for test-time use
fit <- featEx4Train(x_train, y_train, feature = "CSP",
                    params = list(ncomps = 4L))
X_train <- fit$features       # numeric matrix (trials × features)

# Test — deterministic transform using the fitted object
X_test <- featEx4Test(x_test, fit$object, feature = "CSP")
```

### Session-list input

When you have multiple recording sessions, pass lists of trial lists and label
vectors. Each session is fitted/transformed independently using a shared model:

```r
multi_fit <- featEx4Train(
  x = list(session1_trials, session2_trials),
  y = list(session1_labels, session2_labels),
  feature = "logvar", params = list()
)
# Returns a list of per-session results
```

## Data format

- **Trials** (`x`): A list of numeric matrices, each `samples × channels`.
  All trials must have the same number of channels.
- **Labels** (`y`): A factor (or coercible vector) with one entry per trial.

## Dependencies

- **R** ≥ 4.0.0
- **Imports**: `stats`, `geigen`, `gsignal`
- **Suggests**: `infotheo` (for MIBIF feature selection), `testthat` (for tests)

## Testing

```r
devtools::test()
```

## License

MIT — see [LICENSE](LICENSE) for details.
