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
| `bandpower` | Per-channel (log or relative) band power over a frequency-band bank |
| `Hjorth` | Per-channel Hjorth parameters (activity, mobility, complexity) |
| `MVAR` | Per-trial vector-autoregression coefficients (+ residual log-variance) |
| `MSVAR` | Per-class Markov-switching VAR log-likelihood (generative, supervised) |

### Additional utilities

- **Covariance estimation** — Ledoit–Wolf and OAS shrinkage estimators
- **Riemannian geometry** — SPD manifold operations (Riemannian mean, log/exp maps, geodesic filtering, FGDA)
- **Bandpass filtering** — Butterworth filter bank and data-driven frequency-band selection
- **Feature selection** — Fisher score, PCA, and MIBIF
- **Multi-view fusion** — regularized SUMCOR multiset CCA (`mcca_train` /
  `mcca_transform`) for aligning multiple feature views
- **Attention** — parameter-free SimAM re-weighting (`simAM`), also available as
  the `params$simam` switch in the train/test dispatcher
- **Classifiers** — Fisher LDA, elastic-net multiclass classifier, a
  Riemannian Minimum Distance to Mean (MDM) classifier (`mdm_train` /
  `mdm_predict`), and a per-class Markov-switching VAR generative classifier
  (`msvar_train` / `msvar_predict`)

## Installation

```r
# Install dependencies first (especially when installing from local/source tarballs)
install.packages("gsignal")

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

### Classifiers

```r
# Riemannian Minimum Distance to Mean (operates on raw trials)
mdm  <- mdm_train(x_train, y_train, metric = "logeuclid")
pred <- mdm_predict(mdm, x_test)$y_hat

# Markov-switching VAR generative classifier (strongest on low-dimensional,
# e.g. CSP-projected, trials; argmax of the per-class trial log-likelihood)
msv  <- msvar_train(x_train, y_train, M = 2L, p = 1L, seed = 1L)
out  <- msvar_predict(msv, x_test)   # out$y_hat (labels), out$loglik (scores)
```

### Multi-view fusion (MCCA)

```r
# Align several feature views (rows = the same trials) into shared components
V1 <- featEx4Train(x_train, y_train, "logvar", list())$features
V2 <- featEx4Train(x_train, y_train, "bandpower",
                   list(fs = 250, frequency_bands = list(c(8, 13), c(13, 30))))$features
mc          <- mcca_train(list(V1, V2), ncomp = 4L)
fused_train <- mcca_transform(mc, list(V1, V2))$consensus
```

### Preprocessing & attention

```r
# Channel preprocessing (none/center/scale/whiten) and optional parameter-free
# SimAM attention re-weighting apply before any feature, via `params`:
fit <- featEx4Train(x_train, y_train, "TS",
                    params = list(cov_type   = "oas",
                                  preprocess = "whiten",
                                  simam      = "vanilla"))
```

## Data format

- **Trials** (`x`): A list of numeric matrices, each `samples × channels`.
  All trials must have the same number of channels.
- **Labels** (`y`): A factor (or coercible vector) with one entry per trial.

## Dependencies

- **R** ≥ 4.0.0
- **Imports**: `stats`, `gsignal`
- **Suggests**: `infotheo` (for MIBIF feature selection), `testthat` (for tests)

## Testing

```r
devtools::test()
```

## References

The Markov-switching VAR engine behind the `MSVAR` feature and the
`msvar_train()` / `msvar_predict()` classifier is adapted from the Degras Lab
`msvar` implementation:

> Degras, D., Ting, C.-M., & Ombao, H. (2022). Markov-switching state-space
> models with applications to neuroimaging. *Computational Statistics & Data
> Analysis*.

## License

MIT — see [LICENSE](LICENSE) for details.
