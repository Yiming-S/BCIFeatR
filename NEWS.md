# BCIFeatR 0.3.4

Robustness fix from real-data testing (zhou2016, BNCI2014-004).

* The MSVAR EM no longer crashes on degenerate / low-dimensional inputs. When a
  regime collapses (zero responsibility) the M-step previously divided by zero,
  producing `NaN` sufficient statistics that aborted `kmeans()`/`svd()` with
  "infinite or missing values in 'x'" — this hit `CSP -> MSVAR` on the
  3-channel BNCI2014-004 data (7/9 subjects). The M-step now floors the
  responsibility totals, sanitizes non-finite statistics, falls back to a benign
  `A = 0, Q = I` update for empty regimes, keeps cluster-seed covariances SPD,
  and guards the k-means seeding. MSVAR now runs on every subject (it degrades
  gracefully to chance on data with no usable regime dynamics rather than
  erroring). Well-conditioned fits are numerically unchanged.

# BCIFeatR 0.3.3

Documentation and release-hygiene polish (no code-behavior changes).

* Added a `vignette("mi-decoding")` — an end-to-end, simulated motor-imagery
  decoding walkthrough (feature families, classifiers, preprocessing/SimAM,
  MCCA fusion, cross-session input).
* Added runnable `@examples` to many more exported functions (now covers the
  primary entry points, classifiers, covariance/Riemannian utilities, CSP,
  ATM, ACM, and band-pass helpers).
* `inst/CITATION` now derives the version from the package metadata instead of
  hardcoding it.
* snake_case aliases (`feat_ex_train`, `lw_covariance`, `compute_acm`, ...) are
  now direct bindings to their camelCase counterparts, so argument defaults can
  never drift between the two names.
* DESCRIPTION gains `Language: en-US`; README shows the `print()` methods.

# BCIFeatR 0.3.2

Correctness release from a second deep review.

## Bug fixes

* **FBCSP / FBCSSP now band-pass each trial per band at feature-extraction
  time.** Previously the per-band spatial filters (learned on band-passed data)
  were applied to the *broadband* trial, discarding band selectivity; on a weak
  narrow-band signal with strong out-of-band nuisance this collapsed accuracy to
  chance (~0.54 vs ~0.85 once fixed). Note: FBCSP/FBCSSP feature *values* change
  (to the correct ones).
* **`mibif()` no longer silently mis-selects features on balanced designs.**
  With equal rows per session the class vector collapsed to a matrix, which
  `infotheo::discretize()` mis-shaped, zeroing every mutual information and
  returning arbitrary (first) features. Fixed the row-label expansion.
* `pca()` and `fisher()` no longer crash on a constant/zero-variance column
  (e.g. sparse ATM features); the scale is floored as in `logvar_pca`.
* `msvar_train()` validates `M >= 2` and requires at least two trials per class,
  with clear messages instead of cryptic dimension/k-means errors.
* `freqBandSelect()` no longer crashes on single-channel input.

# BCIFeatR 0.3.1

Polish release from an internal improvement review (bug fixes, performance,
consistency, and packaging). No backward-incompatible changes.

## Bug fixes

* `FBCSP` / `FBCSSP` no longer error when a `channels` subset is supplied; the
  subset is now applied consistently at train and test time.
* `mdm_predict()` measures distance with the same geometry used to build the
  class means (Frobenius for `euclid`, log-Euclidean for `logeuclid`,
  affine-invariant for `riemann`), so the nearest-mean rule is coherent for
  every metric (previously the distance was always affine-invariant).
* The `Riemannian` tie-breaking jitter is now batch-invariant: a trial yields
  the same features whether transformed alone or within a larger batch.
* `logvar_pca` no longer crashes on a constant/dead channel.
* `mdm_train(cov_type = "cov")` errors clearly on rank-deficient (p >= n) trials
  instead of producing silent garbage; use `oas`/`lw` shrinkage instead.
* SimAM's `robust` variant is auto-scaled to the data (consistent MAD), so its
  default no longer collapses the signal; `lambda` is a scale-invariant relative
  bandwidth.
* The vendored MSVAR Viterbi decode now runs in log space (no underflow on long
  trials); an uninitialized `LLflag` in the EM loop was fixed.

## New

* S3 `predict()` and `print()` methods for `mdm`, `msvar`, and `mcca` objects
  (e.g. `predict(model, x)`, `predict(model, x, type = "distance"/"loglik")`).
* Runnable `@examples` on the primary entry points; an `inst/CITATION`; a
  `NEWS.md`; and a GitHub Actions R-CMD-check workflow.

## Performance

* ATM time-binning vectorized (~68x faster).
* CSP class-mean covariances cached across one-vs-one pairs (~`nclass - 1`x,
  multiplied per band in FBCSP).
* `logvar_transform()` uses the vectorized column-variance helper (~4x faster).
* `mdm_predict()` precomputes per-class quantities once instead of per trial.

# BCIFeatR 0.3.0

## New features

* New dispatcher feature methods in `featEx4Train()` / `featEx4Test()`:
  `bandpower` (per-channel log/relative band power), `Hjorth`
  (activity/mobility/complexity), `MVAR` (vector-autoregression coefficients),
  and `MSVAR` (per-class Markov-switching VAR log-likelihood).
* New classifiers `mdm_train()` / `mdm_predict()` (Riemannian Minimum Distance
  to Mean) and `msvar_train()` / `msvar_predict()` (generative Markov-switching
  VAR).
* New multi-view fusion `mcca_train()` / `mcca_transform()` (regularized SUMCOR
  multiset CCA) and parameter-free `simAM()` attention (also via `params$simam`).
* New `fit_var()` and `splitTimeRange()`; `freqBandSelect()` is now functional.

## Other

* Removed the unused `geigen` dependency (hard deps are `stats` + `gsignal`).
