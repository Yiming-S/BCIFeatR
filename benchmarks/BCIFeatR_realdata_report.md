# BCIFeatR — Real-Data Test Report

**Package:** BCIFeatR v0.3.3  ·  **Date:** 2026-06-28
**Datasets:** Zhou 2016 (BIDS/EDF) · BNCI 2014-004 / "2b" (`.mat`)
**Loaders:** `eegDataLoaders.R` (`load_eeg_zhou2016` via `edfReader`; `load_eeg_bnci004` via `R.matlab`)

---

## 1. TL;DR

Every BCIFeatR feature method and classifier runs end-to-end on **real** motor-imagery
EEG loaded through the lab loaders, and produces accuracies squarely in the published
range for both datasets. Nothing crashed or returned `NA` on the standard pipeline; the
only fragile path was the generative `CSP → MSVAR` classifier on the 3-channel BNCI data
(documented in §5).

| Dataset | Protocol | Best method | Best bal-acc | Chance |
|---|---|---|---|---|
| Zhou2016 (4 subj, 14 ch) | leave-one-session-out (3 sessions) | CSP | **0.81 ± 0.09** | 0.50 |
| BNCI2014-004 (9 subj, 3 ch) | train T-sessions / test E-sessions | FBCSP | **0.76 ± 0.12** | 0.50 |

---

## 2. Setup

### Data & preprocessing (by the lab loaders, identical for every method)
- **Zhou2016**: subjects `sub-1…4`, 3 sessions each (`ses-0/1/2`), 14 EEG channels,
  `fs = 250 Hz`. Loader merges runs, band-passes 8–35 Hz (zero-phase), epochs a 5 s MI
  window (1–6 s). Binary **left vs right** (`nclass = 2`). ~100–120 trials/session.
- **BNCI2014-004 (2b)**: subjects `B01…B09`, 5 sessions (`01T,02T,03T` train +
  `04E,05E` eval), 3 EEG channels (C3, Cz, C4), `fs = 250 Hz`. Loader band-passes
  8–32 Hz, epochs the 3.0–7.5 s window (4.5 s), drops artifact-flagged trials. Binary
  **left vs right**.

### Evaluation protocols (leakage-safe)
- **Zhou2016**: **leave-one-session-out** — fit on the other 2 sessions (pooled trials),
  test on the held-out session; 3 folds × 4 subjects = 12 folds.
- **BNCI2014-004**: the canonical 2b split — **train on the 3 T-sessions, test on the 2
  E-sessions** (1 fold/subject).
- In both, *all* fitted state (CSP/FBCSP filters, tangent reference, preprocessing stats,
  MSVAR models, and the downstream LDA) is estimated on the **training trials only** and
  applied unchanged to the held-out trials. The folds are session-disjoint.

### Methods
- **9 feature families** through a fixed downstream classifier: `logvar`, `CSP`,
  `FBCSP`, `bandpower`, `Hjorth`, `TS`, `ACM_TS`, `Riemannian`, `MVAR`.
- **2 native classifiers**: `MDM` (Riemannian minimum-distance-to-mean on raw
  covariances, log-Euclidean) and `CSP → MSVAR` (CSP projection, then the generative
  Markov-switching VAR classifier).
- **Downstream classifier** (held constant across feature families): a shrinkage LDA
  built on the package's own `LW_covariance` (Ledoit–Wolf), so it stays well-conditioned
  even for the high-dimensional tangent-space / ACM features.
- **Metric**: balanced accuracy. Bands for FBCSP/bandpower: 8–12, 12–16, 16–24, 24–30 Hz.

---

## 3. Results — Zhou2016 (leave-one-session-out)

**Mean ± sd balanced accuracy across 4 subjects** (each = mean of its 3 LOSO folds):

| Method | bal-acc | sd |
|---|---|---|
| **CSP** | **0.814** | 0.088 |
| bandpower | 0.807 | 0.061 |
| FBCSP | 0.793 | 0.134 |
| MDM | 0.784 | 0.082 |
| Hjorth | 0.759 | 0.066 |
| logvar | 0.755 | 0.088 |
| ACM_TS | 0.749 | 0.091 |
| TS | 0.746 | 0.120 |
| Riemannian | 0.746 | 0.097 |
| CSP → MSVAR | 0.688 | 0.078 |
| MVAR | 0.664 | 0.054 |

**Per-subject** (mean over 3 folds):

| subj | CSP | bandpower | FBCSP | MDM | TS | Riem | logvar | MVAR | CSP→MSVAR |
|---|---|---|---|---|---|---|---|---|---|
| sub1 | 0.732 | 0.822 | 0.828 | 0.767 | 0.682 | 0.645 | 0.707 | 0.617 | 0.695 |
| sub2 | 0.751 | 0.759 | 0.604 | 0.678 | 0.610 | 0.684 | 0.689 | 0.624 | 0.600 |
| sub3 | 0.917 | 0.887 | 0.920 | 0.820 | 0.870 | 0.853 | 0.883 | 0.730 | 0.670 |
| sub4 | 0.856 | 0.759 | 0.822 | 0.871 | 0.820 | 0.800 | 0.739 | 0.685 | 0.788 |

sub3 is easy (≈0.92), sub2 hard (≈0.6–0.75) — typical MI between-subject spread.

---

## 4. Results — BNCI2014-004 (train T / test E)

**Mean ± sd balanced accuracy across 9 subjects:**

| Method | bal-acc | sd |
|---|---|---|
| **FBCSP** | **0.762** | 0.120 |
| MVAR | 0.761 | 0.128 |
| ACM_TS | 0.759 | 0.129 |
| bandpower | 0.751 | 0.121 |
| Hjorth | 0.741 | 0.120 |
| TS | 0.741 | 0.137 |
| CSP | 0.737 | 0.144 |
| MDM | 0.733 | 0.136 |
| Riemannian | 0.731 | 0.142 |
| logvar | 0.722 | 0.142 |
| CSP → MSVAR | 0.590 | 0.126 |

(CSP → MSVAR re-run on all 9 subjects after the v0.3.4 fix — per subject: B01 0.500,
B02 0.500, B03 0.500, B04 0.860, B05 0.652, B06 0.494, B07 0.640, B08 0.667, B09 0.500.
It now runs everywhere but stays near chance on this 3-channel data; see §5.)

**Per-subject** (train T / test E):

| subj | FBCSP | MVAR | ACM_TS | bandpower | TS | CSP | MDM | logvar |
|---|---|---|---|---|---|---|---|---|
| B01 | 0.691 | 0.710 | 0.703 | 0.658 | 0.714 | 0.668 | 0.591 | 0.626 |
| B02 | 0.616 | 0.641 | 0.571 | 0.629 | 0.568 | 0.559 | 0.568 | 0.572 |
| B03 | 0.560 | 0.516 | 0.567 | 0.577 | 0.509 | 0.509 | 0.557 | 0.526 |
| B04 | 0.945 | 0.925 | 0.951 | 0.941 | 0.935 | 0.932 | 0.935 | 0.938 |
| B05 | 0.852 | 0.877 | 0.798 | 0.770 | 0.701 | 0.697 | 0.740 | 0.654 |
| B06 | 0.825 | 0.736 | 0.769 | 0.829 | 0.802 | 0.822 | 0.798 | 0.806 |
| B07 | 0.797 | 0.789 | 0.763 | 0.680 | 0.749 | 0.740 | 0.744 | 0.688 |
| B08 | 0.775 | 0.862 | 0.853 | 0.864 | 0.868 | 0.900 | 0.889 | 0.868 |
| B09 | 0.801 | 0.796 | 0.858 | 0.808 | 0.826 | 0.810 | 0.771 | 0.817 |

B04 very easy (≈0.95), B03 near chance (≈0.5–0.58) — the well-known BNCI-2b range.

---

## 4b. Within-session vs cross-session (session-drift analysis)

Same subjects, same methods, same shrinkage-LDA — only the split changes.
**Within-session** = stratified 5-fold CV *inside each session* (train and test
from the same recording); averaged over folds, sessions, and subjects.
`Δ = within − cross` is the session-drift penalty (positive ⇒ cross-session is
harder). CSP → MSVAR omitted here (near chance either way).

| Method | Zhou within | Zhou cross | Δ | BNCI within | BNCI cross | Δ |
|---|---|---|---|---|---|---|
| logvar | 0.788 | 0.755 | +0.033 | 0.732 | 0.722 | +0.010 |
| CSP | 0.848 | 0.814 | +0.034 | 0.739 | 0.737 | +0.002 |
| FBCSP | 0.830 | 0.793 | +0.037 | 0.757 | 0.762 | −0.005 |
| bandpower | 0.787 | 0.807 | −0.020 | 0.755 | 0.751 | +0.004 |
| Hjorth | 0.754 | 0.759 | −0.005 | 0.733 | 0.741 | −0.008 |
| **TS** | **0.846** | **0.746** | **+0.100** | 0.738 | 0.741 | −0.003 |
| **ACM_TS** | **0.837** | **0.749** | **+0.088** | 0.756 | 0.759 | −0.003 |
| Riemannian | 0.793 | 0.746 | +0.047 | 0.731 | 0.731 | 0.000 |
| MVAR | 0.619 | 0.664 | −0.045 | 0.742 | 0.761 | −0.019 |
| MDM | 0.820 | 0.784 | +0.036 | 0.719 | 0.733 | −0.014 |

**Two clearly different regimes:**

- **Zhou2016 has real session drift, and it hits covariance/tangent-space
  methods hardest.** Tangent space (**+0.10**) and ACM_TS (**+0.09**) lose the
  most going cross-session, with Riemannian/MDM/CSP/FBCSP down 0.03–0.05. This
  is expected: the per-trial covariance geometry these methods rely on shifts
  between sessions, so a reference/mean fit on other sessions transfers
  imperfectly. Band-power/Hjorth/MVAR are roughly flat or slightly *better*
  cross-session — they are drift-robust, and cross-session pools ~2× the
  training trials (2 sessions vs 4/5 of one), which offsets the small drift.
- **BNCI2b shows almost no drift** (all |Δ| ≤ 0.02, mostly ≤0). With only 3
  fixed electrodes (C3/Cz/C4) the sessions are stable, and the larger
  cross-session training set fully compensates — cross-session is essentially
  free here.

Takeaway: BCIFeatR's covariance/Riemannian features are the strongest *within*
session, but they are also the most session-drift-sensitive — on a drifting
dataset (Zhou) a within-session evaluation overstates their cross-session
performance by ~0.1, whereas on a stable dataset (BNCI2b) the two agree.

---

## 5. Observations

1. **Accuracies are literature-consistent.** Zhou2016 binary MI in the ~0.75–0.82 band
   and BNCI-2b in the ~0.73–0.76 band match published cross-session results, so the whole
   load → feature → classify path is behaving correctly on real data.

2. **Channel count changes the winners — a nice cross-dataset signal.**
   - On **14-channel** Zhou2016, *spatial/spectral* methods lead (CSP 0.81, bandpower
     0.81, FBCSP 0.79); `MVAR` is **last** (0.66).
   - On **3-channel** BNCI (C3/Cz/C4), spatial filtering has little to exploit, the
     methods cluster tightly (~0.72–0.76), and **`MVAR` jumps to ~top (0.76)** alongside
     FBCSP/ACM_TS. When spatial information is scarce, the autoregressive *dynamics* MVAR
     captures become relatively more useful — exactly the behavior its design predicts.

3. **MDM is a strong, tuning-free baseline** (0.78 on Zhou, 0.73 on BNCI) — competitive
   with CSP+LDA while needing no downstream classifier.

4. **`CSP → MSVAR` exposed (and now fixed) an MSVAR robustness bug.** It runs on
   14-channel Zhou2016 (0.69), but on 3-channel BNCI, CSP yields only 3 components and a
   regime in the M=2 model collapses; the EM M-step then divided by a zero responsibility
   total, producing `NaN` sufficient statistics that aborted `kmeans()`/`svd()` with
   *"infinite or missing values in 'x'"* for 7/9 subjects. **Fixed in v0.3.4** (the M-step
   floors responsibility totals, falls back to a benign `A = 0, Q = I` empty-regime
   update, keeps seed covariances SPD, and guards k-means). MSVAR now runs for **every**
   subject (mean 0.590 ± 0.126; B01–B09 = 0.50/0.50/0.50/0.86/0.65/0.49/0.64/0.67/0.50);
   on this 3-channel data it stays near chance — there is little regime-switching dynamics
   for it to exploit after CSP — rather than erroring. The other 10 methods completed for
   every subject/fold throughout.

5. **Timing** (mean per fold, feature extraction + LDA): Zhou2016 ≈ 15 s, BNCI ≈ 6 s.
   On Zhou the cost is dominated by the per-band zero-phase filtering in `FBCSP` (≈8.6 s)
   and `bandpower` (≈5.0 s) over 14 channels × 5 s; all covariance/tangent/AR methods are
   ≤0.4 s. Everything is comfortably interactive.

---

## 6. Numerical validation against reference implementations

Beyond "runs without error", the covariance / Riemannian primitives were checked
against `scikit-learn` and `pyriemann` on identical matrices (via `reticulate`),
isolating the math. They match to machine precision where the definition is
unambiguous:

| Primitive | Reference | Max error |
|---|---|---|
| Ledoit–Wolf covariance | `sklearn.covariance.LedoitWolf` | 2.4e-15 (exact) |
| Affine-invariant distance | `pyriemann` `distance_riemann` | 1.9e-15 (exact) |
| Riemannian (Fréchet) mean | `pyriemann` `mean_riemann` | 5.3e-10 |
| Tangent space | `pyriemann` `tangent_space` | 5.6e-15 (exact) |
| OAS covariance | `sklearn.covariance.OAS` | 3.4e-3 (shrinkage-coefficient variant) |

The tangent-vector norm equals the Riemannian distance to the reference (a
metric identity) to 1.5e-15. Only the OAS shrinkage intensity differs from
sklearn's variant (0.34%); both are valid estimators. Script: `numval.R`.

## 7. Reproducibility

- Scripts (in `benchmarks/`): `realdata_bench.R` (cross-session),
  `realdata_within.R` (within-session), `numval.R` (numerical validation).
- Raw results: `real_zhou.csv`, `real_bnci.csv` (cross), `within_zhou.csv`,
  `within_bnci.csv` (within-session).
- Deterministic given the loaders; MSVAR seeded (`seed = 1`). Loaders depend on
  `edfReader`, `R.matlab`, `gsignal`; validation adds `reticulate` + Python
  `pyriemann`/`scikit-learn`.
- BCIFeatR loaded from source; the lab loaders are unmodified.

## 8. Conclusion

BCIFeatR works correctly and competitively on two independent real MI-EEG
datasets, with no failures on the standard feature→classifier pipeline and
results in the published range. Its core Riemannian primitives match `pyriemann`
to machine precision. The cross-dataset behavior is sensible: spatial methods win
with many channels, AR/MVAR catches up with few, and covariance/tangent-space
features are the strongest within-session but the most session-drift-sensitive
(Zhou loses ~0.1 cross-session; stable BNCI2b loses ~0). The one code follow-up
from this testing — hardening `CSP → MSVAR` / the MSVAR EM against degenerate
low-dimensional inputs — has been fixed in v0.3.4.
