## Submission summary

BCIFeatR is a feature-extraction toolkit for EEG-based brain-computer interfaces
(CSP/FBCSP/FBCSSP, Riemannian tangent space, band-power, Hjorth, AR/MSVAR
features, MDM and MSVAR classifiers, MCCA fusion, SimAM attention).

## Test environments

* local: macOS, R 4.4.2
* GitHub Actions (r-lib/actions): ubuntu-latest, macOS-latest, windows-latest
  (R release), via .github/workflows/R-CMD-check.yaml

## R CMD check results

0 errors | 0 warnings | 1 note

* The only local NOTE ("checking for future file timestamps ... unable to
  verify current time") is environmental (no network access to the time server
  during the check) and does not occur on CRAN's build machines.
* On first submission the standard "New submission" NOTE is expected.

## Downstream dependencies

There are currently no reverse dependencies.

## Notes

* Runtime dependencies are `stats` and `gsignal`; `infotheo`, `knitr`,
  `rmarkdown`, and `testthat` are Suggests (all uses are guarded).
* The Markov-switching VAR engine is adapted from the authors' `msvar`
  implementation; the method reference (Degras, Ting & Ombao, 2022) is cited in
  `inst/CITATION` and the relevant help pages.
