# CalibrationCurves 3.1.0

## Major changes

* **Archived timeROC dependency**: The `timeROC` package (v0.4) is being archived
  on CRAN. Its required functions have been integrated directly into
  CalibrationCurves (`R/timeROC_archived.R`) with proper attribution. `timeROC`
  has been removed from `Imports`; `pec` has been added instead (required for
  `ipcw()`).

## Bug fixes

* **`valProbCluster()`: `cl.level` parameter was ignored (Issue #22)**: The
  confidence level was hardcoded to 0.95 throughout the clustering pipeline.
  `cl.level` is now correctly passed to `CGC()`, `MAC2()`, and `MIXC()` and
  propagated to `metaprop()`, `metagen()`, and `rma.mv()` calls. Hardcoded
  "95%" plot labels have been replaced with dynamic labels via the new
  `ci_pi_labels()` helper function.

* **`valProbSurvival()`: crash near max follow-up (Issue #24)**: Uno's
  time-dependent AUC was evaluated at `max(fit$y) - 0.01` instead of the
  user-specified `timeHorizon`, causing "incorrect number of dimensions" errors
  when the risk set was depleted. Fixed to use `times = timeHorizon`.

## Enhancements

* **New `"default"` approach in `valProbCluster()`**: Combines MAC2 (splines)
  for the overall calibration curve, confidence intervals, and prediction
  intervals, with MIXC for cluster-specific curves. This is now the default
  when `approach` is not specified. The returned object contains both the MAC2
  overall results (`results$overall`) and the MIXC cluster results
  (`results$clusters`).

* Unified the legend title to "Heterogeneity" across all `valProbCluster()`
  approaches.

* Updated plot font to sans-serif and increased base size to 11 for improved
  readability across MIXC, MAC2, and CGC approaches.

## Minor changes

* Added input validation for `cl.level` in `valProbCluster()`.
* New helper function `ci_pi_labels()` for formatting confidence/prediction
  interval labels dynamically.
* Added documentation of the `"default"` approach to the package vignette.

# CalibrationCurves 3.0.0

* Added `valProbCluster()` for calibration of clustered data with three
  approaches: MIXC, MAC2, and CGC.
* Added `valProbSurvival()` for calibration of survival/time-to-event data.
* Added `genCalCurve()` for generalized calibration across the exponential
  family.
