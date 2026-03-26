# CalibrationCurves v3.1.0  CRAN Submission Summary

## Changes

### 1. Archived timeROC dependency
The `timeROC` package (v0.4) is being archived on CRAN. Its required functions have been copied verbatim into `R/timeROC_archived.R` with proper attribution. `timeROC` has been removed from `Imports`; `pec` was added instead (required for `ipcw()`).

### 2. Bug fix: `valProbCluster` ignoring `cl.level` (Issue #22)
The `cl.level` parameter was hardcoded to 0.95 throughout the clustering pipeline. Fixed by:
- Passing `cl.level` correctly to `CGC`, `MAC2`, and `MIXC`.
- Propagating it to `metaprop()`, `metagen()`, and `rma.mv()` calls.
- Replacing hardcoded "95%" plot labels with a new `ci_pi_labels()` helper that formats labels dynamically.
- Adding input validation for `cl.level`.

### 3. Bug fix: `valProbSurvival` crash near max follow-up (Issue #24, reported by lisch7)
Uno's time-dependent AUC was evaluated at `max(fit$y) - 0.01` instead of `timeHorizon`, causing dimension errors when the risk set was depleted. Fixed to use `times = timeHorizon`.

---

**Date**: March 2026 | **Status**: Ready for CRAN submission
