# CRAN Submission Checklist - CalibrationCurves 3.1.0

## Pre-Submission Steps

### Quality Checks
- [ ] Run `R CMD check --as-cran` and ensure no ERRORs or WARNINGs
- [ ] Test package on multiple R versions (if possible)
- [ ] Verify all examples run without errors
- [ ] Check for any NOTEs in the check output

### Testing
- [ ] Test `valProbCluster()` with various `cl.level` values:
  - [ ] cl.level = 0.80
  - [ ] cl.level = 0.90
  - [ ] cl.level = 0.95
  - [ ] cl.level = 0.99
- [ ] Verify plot labels display correct percentages for each cl.level
- [ ] Test `valProbSurvival()` with timeHorizon near max follow-up
- [ ] Verify no "incorrect number of dimensions" errors occur

### Documentation
- [ ] Update README.md if needed (versions, new features)
- [ ] Verify all DESCRIPTION fields are correct
- [ ] Check that NAMESPACE is properly updated
- [ ] Ensure NEWS.md is properly formatted

## Files Ready for CRAN

### Core Changes
- `DESCRIPTION` - Updated with new version (3.1.0) and revised dependencies
- `NAMESPACE` - Updated imports (removed timeROC, added pec)
- `NEWS.md` - Complete changelog for version 3.1.0
- `CRAN_SUBMISSION_SUMMARY.md` - Detailed technical summary (for reference)

### Modified R Files
- `R/timeROC_archived.R` (NEW) - Integrated timeROC functions
- `R/valProbCluster.R` - Fixed cl.level parameter passing
- `R/valProbSurvival.R` - Fixed timeHorizon issue
- `R/CGC.R` - Dynamic CI/PI labels, meta-analysis updates
- `R/MAC2.R` - Dynamic CI/PI labels, meta-analysis updates
- `R/MIXC.R` - Dynamic CI/PI labels
- `R/helperFunctions.R` - New ci_pi_labels() function

## Submission Steps

1. **Update Version**
   - [ ] Increment version in DESCRIPTION to 3.1.0 (if not already done)
   - [ ] Ensure Date field is current

2. **Run R CMD check**
   ```bash
   R CMD check --as-cran CalibrationCurves_3.1.0.tar.gz
   ```

3. **Create Submission Package**
   - [ ] Build source package: `R CMD build .`

4. **Submit to CRAN**
   - [ ] Go to https://cran.r-project.org/submit.html
   - [ ] Upload CalibrationCurves_3.1.0.tar.gz
   - [ ] Write submission comment: Include bug fix details for Issues #22 and #24, timeROC integration

5. **Example Submission Comment**
   ```
   Dear CRAN Team,

   Please find attached CalibrationCurves 3.1.0 for submission to CRAN.

   This release includes the following major changes:

   1. **Integration of Archived timeROC Package**: The timeROC package (v0.4) is being archived on CRAN. This release incorporates the necessary functions from timeROC directly into CalibrationCurves to ensure continuity of functionality. (Fixes potential future dependency issues)

   2. **Bug Fix - Issue #22**: The `cl.level` parameter in `valProbCluster()` was hardcoded to 0.95, ignoring user-specified confidence levels. This is now fixed throughout the entire pipeline, including meta-analysis and visualization.

   3. **Bug Fix - Issue #24**: The `valProbSurvival()` function crashed when `timeHorizon` was set near the maximum follow-up time. This is now fixed by calculating Uno's AUC at the user-specified timeHorizon rather than the extreme tail.

   All changes are fully backward compatible with default parameter values preserved.

   Thank you,
   [Author Name]
   ```

## Post-Submission

- [ ] Monitor CRAN email for review feedback
- [ ] Be prepared to address any WARNINGs or NOTEs
- [ ] Update version to next development version after acceptance

## Additional Notes

- **Backward Compatibility**: All changes maintain backward compatibility. Default `cl.level = 0.95` preserves original behavior.
- **Test Data**: Use the provided survival test data in the package for testing Issue #24 fixes: `data(testDataSurvival)`
- **Dependencies**: Verify that `pec` package is accessible and properly specified in DESCRIPTION

---

**Package**: CalibrationCurves  
**New Version**: 3.1.0  
**Previous Version**: 3.0.0  
**Date Prepared**: March 24, 2026
