## Update
This is an update. In this version:
- R/val.prob.ci.2.R (lines 207 - 246): the previous code did not return the correct output when length(unique(p)) == 1. This issue has been addressed. 
- man/val.prob.ci.2.Rd (details and references): we explain what the function returns in case of an uninformative model and added the reference to Edlinger et al. (2021).
- R/valProbggPlot.R: this is a new function, which is the 'ggplot' version of the val.prob.ci.2 function.
- R/print.ggplotCalibrationCurve.R: added a print function for the output of the valProbggPlot function.
- man/valProbggPlot.Rd: help-file for the valProbggPlot function.
- man/print.ggplotCalibrationCurve.Rd: help-file for the print.ggplotCalibrationCurve function.
- NAMESPACE: has been adjusted accordingly.
