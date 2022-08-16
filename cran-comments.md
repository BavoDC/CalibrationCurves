## Resubmission
This is a resubmission. In this version:
- man/dot-rcspline.plot.Rd: removed \var in \eqn (resulted in a note in the pretest on Debian).
- R/ci.auc.R (line 14): replaced cat() by warning().
- R/rcspline.plot.noprint.R: 
  * line 203: replaced cat() by message().
  * Graphical settings changed by rcspline.plot.noprint.R are immediately restored when the function is exited.
- R/val.prob.ci.2.R (line 302): removed the cat() statement.

In the previous resubmission, I:
- Changed 'T' and 'F' to 'TRUE' and 'FALSE' in man/val.prob.ci.2.Rd.
- Adjusted print.CalibrationCurve.Rd. Made a custom .Rd file (previous one was inherited from the base print function) and added the \value field.
- Function val.prob.ci.2 no longer writes messages to the console.
- Graphical settings changed by val.prob.ci.2 are immediately restored when the function is exited.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
