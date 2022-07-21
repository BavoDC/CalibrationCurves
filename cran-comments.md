## Resubmission
This is a resubmission. In this version I have:
- Changed 'T' and 'F' to 'TRUE' and 'FALSE' in man/val.prob.ci.2.Rd.
- Adjusted print.CalibrationCurve.Rd. Made a custom .Rd file (previous one was inherited from the base print function) and added the \\value field.
- Function val.prob.ci.2 no longer writes messages to the console.
- Graphical settings changed by val.prob.ci.2 are immediately restored when the function is exited.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
