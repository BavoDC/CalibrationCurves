CalibrationCurves: Validating predicted probabilities against binary events.
====
Package to generate logistic and flexible calibration curves and related statistics. Based on the val.prob function from Frank Harrell's [rms](https://cran.r-project.org/web/packages/rms/index.html) package.

## Installation

### On current R (>= 3.0.0)
* Development version from Github:

```
library("devtools"); install_github("BavoDC/CalibrationCurves", dependencies = TRUE)
```

(This requires `devtools` >= 1.6.1, and installs the "master" (development) branch.)
This approach builds the package from source, i.e. `make` and compilers must be installed on your system -- see the R FAQ for your operating system; you may also need to install dependencies manually.

## Contact
If you have questions, remarks or suggestions regarding the package, you can contact me at [bavo.decock@kuleuven.be](mailto:bavo.decock@kuleuven.be).

## Citation
If you use this package, please cite:

- De Cock, B., Nieboer, D., Van Calster, B., Steyerberg, E.W., Vergouwe, Y. (2016). _CalibrationCurves: Calibration performance_. R package version 0.1.2. [https://github.com/BavoDC/CalibrationCurves](https://github.com/BavoDC/CalibrationCurves)
- Van Calster, B., Nieboer, D., Vergouwe, Y., De Cock, B., Pencina, M.J., Steyerberg, E.W. (2016). A calibration hierarchy for risk models was defined: from utopia to empirical data. _Journal of Clinical Epidemiology_, __74__, pp. 167-176
