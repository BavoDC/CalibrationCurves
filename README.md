CalibrationCurves: assessing the agreement between observed outcomes and predictions.
====
Package to generate (generalized) calibration curves and related statistics. The function for the logistic/flexible calibration curves are based on the val.prob function from Frank Harrell's [rms](https://cran.r-project.org/package=rms) package.

<p align="left">
  <img src="vignettes/CalibrationCurves.png" width="25%">
</p>

## Installation

### On current R (>= 3.0.0)
* You can install the latest development version from Github using the code below

```
library("devtools")
install_github("BavoDC/CalibrationCurves", dependencies = TRUE, build_vignettes = TRUE, ref = "master")
```

This requires `devtools` >= 1.6.1, and installs the "master" branch. This approach builds the package from source.

## Documentation
The basic functionality of the package is explained and demonstrated in the vignette, which you can access using
```
vignette("CalibrationCurves")
```

or via the [homepage](https://bavodc.github.io/websiteCalibrationCurves/articles/CalibrationCurves.html) of the package. 


## Contact
If you have questions, remarks or suggestions regarding the package, you can contact me at [bavo.campo@kuleuven.be](mailto:bavo.campo@kuleuven.be) (all emails to [bavo.decock@kuleuven.be](mailto:bavo.decock@kuleuven.be) are
forwarded to this one).

## Citation
If you use this package, please cite: <br />
- Barreñada, L., De Cock Campo, B., Wynants, L., Van Calster, B. (2025). Clustered Flexible Calibration Plots for
Binary Outcomes Using Random Effects Modeling. arXiv:2503.08389, available at [https://arxiv.org/abs/2503.08389](https://arxiv.org/abs/2503.08389). <br />
- De Cock Campo, B. (2023). Towards reliable predictive analytics: a generalized calibration framework.
arXiv:2309.08559, available at [https://arxiv.org/abs/2309.08559](https://arxiv.org/abs/2309.08559). <br />
- De Cock, B., Nieboer, D., Van Calster, B., Steyerberg, E.W., Vergouwe, Y. (2023). _The CalibrationCurves package: assessing the agreement between observed outcomes and predictions_. R package version 2.0.3, [doi:10.32614/CRAN.package.CalibrationCurves](https://doi.org/10.32614/CRAN.package.CalibrationCurves), available at [https://cran.r-project.org/package=CalibrationCurves](https://cran.r-project.org/package=CalibrationCurves) <br />
- Van Calster, B., Nieboer, D., Vergouwe, Y., De Cock, B., Pencina, M.J., Steyerberg, E.W. (2016). A calibration hierarchy for risk models was defined: from utopia to empirical data. _Journal of Clinical Epidemiology_, __74__, pp. 167-176 <br />
