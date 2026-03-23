#' @keywords internal
"_PACKAGE"

## Base R imports ----
#' @importFrom grDevices rgb
#' @importFrom graphics abline arrows clip legend lines locator par plot points
#'   polygon segments text persp
#' @importFrom methods existsFunction
#' @importFrom stats approx confint loess na.omit pchisq plogis predict
#'   quantile supsmu pnorm qnorm qt binomial coef lm glm fitted glm.control
#'   vcov setNames update loess.control model.matrix optimize qlogis sd var
#' @importFrom utils getFromNamespace sessionInfo globalVariables

## External package imports ----
#' @import rms
#' @import ggplot2
#' @importFrom Hmisc cut2 groupn label rcspline.eval
#' @importFrom survival basehaz coxph coxph.control concordance survfit
#' @importFrom timeROC timeROC
#' @importFrom riskRegression predictRisk Score
#' @importFrom meta metaprop metagen
#' @importFrom metafor escalc rma.mv
#' @importFrom zoo na.approx
#' @importFrom lme4 glmer VarCorr
#' @importFrom merTools predictInterval
#' @importFrom dplyr filter mutate group_by ungroup ntile summarise n
#' @importFrom magrittr %>%
NULL
