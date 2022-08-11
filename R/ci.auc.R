ci.auc <- function(crit,pred,conf.level=0.95,method="pepe"){
  tmp <- cbind.data.frame(crit = crit, pred = pred)

  # healthy
  nondis <- tmp[which(tmp$crit == 0), ]

  # disease
  dis <- tmp[which(tmp$crit == 1), ]

  # ci auc
  if (!grepl("bootstrap", method)) {
    result <- auc.nonpara.mw(dis$pred, nondis$pred, conf.level, method)
  } else if (grepl("bootstrap", method)) {
    warning(
      "Bootstrap-based methods are not supported by this package. Method will be set to 'pepe'. \n\n"
    , immediate. = TRUE)
    result <- auc.nonpara.mw(dis$pred, nondis$pred, conf.level, method = "pepe")
  }

  return(result)
}






