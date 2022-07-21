BT.samples <- function(y, p, to.pred){
  Df = cbind.data.frame(y, p)

  # REPEAT TO PREVENT BT SAMPLES WITH NA'S
  repeat {
    BT.sample  = Df[sample(1:nrow(Df), replace = T), ]
    loess.BT   = loess(y ~ p, BT.sample)
    pred.loess = predict(loess.BT, to.pred, type = "fitted")
    if (!any(is.na(pred.loess)))
      break
  }
  return(pred.loess)
}
