BT.samples <- function(y,p,to.pred){
  data.1 <- cbind.data.frame(y,p)

  # REPEAT TO PREVENT BT SAMPLES WITH NA'S
  repeat{
    BT.sample.rows <- sample(1:nrow(data.1),replace=T)
    BT.sample <- data.1[BT.sample.rows,]
    loess(y~p,BT.sample) ->loess.BT
    predict(loess.BT,to.pred,type="fitted") ->pred.loess
    if(!any(is.na(pred.loess))){break}
  }
  return(pred.loess)
}
