auc.logit.ci <- function(crit,pred,auc,conf.level=0.95){

  tmp <- cbind.data.frame(crit=crit,pred=pred)

  ntot <- length(crit)
  tmp$rank.pred <- rank(tmp$pred)
  tmp <- tmp[with(tmp,order(crit,rank.pred)),]

  # healthy
  nondis <- tmp[which(tmp$crit==0),]
  m <- nrow(nondis)

  # disease
  dis <- tmp[which(tmp$crit==1),]
  n <- nrow(dis)


  s210 <- (1/((m-1)*n^2))*((m*mean((nondis$rank.pred-1:m)^2))-(m*(mean(nondis$rank.pred)-0.5*(m+1))^2))

  s201 <- (1/((n-1)*m^2))*((n*mean((dis$rank.pred-1:n)^2))-(n*(mean(dis$rank.pred)-0.5*(n+1))^2))

  s <- sqrt((m*s201+n*s210)/ntot)
  vard <- (ntot*s^2)/(m*n)
  lcl <- plogis(qlogis(auc)+qnorm((1-conf.level)/2)*(sqrt(vard))/(auc*(1-auc)))
  ucl <- plogis(qlogis(auc)-qnorm((1-conf.level)/2)*(sqrt(vard))/(auc*(1-auc)))

  return(c(auc,lcl,ucl))
}






