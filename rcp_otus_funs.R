colMedian <- function(x) {
  meds <- numeric(ncol(x))
  for (i in 1:ncol(x)){
    meds[i] <- median(x[,i])
  }
  return(meds)
}

colMin <- function(x) {
  mins <- numeric(ncol(x))
  for (i in 1:ncol(x)){
    mins[i] <- min(x[,i])
  }
  return(mins)
}

ParPlotIndiv.regimix = function(fm, covar.data, newdata, variable, RCPs) {
  var.range = seq(min(covar.data[,variable]),max(covar.data[,variable]),length.out=1000)
  newdata[,grep(variable, colnames(newdata))] = predict(poly(covar.data[,variable],2),newdata=var.range)
  # do prediction
  preds = predict.regimix(fm, newdata=newdata, nboot=0)
  # plot it
  for (i in RCPs) {
    matplot(var.range, preds[,i], type='l', ylab=paste0("RCP",i), 
            main=paste0("Partial effect of ",variable), xlab=variable)
  }
}

ParPlotMany.regimix = function(fm, covar.data, newdata, variable, RCPs) {
  var.range = seq(min(covar.data[,variable]),max(covar.data[,variable]),length.out=1000)
  newdata[,grep(variable, colnames(newdata))] = predict(poly(covar.data[,variable],2),newdata=var.range)
  # do prediction
  preds = predict.regimix(fm, newdata=newdata, nboot=0)
  # plot it
  matplot(var.range, preds[,RCPs], type='l', lty=rep(1:3,length.out=length(RCPs)), col=RCPs, lwd=3,  
          ylab=paste0("Probability"), main=paste0("Partial effect of ",variable), xlab=variable)
  legend("topleft", legend=paste0("RCP",RCPs), lty=rep(1:3,length.out=length(RCPs)),
         col=RCPs, lwd=3, horiz=FALSE)
}