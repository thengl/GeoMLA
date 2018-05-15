## Functions to use RF to make spatial predictions
## tom.hengl@gmail.com

## General Settings
## Legend for plots:
leg = c("#0000ff", "#0028d7", "#0050af", "#007986", "#00a15e", "#00ca35", "#00f20d", "#1aff00", "#43ff00", "#6bff00", "#94ff00", "#bcff00", "#e5ff00", "#fff200", "#ffca00", "#ffa100", "#ff7900", "#ff5000", "#ff2800", "#ff0000")
axis.ls = list(at=c(4.8,5.7,6.5,7.6), labels=round(expm1(c(4.8,5.7,6.5,7.6))))
## 1 s.d. quantiles
quantiles = c((1-.682)/2, 0.5, 1-(1-.682)/2)

## convert standard deviation to 80% prob prediction range
sd.range = function(m, q, t=1){ 
  l = (m-t*q)
  u = (m+t*q)
  out = u-ifelse(l<0, 0, l)
  return(out)
}

source_https <- function(url, ...) {
  # load package
  require(RCurl)
  # download:
  cat(getURL(url, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl")), file = basename(url))
  source(basename(url))
  unlink(basename(url))
}

## derive Scaled Shannon Entropy (100 is a maximum error; 0 is perfect prediction)
entropy_index <- function(x){
    require(entropy)
    v <- unlist(plyr::alply(x, 1, .fun=function(i){entropy.empirical(unlist(i))})) 
    SSI <- round(v/entropy.empirical(rep(1/ncol(x),ncol(x)))*100)
    return(SSI)
}

pfun.line <- function(x,y, ...){
  panel.xyplot(x,y, ...)  
  panel.abline(0,1,lty=2,lw=1,col="black")
}

## standard variogram plot:
plot_vgm <- function(formulaString, rmatrix, predictionDomain, r1, r2, main){
  require(GSIF)
  require(gstat)
  v <- GSIF::fit.vgmModel(formulaString, rmatrix, predictionDomain, dimensions="2D")
  plot(x=v$svgm$dist, y=v$svgm$gamma, pch="+", col = "grey18", xlab='distance', cex=1.1, ylab='gamma', ylim = c(0, max(v$svgm$gamma)), main=main)
  vline <- gstat::variogramLine(v$vgm, maxdist=max(v$svgm$dist), n=length(v$svgm$dist))
  lines(x=vline$dist, y=vline$gamma, col="darkgrey", lwd=2)
  ## validation residuals
  v.r1 <- gstat::variogram(as.formula(paste(r1, "~1")), v$observations[r1]) 
  points(x=v.r1$dist, y=v.r1$gamma, pch="+", col = "red", cex=1.4)
  v.r2 <- gstat::variogram(as.formula(paste(r2, "~1")), v$observations[r2])
  points(x=v.r2$dist, y=v.r2$gamma, pch="+", col = "blue", cex=1.4)
}

cv_numeric <- function(varn, points, covs, nfold=5, idcol, method="ranger", cpus=1, Nsub=1e4, OK=FALSE, spcT=TRUE, Log=FALSE, LLO=TRUE, pars.ranger, predDist=NULL){
  points = points[!is.na(points@data[,varn]),]
  if(missing(idcol)) { 
    points$SOURCEID = row.names(points@data)
    idcol = "SOURCEID"
  }
  message(paste("Running ", nfold, "-fold cross validation with model re-fitting method ", method," ...", sep=""))
  if(nfold > nrow(points@data)){ 
    stop("'nfold' argument must not exceed total number of points") 
  }
  if(sum(duplicated(points@data[,idcol]))>0.5*nrow(points@data)){
    if(LLO==TRUE){
      ## TH: Leave whole locations out
      ul <- paste(unique(points@data[,idcol]))
      sel.ul <- dismo::kfold(ul, k=nfold)
      sel <- lapply(1:nfold, function(o){ data.frame(row.names=which(points@data[,idcol] %in% ul[sel.ul==o]), x=rep(o, length(which(points@data[,idcol] %in% ul[sel.ul==o])))) })
      sel <- do.call(rbind, sel)
      sel <- sel[order(as.numeric(row.names(sel))),]
      message(paste0("Subsetting observations by unique location"))
    } else {
      sel <- dismo::kfold(points@data, k=nfold, by=points@data[,idcol])
      message(paste0("Subsetting observations by '", idcol, "'"))
    }
  } else {
    sel <- dismo::kfold(points@data, k=nfold)
    message(paste0("Simple subsetting of observations using kfolds"))
  }
  if(missing(cpus)){
    cpus <- parallel::detectCores() 
  }
  if(cpus>1){
    require(snowfall)
    snowfall::sfInit(parallel=TRUE, cpus=ifelse(nfold>cpus, cpus, nfold))
    snowfall::sfExport("predict_parallelP","idcol","points","covs","sel","varn","method","Nsub","OK","spcT","pars.ranger")
    snowfall::sfLibrary(package="plyr", character.only=TRUE)
    snowfall::sfLibrary(package="GSIF", character.only=TRUE)
    if(method=="ranger"){
      snowfall::sfLibrary(package="ranger", character.only=TRUE)
      out <- snowfall::sfLapply(1:nfold, function(j){predict_parallelP(j, sel=sel, idcol=idcol, varn=varn, points=points, covs=covs, method=method, cpus=cpus, Nsub=Nsub, OK=OK, spcT=spcT, pars.ranger=pars.ranger, predDist=predDist)})
    } else {
      snowfall::sfLibrary(package="geoR", character.only=TRUE)
      out <- snowfall::sfLapply(1:nfold, function(j){predict_parallelP(j, sel=sel, idcol=idcol, varn=varn, points=points, covs=covs, method=method, cpus=cpus, Nsub=Nsub, OK=OK, spcT=spcT, pars.ranger=pars.ranger, predDist=predDist)})
    }
    snowfall::sfStop()
  } else {
    out <- lapply(1:nfold, function(j){predict_parallelP(j, sel=sel, idcol=idcol, varn=varn, points=points, covs=covs, method=method, cpus=cpus, Nsub=Nsub, OK=OK, spcT=spcT, pars.ranger=pars.ranger, predDist=predDist)})
  }
  ## calculate mean accuracy:
  out <- plyr::rbind.fill(out)
  out$z_score = (out$Observed - out$Predicted)/out$sdPE
  ME = mean(out$Observed - out$Predicted, na.rm=TRUE)
  ##>> MN: MAE - MAD ----
  # maybe remove, MAE is not in the paper. Or better use MADE median absolute deviation error?
  MAE = mean(abs(out$Observed - out$Predicted), na.rm=TRUE)
  RMSE = sqrt(mean((out$Observed - out$Predicted)^2, na.rm=TRUE))
  MZS = mean(out$z_score, na.rm=TRUE)
  ZSV = sd(out$z_score, na.rm=TRUE)
  ## Errors of errors:
  MAE.SE = mean(abs(out$Observed - out$Predicted) - out$sdPE, na.rm=TRUE)
  ## https://en.wikipedia.org/wiki/Coefficient_of_determination
  #R.squared = 1-sum(( (out$Observed - out$Predicted) - mean(out$Observed - out$Predicted) )^2)/(var(out$Observed, na.rm=TRUE)*(sum(!is.na(out$Observed))-1))
  R.squared = 1 - (t(out$Observed - out$Predicted) %*% (out$Observed - out$Predicted)) / (t(out$Observed - mean(out$Observed)) %*% (out$Observed - mean(out$Observed)))
  ccc = DescTools::CCC(out$Observed, out$Predicted, ci = "z-transform", conf.level = 0.95, na.rm=TRUE)$rho.c
  if(Log==TRUE){
    ## If the variable is log-normal then logR.squared is probably more correct
    ##>> MN: logRMSE ----
    # you mention the unbiased backtransform in the paper (Eq. 12 and 13), 
    # predictions on logscale are then not achieved by log1p(out$Predicted)
    # Moreover, I do not know if logRMSE, logR2 should be preferred, maybe rather report both
    # or just backtransfored (In the end most people would prefer to know model 
    # performance on original scale not on rather abstract logscale). 
    logRMSE = sqrt(mean((log1p(out$Observed) - log1p(out$Predicted))^2, na.rm=TRUE))
    #logR.squared = 1-sum((log1p(out$Observed) - log1p(out$Predicted))^2, na.rm=TRUE)/(var(log1p(out$Observed), na.rm=TRUE)*sum(!is.na(out$Observed)))
    logR.squared = 1 - (t(log1p(out$Observed) - log1p(out$Predicted)) %*% (log1p(out$Observed) - log1p(out$Predicted))) / (t(log1p(out$Observed) - mean(log1p(out$Observed)) %*% (log1p(out$Observed) - mean(log1p(out$Observed)))))
    cv.r <- list(out, data.frame(ME=ME, MAE=MAE, RMSE=RMSE, MAE.SE=MAE.SE, MZS=MZS, ZSV=ZSV, R.squared=R.squared, logRMSE=logRMSE, logR.squared=logR.squared, CCC_est=ccc$est, CCC_lwr=ccc$lwr.ci, CCC_upr=ccc$upr.ci)) 
  } else {
    cv.r <- list(out, data.frame(ME=ME, MAE=MAE, RMSE=RMSE, MAE.SE=MAE.SE, MZS=MZS, ZSV=ZSV, R.squared=R.squared, CCC_est=ccc$est, CCC_lwr=ccc$lwr.ci, CCC_upr=ccc$upr.ci))
  }
  message("DONE")
  names(cv.r) <- c("CV_residuals", "Summary")
  return(cv.r)
}

predict_parallelP <- function(j, sel, idcol, varn, points, covs, method, cpus, Nsub, OK=FALSE, spcT=TRUE, pars.ranger, predDist=NULL){ 
  message(paste0("Running ", j, " iteration..."))
  # For method == "geoR" spcT has no meaning, the same for method == "ranger" & OK = T
  s.train <- points[!sel==j,]
  s.test <- points[sel==j,]
  if(!Nsub>nrow(s.train)){ 
    s.train = s.train[sample.int(nrow(s.train), Nsub),]
  }
  if(method=="ranger"){
    require(ranger)
    require(GSIF)
    dist0 <- GSIF::buffer.dist(s.train[varn], covs, as.factor(1:nrow(s.train)))
    dn0 <- paste(names(dist0), collapse="+")
    ovT <- over(s.train[varn], dist0)
    ovV <- over(s.test[varn], dist0)
  }
  ## OK indicates whether geographical distances should be used 
  if(OK==FALSE){
    if(spcT==TRUE){
      covs = GSIF::spc(covs, as.formula(paste0("~ ", paste(names(covs), collapse = "+"))))@predicted
      fm0 <- as.formula(paste(varn, " ~ ", dn0, " + ", paste(names(covs)[-length(covs)], collapse = "+")))
    } else {
      fm0 <- as.formula(paste(varn, " ~ ", dn0, " + ", paste(names(covs), collapse = "+")))
    }
  }
  if(method=="ranger"){
    if(OK==FALSE){
      rmatrix = do.call(cbind, list(s.train@data, ovT, over(s.train[varn], covs)))
      rmatrix.test = do.call(cbind, list(s.test@data, ovV, over(s.test[varn], covs)))
    } else {
      fm0 <- as.formula(paste(varn, " ~ ", dn0))
      rmatrix = do.call(cbind, list(s.train@data, ovT))
      rmatrix.test = do.call(cbind, list(s.test@data, ovV))
    }
    #gm <- ranger(fm0, rmatrix[complete.cases(rmatrix),], keep.inbag = TRUE)
    if(missing(pars.ranger)){
      gm <- ranger(fm0, rmatrix[complete.cases(rmatrix),])
    } else {
      pars.ranger$mtry = ifelse(pars.ranger$mtry >= length(all.vars(fm0)), length(all.vars(fm0))-1, pars.ranger$mtry)
      ##  mtry can not be larger than number of variables in data
      rmatrix0 = rmatrix[complete.cases(rmatrix),]
      gm <- ranger(fm0, rmatrix0, mtry=pars.ranger$mtry, min.node.size=pars.ranger$min.node.size, num.trees = pars.ranger$num.trees, sample.fraction=pars.ranger$sample.fraction, seed=pars.ranger$seed, quantreg = TRUE)
      gm0 <- ranger(fm0, rmatrix0, mtry=pars.ranger$mtry, min.node.size=pars.ranger$min.node.size, num.trees = pars.ranger$num.trees, sample.fraction=pars.ranger$sample.fraction, seed=pars.ranger$seed)
      ## model to predict mean values
    }
    sel.t = complete.cases(rmatrix.test)
    ## predict s.d. of the prediction error
    x.pred <- predict(gm, rmatrix.test[sel.t,], type="quantiles", quantiles = c((1-.682)/2, 1-(1-.682)/2))$predictions
    pred <- data.frame(predictions=predict(gm0, rmatrix.test[sel.t,])$predictions, se=(x.pred[,2]-x.pred[,1])/2)
    # compute predictive distribution of quantiles given in parameter predDist
    if( !is.null(predDist) ){
      pred.dist <- data.frame( predict(gm, rmatrix.test[sel.t,], type="quantiles", quantiles = predDist )$predictions )
      names(pred.dist) <- paste0("Pred.Quantile.", predDist)
    }
  }
  if(method=="geoR"){
    require(geoR)
    sel.d = complete.cases(over(y=covs, x=s.train))
    x.geo <- as.geodata(s.train[sel.d,varn])
    Range = sqrt(areaSpatialGrid(covs))/3
    locs = s.test@coords
    if(OK==FALSE){
      ## add covariates:
      x.geo$covariate = over(x=s.train[sel.d,varn], y=covs)
      t = as.formula(paste(" ~ ", paste(names(covs), collapse = "+")))
      x.vgm <- likfit(x.geo, trend = t, ini=c(var(log1p(x.geo$data)), Range), fix.psiA = FALSE, fix.psiR = FALSE)
      ov0 = over(s.test[varn], covs)
      sel.t = complete.cases(ov0)
      ov0 = ov0[sel.t,]
      l = as.formula(paste(" ~ ", paste0("ov0$", names(ov0), collapse = "+")))
      KC = krige.control(trend.d = t, trend.l = l, obj.model = x.vgm)
      pred <- krige.conv(x.geo, locations=locs[sel.t,], krige=KC)
    } else {
      x.vgm <- likfit(x.geo, lambda=0, messages=FALSE, ini=c(var(log1p(x.geo$data)), Range), cov.model="exponential")
      pred <- krige.conv(x.geo, locations=locs, krige=krige.control(obj.model=x.vgm))
      sel.t = 1:length(pred$predict)
    }
    names(pred)[1] = "predictions"
    pred$se = sqrt(pred$krige.var)
  }
  obs.pred <- as.data.frame(list(s.test@data[sel.t,varn], pred$predictions, pred$se), col.names=c("Observed", "Predicted", "sdPE"))
  obs.pred[,idcol] <- s.test@data[sel.t,idcol]
  obs.pred$fold = j
  if( !is.null(predDist) & method=="ranger" ){ obs.pred <- cbind(obs.pred, pred.dist) }
  return(obs.pred)
}
