## Functions to use RF to make spatial predictions
## tom.hengl@gmail.com

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

cv_numeric <- function(varn, points, covs, nfold=5, idcol, method="ranger", cpus=1, Nsub=1e4, OK=FALSE, spcT=TRUE, Log=FALSE){
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
    sel <- dismo::kfold(points@data, k=nfold, by=points@data[,idcol])
  } else {
    sel <- dismo::kfold(points@data, k=nfold)
  }
  if(missing(cpus)){ 
    cpus <- parallel::detectCores(all.tests = FALSE, logical = FALSE) 
  }
  if(cpus>1){
    require(snowfall)
    snowfall::sfInit(parallel=TRUE, cpus=ifelse(nfold>cpus, cpus, nfold))
    snowfall::sfExport("predict_parallelP","idcol","points","covs","sel","varn","method","Nsub","OK","spcT")
    snowfall::sfLibrary(package="plyr", character.only=TRUE)
    snowfall::sfLibrary(package="GSIF", character.only=TRUE)
    if(method=="ranger"){
      snowfall::sfLibrary(package="ranger", character.only=TRUE)
      out <- snowfall::sfLapply(1:nfold, function(j){predict_parallelP(j, sel=sel, idcol=idcol, varn=varn, points=points, covs=covs, method=method, cpus=cpus, Nsub=Nsub, OK=OK, spcT=spcT)})
    } else {
      snowfall::sfLibrary(package="geoR", character.only=TRUE)
      out <- snowfall::sfLapply(1:nfold, function(j){predict_parallelP(j, sel=sel, idcol=idcol, varn=varn, points=points, covs=covs, method=method, cpus=cpus, Nsub=Nsub, OK=OK, spcT=spcT)})
    }
    snowfall::sfStop()
  } else {
    out <- lapply(1:nfold, function(j){predict_parallelP(j, sel=sel, idcol=idcol, varn=varn, points=points, covs=covs, method=method, cpus=cpus, Nsub=Nsub, OK=OK, spcT=spcT)})
  }
  ## calculate mean accuracy:
  out <- plyr::rbind.fill(out)
  ME = mean(out$Observed - out$Predicted, na.rm=TRUE) 
  MAE = mean(abs(out$Observed - out$Predicted), na.rm=TRUE)
  RMSE = sqrt(mean((out$Observed - out$Predicted)^2, na.rm=TRUE))
  ## Errors of errors:
  MAE.SE = mean(abs(out$Observed - out$Predicted)-out$SE, na.rm=TRUE)
  ## https://en.wikipedia.org/wiki/Coefficient_of_determination
  #R.squared = 1-sum((out$Observed - out$Predicted)^2, na.rm=TRUE)/(var(out$Observed, na.rm=TRUE)*sum(!is.na(out$Observed)))
  R.squared = 1-var(out$Observed - out$Predicted, na.rm=TRUE)/var(out$Observed, na.rm=TRUE)
  if(Log==TRUE){
    ## If the variable is log-normal then logR.squared is probably more correct
    logRMSE = sqrt(mean((log1p(out$Observed) - log1p(out$Predicted))^2, na.rm=TRUE))
    #logR.squared = 1-sum((log1p(out$Observed) - log1p(out$Predicted))^2, na.rm=TRUE)/(var(log1p(out$Observed), na.rm=TRUE)*sum(!is.na(out$Observed)))
    logR.squared = 1-var(log1p(out$Observed) - log1p(out$Predicted), na.rm=TRUE)/var(log1p(out$Observed), na.rm=TRUE)
    cv.r <- list(out, data.frame(ME=ME, MAE=MAE, RMSE=RMSE, MAE.SE=MAE.SE, R.squared=R.squared, logRMSE=logRMSE, logR.squared=logR.squared)) 
  } else {
    cv.r <- list(out, data.frame(ME=ME, MAE=MAE, RMSE=RMSE, MAE.SE=MAE.SE, R.squared=R.squared))
  }
  names(cv.r) <- c("CV_residuals", "Summary")
  return(cv.r)
}

predict_parallelP <- function(j, sel, idcol, varn, points, covs, method, cpus, Nsub=1e4, OK=FALSE, spcT=TRUE){ 
  s.train <- points[!sel==j,]
  s.test <- points[sel==j,]
  if(missing(Nsub)){ Nsub = length(all.vars(formulaString))*50 }
  if(!Nsub>nrow(s.train)){ 
    s.train = s.train[sample.int(nrow(s.train), Nsub),]
  }
  if(method=="ranger"){
    dist0 <- GSIF::buffer.dist(s.train[varn], covs, as.factor(1:nrow(s.train)))
    dn0 <- paste(names(dist0), collapse="+")
    ovT <- over(s.train[varn], dist0)
    ovV <- over(s.test[varn], dist0)
  }
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
    gm <- ranger(fm0, rmatrix[complete.cases(rmatrix),], keep.inbag = TRUE)
    sel.t = complete.cases(rmatrix.test)
    pred <- predict(gm, rmatrix.test[sel.t,], type = "se")
  }
  if(method=="geoR"){
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
      sel.t = 1:length(pred$predictions)
    }
    names(pred)[1] = "predictions"
    pred$se = sqrt(pred$krige.var)
  }
  obs.pred <- as.data.frame(list(s.test@data[sel.t,varn], pred$predictions, pred$se), col.names=c("Observed", "Predicted", "SE"))
  obs.pred[,idcol] <- s.test@data[sel.t,idcol]
  obs.pred$fold = j
  return(obs.pred)
}