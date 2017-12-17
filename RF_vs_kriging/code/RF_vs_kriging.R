## Comparison RF vs kriging, see slides at: https://github.com/ISRICWorldSoil/GSIF_tutorials/blob/master/geul/5H_Hengl.pdf
## By: tom.hengl@gmail.com, contributions by: Madlene Nussbaum <madlene.nussbaum@bfh.ch> and Marvin Wright <marv@wrig.de>
## Cite as: Hengl et al., "Random Forest as a Generic Framework for Predictive Modeling of Spatial and Spatio-temporal Variables", to be submitted to PeerJ Computer Science
## Licence: GNU GPL

list.of.packages <- c("plyr", "parallel", "randomForest", "quantregForest", "plotKML", "GSIF", "RCurl", "raster", "rgdal", "geoR", "gstat", "scales", "gdistance", "entropy", "lattice", "gridExtra", "intamap", "maxlike", "spatstat")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = TRUE)

setwd("~/git/GeoMLA/RF_vs_kriging")

library(GSIF)
library(rgdal)
library(raster)
library(gstat)
library(plyr)
library(randomForest)
library(plotKML)
library(scales)
library(RCurl)
library(parallel)
library(geoR)
library(lattice)
library(gridExtra)
library(intamap)
library(maxlike)
library(spatstat)
library(entropy)
library(gdistance)

## RANGER connected packages (best install from github):
## speed up of computing of "se" in ranger https://github.com/imbs-hl/ranger/pull/231
#devtools::install_github("imbs-hl/ranger")
library(ranger)
#devtools::install_github("PhilippPro/quantregRanger")
library(quantregRanger)
## http://philipppro.github.io/Tuning_random_forest/
#devtools::install_github("PhilippPro/tuneRF")
library(tuneRF)

source('code/BGUP_functions.R')

## General Settings
## Legend for plots:
leg = c("#0000ff", "#0028d7", "#0050af", "#007986", "#00a15e", "#00ca35", "#00f20d", "#1aff00", "#43ff00", "#6bff00", "#94ff00", "#bcff00", "#e5ff00", "#fff200", "#ffca00", "#ffa100", "#ff7900", "#ff5000", "#ff2800", "#ff0000")
axis.ls = list(at=c(4.8,5.7,6.5,7.4), labels=round(expm1(c(4.8,5.7,6.5,7.4))))
## 1 s.d. quantiles
quantiles = c((1-.682)/2, 0.5, 1-(1-.682)/2)


## ** Example Meuse data set ----------------------------------------------
demo(meuse, echo=FALSE)

## Zinc predicted using ordinary kriging ----
zinc.geo <- as.geodata(meuse["zinc"])
ini.v <- c(var(log1p(zinc.geo$data)),500)
zinc.vgm <- likfit(zinc.geo, lambda=0, ini=ini.v, cov.model="exponential")
locs = meuse.grid@coords
zinc.ok <- krige.conv(zinc.geo, locations=locs, krige=krige.control(obj.model=zinc.vgm))
meuse.grid$zinc_ok = zinc.ok$predict
## Prediction error:
meuse.grid$zinc_ok_var = sqrt(zinc.ok$krige.var)

## Zinc predicted using RF and buffer distances only ----
grid.dist0 <- GSIF::buffer.dist(meuse["zinc"], meuse.grid[1], as.factor(1:nrow(meuse)))
dn0 <- paste(names(grid.dist0), collapse="+")
fm0 <- as.formula(paste("zinc ~ ", dn0))
ov.zinc <- over(meuse["zinc"], grid.dist0)
## Get a good estimate of "mtry" (fine-tuning):
rm.zinc <- cbind(meuse@data["zinc"], ov.zinc)
rt.zinc <- makeRegrTask(data = rm.zinc, target = "zinc")
estimateTimeTuneRF(rt.zinc)
# Make reproducible tuning
set.seed(1)
t.zinc <- tuneRF(rt.zinc, num.trees = 150, build.final.model = FALSE)
t.zinc
## With seed = 1. 
# Recommended parameter settings: 
#   mtry min.node.size sample.fraction
# 1   98             4       0.9259533
# Results: 
#   mse exec.time
# 1 60471.79    0.2716
pars.zinc = list(mtry= t.zinc$recommended.pars$mtry, min.node.size=t.zinc$recommended.pars$min.node.size, sample.fraction=t.zinc$recommended.pars$sample.fraction, num.trees=150, seed = 1)
m.zinc <- quantregRanger(fm0, rm.zinc, params.ranger = pars.zinc)
m.zinc
zinc.rfd <- predict(m.zinc, grid.dist0@data, quantiles)
meuse.grid$zinc_rfd = zinc.rfd[,2]
## Prediction error:
meuse.grid$zinc_rfd_var = (zinc.rfd[,3]-zinc.rfd[,1])/2


## RF and coordinates as covariates only (to be fair, same tuning) ----
fmc <- zinc ~ x + y
rm.zinc.coord <- cbind(meuse@data["zinc"], meuse@coords)
rt.zinc.coord <- makeRegrTask(data = rm.zinc.coord, target = "zinc")
set.seed(1)
t.zinc.coord <- tuneRF(rt.zinc.coord, num.trees = 150, build.final.model = FALSE)
pars.zinc.coord = list(mtry= t.zinc.coord$recommended.pars$mtry, min.node.size=t.zinc.coord$recommended.pars$min.node.size, sample.fraction=t.zinc.coord$recommended.pars$sample.fraction, num.trees=150, seed = 1)
m.zinc.coord <- quantregRanger(fmc, rm.zinc.coord, params.ranger = pars.zinc.coord)
m.zinc.coord
## R-squared still quite high 0.538!
zinc.rfc <- predict(m.zinc.coord, meuse.grid@coords, quantiles)
meuse.grid$zinc_rfc = zinc.rfc[,2]
meuse.grid$zinc_rfc_var = (zinc.rfc[,3]-zinc.rfc[,1])/2

## Plot predictions UK and RF next to each other:
pdf(file = "results/meuse/Fig_comparison_OK_RF_zinc_meuse.pdf", width=9.5, height=8)
var.max = max(c(meuse.grid$zinc_rfd_var, meuse.grid$zinc_ok_var))
var.min = min(c(meuse.grid$zinc_rfd_var, meuse.grid$zinc_ok_var))
par(mfrow=c(2,3), oma=c(0,0,0,1), mar=c(0,0,4,3))
plot(log1p(raster(meuse.grid["zinc_ok"])), col=leg, zlim=c(4.7,7.4), main="Ordinary Kriging (OK)", axes=FALSE, box=FALSE, axis.args=axis.ls)
points(meuse, pch="+")
plot(log1p(raster(meuse.grid["zinc_rfd"])), col=leg, zlim=c(4.7,7.4), main="Random Forest (RF), buffers", axes=FALSE, box=FALSE, axis.args=axis.ls)
points(meuse, pch="+")
plot(log1p(raster(meuse.grid["zinc_rfc"])), col=leg, zlim=c(4.7,7.4), main="Random Forest (RF), coordinates", axes=FALSE, box=FALSE, axis.args=axis.ls)
points(meuse, pch="+")
plot(raster(meuse.grid["zinc_ok_var"]), col=rev(bpy.colors()), zlim=c(var.min,var.max), main="OK prediction error", axes=FALSE, box=FALSE)
points(meuse, pch="+")
plot(raster(meuse.grid["zinc_rfd_var"]), col=rev(bpy.colors()), zlim=c(var.min,var.max), main="RF prediction error, buffers", axes=FALSE, box=FALSE)
points(meuse, pch="+")
plot(raster(meuse.grid["zinc_rfc_var"]), col=rev(bpy.colors()), zlim=c(var.min,var.max), main="RF prediction error, coords", axes=FALSE, box=FALSE)
points(meuse, pch="+")
dev.off()
## TH: RF creates smoother predictions than OK

## Cross-validation Meuse data set ----
cv.RF = cv_numeric(varn="zinc", points=meuse, covs=meuse.grid, cpus=1, method="ranger", OK=TRUE, spcT=FALSE, pars.ranger=pars.zinc)
cv.OK = cv_numeric(varn="zinc", points=meuse, covs=meuse.grid, cpus=1, method="geoR", OK=TRUE, spcT=FALSE)
cv.RF$Summary$RMSE^2/var(meuse$zinc, na.rm = TRUE); cv.RF$Summary$R.squared; cv.RF$Summary$ZSV; cv.RF$Summary$MAE.SE
cv.OK$Summary$RMSE^2/var(meuse$zinc, na.rm = TRUE); cv.OK$Summary$R.squared; cv.OK$Summary$ZSV; cv.OK$Summary$MAE.SE

## plot results
lim.zinc = range(meuse$zinc, na.rm = TRUE)
pdf(file = "results/meuse/Fig_correlation_plots_OK_RF_zinc_meuse.pdf", width=9, height=5)
par(oma=c(0,0,0,1), mar=c(0,0,0,2))
pfun.loess <- function(x,y, ...){ 
  panel.xyplot(x,y, ...)
  panel.abline(0,1,lty=2,lw=1,col="grey60")
  panel.loess(x,y,span=0.5, col = "grey30", lwd = 1.7)
}
plt.RF = xyplot(cv.RF[[1]]$Observed~cv.RF[[1]]$Predicted, asp=1, par.settings=list(plot.symbol = list(col=alpha("black", 0.6), fill=alpha("red", 0.6), pch=21, cex=0.9)), scales=list(x=list(log=TRUE, equispaced.log=FALSE), y=list(log=TRUE, equispaced.log=FALSE)), ylab="observed zinc [mg/kg]", xlab="znic predicted by RF [mg/kg]", panel = pfun.loess, xlim=lim.zinc, ylim=lim.zinc)
plt.OK = xyplot(cv.OK[[1]]$Observed~cv.OK[[1]]$Predicted, asp=1, par.settings=list(plot.symbol = list(col=alpha("black", 0.6), fill=alpha("red", 0.6), pch=21, cex=0.9)), scales=list(x=list(log=TRUE, equispaced.log=FALSE), y=list(log=TRUE, equispaced.log=FALSE)), ylab="observed zinc [mg/kg]", xlab="zinc predicted by geoR [mg/kg]", panel = pfun.loess, xlim=lim.zinc, ylim=lim.zinc)
grid.arrange(plt.OK, plt.RF, ncol=2)
dev.off()

## RF with combined covariates ----
meuse.grid$SW_occurrence = readGDAL("data/meuse/Meuse_GlobalSurfaceWater_occurrence.tif")$band1[meuse.grid@grid.index]
meuse.grid$AHN = readGDAL("data/meuse/ahn.asc")$band1[meuse.grid@grid.index]
meuse.grid$LGN5 = as.factor(readGDAL("data/meuse/lgn5.asc")$band1[meuse.grid@grid.index])
grids.spc = spc(meuse.grid, as.formula("~ SW_occurrence + AHN + ffreq + dist"))
## fit hybrid RF model:
fm1 <- as.formula(paste("zinc ~ ", dn0, " + ", paste(names(grids.spc@predicted), collapse = "+")))
ov.zinc1 <- over(meuse["zinc"], grids.spc@predicted)
m1.zinc <- quantregRanger(fm1, do.call(cbind, list(meuse@data["zinc"], ov.zinc, ov.zinc1)), params.ranger = c(list(importance = "impurity"), pars.zinc))
m1.zinc
zinc.rfd1 <- predict(m1.zinc, cbind(grid.dist0@data, grids.spc@predicted@data), quantiles)
meuse.grid$zinc_rfd1 = zinc.rfd1[,2]
meuse.grid$zinc_rfd1_var = (zinc.rfd1[,3]-zinc.rfd1[,1])/2
xl <- as.list(m1.zinc$variable.importance)
print(t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:10]])))
m2.zinc <- quantregRanger(paste("zinc ~ ", paste(names(grids.spc@predicted), collapse = "+")), do.call(cbind, list(meuse@data["zinc"], ov.zinc1)))
m2.zinc
meuse.grid$zinc_rfd2 = predict(m2.zinc, grids.spc@predicted@data)[,2]

pdf(file = "results/meuse/Fig_RF_covs_bufferdist_zinc_meuse.pdf", width=9, height=5)
par(mfrow=c(1,2), oma=c(0,0,0,1), mar=c(0,0,4,3))
plot(log1p(raster(meuse.grid["zinc_rfd2"])), col=leg, zlim=c(4.7,7.4), main="Random Forest (RF) covs only", axes=FALSE, box=FALSE, axis.args=axis.ls)
points(meuse, pch="+")
plot(log1p(raster(meuse.grid["zinc_rfd1"])), col=leg, zlim=c(4.7,7.4), main="Random Forest (RF) covs + buffer dist.", axes=FALSE, box=FALSE, axis.args=axis.ls)
points(meuse, pch="+")
dev.off()


## SIC 1997 data set ----
## measurements made in Switzerland on the 8th of May 1986
sic97.sp = readRDS("data/rainfall/sic97.rds")
swiss1km = readRDS("data/rainfall/swiss1km.rds")
ov2 = over(y=swiss1km, x=sic97.sp)
sel.d = which(!is.na(ov2$DEM))
#plot(stack(swiss1km[1:2]))
## linear geostatistical model from: https://www.jstatsoft.org/article/view/v063i12/v63i12.pdf 
sic97.geo <- as.geodata(sic97.sp[sel.d,"rainfall"])
## add covariates:
sic97.geo$covariate = ov2[sel.d,c("CHELSA_rainfall","DEM")]
#plot(variog4(sic97.geo, lambda=0.5, max.dist=65000, messages=FALSE), lwd=2)
sic.t = ~ CHELSA_rainfall + DEM
rain.vgm <- likfit(sic97.geo, trend = sic.t, ini=c(var(log1p(sic97.geo$data)),8000), fix.psiA = FALSE, fix.psiR = FALSE)
rain.vgm
locs2 = swiss1km@coords
KC = krige.control(trend.d = sic.t, trend.l = ~ swiss1km$CHELSA_rainfall + swiss1km$DEM, obj.model = rain.vgm)
rain.uk <- krige.conv(sic97.geo, locations=locs2, krige=KC)
swiss1km$rainfall_UK = rain.uk$predict
swiss1km$rainfall_UK_var = rain.uk$krige.var
#plot(raster(swiss1km["rainfall_UK"]))
#plot(sqrt(raster(swiss1km["rainfall_UK_var"])))

## SIC 1997 Random Forest example ----
swiss.dist0 <- buffer.dist(sic97.sp["rainfall"], swiss1km[1], as.factor(1:nrow(sic97.sp))) ## takes 2-3 mins!
## Deriving buffer distances is computationally very intensive
ov.swiss = over(sic97.sp["rainfall"], swiss.dist0)
sw.dn0 <- paste(names(swiss.dist0), collapse="+")
sw.fm1 <- as.formula(paste("rainfall ~ ", sw.dn0, " + CHELSA_rainfall + DEM"))
ov.rain <- over(sic97.sp["rainfall"], swiss1km[1:2])
sw.rm = do.call(cbind, list(sic97.sp@data["rainfall"], ov.rain, ov.swiss))
## fine-tune RF:
#rt.rain <- makeRegrTask(data = sw.rm[complete.cases(sw.rm[,all.vars(sw.fm1)]),], target = "rainfall")
#estimateTimeTuneRF(rt.rain)
## Too time-consuming
#t.rain <- tuneRF(rt.rain, num.trees = 100, build.final.model = FALSE)
#t.rain
#pars.rain = list(mtry= t.rain$recommended.pars$mtry, min.node.size=t.rain$recommended.pars$min.node.size, sample.fraction=t.rain$recommended.pars$sample.fraction, num.trees=100)
m1.rain <- quantregRanger(sw.fm1, sw.rm[complete.cases(sw.rm),], params.ranger = list(importance = "impurity", mtry=140, num.trees=150))
m1.rain
## 0.81
rain.rfd1 <- predict(m1.rain, cbind(swiss.dist0@data, swiss1km@data), quantiles)
swiss1km$rainfall_rfd1 = rain.rfd1[,2]
swiss1km$rainfall_rfd1_var = (rain.rfd1[,3]-rain.rfd1[,1])/2
xl1 <- as.list(m1.rain$variable.importance)
print(t(data.frame(xl1[order(unlist(xl1), decreasing=TRUE)[1:15]])))

## Plot predictions next to each other:
rain.max = max(swiss1km$rainfall_rfd1, na.rm = TRUE)
rainv.max = max(swiss1km$rainfall_rfd1_var, na.rm = TRUE)
swiss1km$rainfall_UK = ifelse(swiss1km$rainfall_UK<0, 0, ifelse(swiss1km$rainfall_UK>rain.max, rain.max, swiss1km$rainfall_UK))
pdf(file = "results/rainfall/Fig_Swiss_rainfall_UK_vs_RF.pdf", width=12, height=8)
par(mfrow=c(2,2), oma=c(0,0,0,0.5), mar=c(0,0,1.5,1))
plot(raster(swiss1km["rainfall_UK"]), col=leg, main="Universal kriging (UK)", axes=FALSE, box=FALSE, zlim=c(0, rain.max))
points(sic97.sp, pch="+")
plot(raster(swiss1km["rainfall_rfd1"]), col=leg, main="Random Forest (RF)", axes=FALSE, box=FALSE, zlim=c(0, rain.max))
points(sic97.sp, pch="+")
plot(sqrt(raster(swiss1km["rainfall_UK_var"])), col=rev(bpy.colors()), main="Universal kriging (UK) prediction error", axes=FALSE, box=FALSE, zlim=c(0,rainv.max))
points(sic97.sp, pch="+")
plot(raster(swiss1km["rainfall_rfd1_var"]), col=rev(bpy.colors()), main="Random Forest (RF) prediction error", axes=FALSE, box=FALSE, zlim=c(0,rainv.max))
points(sic97.sp, pch="+")
dev.off()

## Cross-validation SIC 1997 ----
## (computationally intensive - takes few minutes!!)
cv.RF2 = cv_numeric(varn="rainfall", points=sic97.sp, covs=swiss1km[c("CHELSA_rainfall","DEM")], cpus=1, method="ranger", spcT=FALSE, pars.ranger = list(mtry=140, num.trees=150))
cv.UK = cv_numeric(varn="rainfall", points=sic97.sp, covs=swiss1km[c("CHELSA_rainfall","DEM")], cpus=1, method="geoR", spcT=FALSE)
cv.RF2$Summary$RMSE^2/var(sic97.sp$rainfall); cv.RF2$Summary$R.squared; cv.RF2$Summary$ZSV; cv.RF2$Summary$MAE.SE
cv.UK$Summary$RMSE^2/var(sic97.sp$rainfall); cv.UK$Summary$R.squared;  cv.UK$Summary$ZSV; cv.UK$Summary$MAE.SE
sqrt(var(sic97.sp$rainfall)/nrow(sic97.sp))
sd(sic97.sp$rainfall)

## Plot residuals vgm:
x2 = plyr::join_all(list(data.frame(r1=cv.RF2$CV_residuals$Observed-cv.RF2$CV_residuals$Predicted, id=cv.RF2$CV_residuals$SOURCEID), data.frame(r2=cv.UK$CV_residuals$Observed-cv.UK$CV_residuals$Predicted, id=cv.UK$CV_residuals$SOURCEID)))
x2 = plyr::join(x2, data.frame(id=row.names(sic97.sp@data), rain=sic97.sp$rainfall, x=sic97.sp@coords[,1], y=sic97.sp@coords[,2]))
plot_vgm(rain~1, x2, swiss1km, r1="r1", r2="r2", main="Rainfall (SIC 1997)")

## plot results
lim.rain = c(20,max(sic97.sp$rainfall, na.rm = TRUE))
pdf(file = "results/rainfall/Fig_correlation_plots_OK_RF_rain_SIC97.pdf", width=9, height=5)
par(oma=c(0,0,0,1), mar=c(0,0,0,2))
plt.RF2 = xyplot(cv.RF2[[1]]$Observed~cv.RF2[[1]]$Predicted, asp=1, par.settings=list(plot.symbol = list(col=alpha("black", 0.6), fill=alpha("red", 0.6), pch=21, cex=0.9)), scales=list(x=list(log=TRUE, equispaced.log=FALSE), y=list(log=TRUE, equispaced.log=FALSE)), ylab="observed", xlab="predicted (machine learning)", panel = pfun.line, xlim=lim.rain, ylim=lim.rain)
plt.UK = xyplot(cv.UK[[1]]$Observed~cv.UK[[1]]$Predicted, asp=1, par.settings=list(plot.symbol = list(col=alpha("black", 0.6), fill=alpha("red", 0.6), pch=21, cex=0.9)), scales=list(x=list(log=TRUE, equispaced.log=FALSE), y=list(log=TRUE, equispaced.log=FALSE)), ylab="observed", xlab="predicted (geoR)", panel = pfun.line, xlim=lim.rain, ylim=lim.rain)
grid.arrange(plt.UK, plt.RF2, ncol=2)
dev.off()

## Intamap example ----
data(sic2004)
coordinates(sic.val) <- ~x+y
sic.val$value <- sic.val$joker
#sic.val$value <- sic.val$dayx
#writeOGR(sic.val, "results/sic2004/sic.val.shp", "sic.val", "ESRI Shapefile")
coordinates(sic.test) <- ~x+y
pred.sic2004 <- interpolate(sic.val, sic.test, maximumTime = 90)
#R 2017-12-06 15:00:58 interpolating 200 observations, 808 prediction locations
#[1] "estimated time for copula 71.2993216100125"
spplot(pred.sic2004$predictions[1])
#plot(sic.test$dayx~pred.sic2004$predictions$mean, asp=1)
#sd(sic.test$dayx-pred.sic2004$predictions$mean)
## 12.4
sd(sic.test$joker-pred.sic2004$predictions$mean)
## 104

## RFsp
library(tuneRF)
library(quantregRanger)
bbox=sic.val@bbox
bbox[,"min"]=bbox[,"min"]-4000
bbox[,"max"]=bbox[,"max"]+4000
de2km = plotKML::vect2rast(sic.val, cell.size=2000, bbox=bbox)
de2km$mask = 1
de2km = as(de2km["mask"], "SpatialPixelsDataFrame")
plot(de2km); points(sic.val)
hist(sic.val$joker, breaks=45, col="grey")
which(sic.val$joker>500) ## only 2 points with very high values
which(sic.test$joker>500)
de.dist0 <- GSIF::buffer.dist(sic.val["joker"], de2km, as.factor(1:nrow(sic.val@data)))
ov.de = over(sic.val["joker"], de.dist0)
de.dn0 <- paste(names(de.dist0), collapse="+")
de.fm1 <- as.formula(paste("joker ~ ", de.dn0))
de.rm = do.call(cbind, list(sic.val@data["joker"], ov.de))
## fine-tuning recommended since the variable has highly skewed distribution with 2 hot spots:
rt.gamma <- makeRegrTask(data = de.rm[complete.cases(de.rm[,all.vars(de.fm1)]),], target = "joker")
estimateTimeTuneRF(rt.gamma)
## 8M
t.gamma <- tuneRF(rt.gamma, build.final.model = FALSE)
t.gamma
pars.gamma = list(mtry=t.gamma$recommended.pars$mtry, min.node.size=t.gamma$recommended.pars$min.node.size, sample.fraction=t.gamma$recommended.pars$sample.fraction)
m1.gamma <- quantregRanger(de.fm1, de.rm[complete.cases(de.rm),], params.ranger = pars.gamma)
m1.gamma
## R squared (OOB): 0.11
gamma.rfd1 <- predict(m1.gamma, de.dist0@data, quantiles)
de2km$gamma_rfd1 = gamma.rfd1[,2]
de2km$gamma_rfd1_var = (gamma.rfd1[,3]-gamma.rfd1[,1])/2
plot(de2km["gamma_rfd1"])
points(sic.val, pch="+")
#plot(de2km["gamma_rfd1_var"])
#pred2km.sic2004 <- interpolate(sic.val, de2km, methodName = "automap")
#de2km$gamma_ok = pred2km.sic2004$predictions$var1.pred
#plot(de2km["gamma_ok"])
#points(sic.val, pch="+")
ov.test <- over(sic.test, de2km["gamma_rfd1"])
#plot(sic.test$dayx~ov.test$gamma_rfd1, asp=1)
sd(sic.test$joker-ov.test$gamma_rfd1, na.rm=TRUE)

## Weighted regression RF ----
carson <- read.csv(file="data/NRCS/carson_CLYPPT.csv")
summary(carson$CLYPPT)
carson$DEPTH.f = ifelse(is.na(carson$DEPTH), 20, carson$DEPTH)
carson1km <- readRDS("data/NRCS/carson_covs1km.rds")
coordinates(carson) <- ~X+Y
proj4string(carson) = carson1km@proj4string
rm.carson <- cbind(as.data.frame(carson), over(carson["CLYPPT"], carson1km))
fm.clay <- as.formula(paste("CLYPPT ~ DEPTH.f + ", paste(names(carson1km), collapse = "+")))
rm.carson <- rm.carson[complete.cases(rm.carson[,all.vars(fm.clay)]),]
rm.carson.s <- rm.carson[sample.int(size=1250, nrow(rm.carson)),]
m.clay <- quantregRanger(fm.clay, rm.carson.s, params.ranger = list(num.trees=150, mtry=25, case.weights=1/(rm.carson.s$CLYPPT.sd^2)))
m.clay
carson1km$DEPTH.f = 30
clay.rfd <- predict(m.clay, carson1km@data, quantiles)
carson1km$clay_rfd = ifelse(clay.rfd[,2]<10, 10, ifelse(clay.rfd[,2]>35, 35, clay.rfd[,2]))
summary(carson1km$clay_rfd)
carson1km$clay_rfd_var = (clay.rfd[,3]-clay.rfd[,1])/2
#plot(raster(carson1km["clay_rfd"]), col=plotKML::SAGA_pal[[1]])
m.clay0 <- quantregRanger(fm.clay, rm.carson.s, params.ranger = list(num.trees=150, mtry=25))
m.clay0
clay0.rfd <- predict(m.clay0, carson1km@data, quantiles)
carson1km$clay0_rfd = ifelse(clay0.rfd[,2]<10, 10, ifelse(clay0.rfd[,2]>35, 35, clay0.rfd[,2]))
carson1km$clay0_rfd_var = (clay0.rfd[,3]-clay0.rfd[,1])/2

pdf(file = "results/NRCS/Fig_clay_RF_weighted.pdf", width=10, height=9)
par(mfrow=c(2,2), oma=c(0,0,0,0.5), mar=c(0,0,1.5,1))
par(oma=c(0,0,0,0.5), mar=c(0,0,3.5,1))
plot(raster(carson1km["clay_rfd"]), col=leg, main="RF predictions clay (weigthed)", axes=FALSE, box=FALSE, zlim=c(10,35))
points(rm.carson.s$X, rm.carson.s$Y, pch="+", cex=.8)
plot(raster(carson1km["clay0_rfd"]), col=leg, main="RF predictions clay", axes=FALSE, box=FALSE, zlim=c(10,35))
points(rm.carson.s$X, rm.carson.s$Y, pch="+", cex=.8)
plot(raster(carson1km["clay_rfd_var"]), col=rev(bpy.colors()), main="Prediction error (RF weigthed)", axes=FALSE, box=FALSE, zlim=c(0,28))
points(rm.carson.s$X, rm.carson.s$Y, pch="+", cex=.8)
plot(raster(carson1km["clay0_rfd_var"]), col=rev(bpy.colors()), main="Prediction error (RF)", axes=FALSE, box=FALSE, zlim=c(0,28))
points(rm.carson.s$X, rm.carson.s$Y, pch="+", cex=.8)
dev.off()



## Ebergotzen binomial variable ----
data(eberg)
eb.s = sample.int(nrow(eberg), 1200)
eberg = eberg[eb.s,]
coordinates(eberg) <- ~X+Y
proj4string(eberg) <- CRS("+init=epsg:31467")
## Binomial variable:
eberg$Parabraunerde <- ifelse(eberg$TAXGRSC=="Parabraunerde", "TRUE", "FALSE")
data(eberg_grid)
gridded(eberg_grid) <- ~x+y
proj4string(eberg_grid) <- CRS("+init=epsg:31467")
eberg_spc <- spc(eberg_grid, ~ PRMGEO6+DEMSRT6+TWISRT6+TIRAST6)
eberg_grid@data <- cbind(eberg_grid@data, eberg_spc@predicted@data)
summary(as.factor(eberg$Parabraunerde))
## overlay and create a regression-matrix:
ov.eberg <- over(eberg, eberg_grid)
sel.eberg = !is.na(ov.eberg$DEMSRT6)
rm.eberg <- cbind(ov.eberg, eberg@data)
summary(rm.eberg$TAXGRSC)
summary(as.factor(rm.eberg$Parabraunerde))
eberg.dist0 <- buffer.dist(eberg[sel.eberg,"Parabraunerde"], eberg_grid[2], as.factor(1:sum(sel.eberg))) ## takes 2-3 mins
ov.eberg2 = over(eberg[sel.eberg,"Parabraunerde"], eberg.dist0)
eb.dn0 <- paste(names(eberg.dist0), collapse="+")
eb.fm1 <- as.formula(paste("Parabraunerde ~ ", eb.dn0, "+", paste0("PC", 1:10, collapse = "+")))
ov.eberg3 <- over(eberg[sel.eberg,"Parabraunerde"], eberg_grid[paste0("PC", 1:10)])
rm.eberg2 = do.call(cbind, list(eberg@data[sel.eberg,c("Parabraunerde","TAXGRSC")], ov.eberg2, ov.eberg3))
m1.Parabraunerde <- ranger(eb.fm1, rm.eberg2[complete.cases(rm.eberg2),], importance = "impurity", probability = TRUE)
m1.Parabraunerde
## Seems to be quite an accurate model
xl1.P <- as.list(ranger::importance(m1.Parabraunerde))
print(t(data.frame(xl1.P[order(unlist(xl1.P), decreasing=TRUE)[1:10]])))
## Predict:
pr.Parabraunerde = predict(m1.Parabraunerde, cbind(eberg.dist0@data, eberg_grid@data[paste0("PC", 1:10)]))
eberg_grid$Parabraunerde_TRUE = pr.Parabraunerde$predictions[,2]
eberg_grid$Parabraunerde_FALSE = pr.Parabraunerde$predictions[,1]
eberg_grid$SSE = entropy_index(eberg_grid@data[,c("Parabraunerde_TRUE","Parabraunerde_FALSE")])

pdf(file = "results/eberg/Fig_Parabraunerde_RF.pdf", width=7, height=6.5)
#par(mfrow=c(1,2), oma=c(0,0,0,0.5), mar=c(0,0,1.5,1))
par(oma=c(0,0,0,0.5), mar=c(0,0,3.5,1))
plot(raster(eberg_grid["Parabraunerde_TRUE"]), col=SAGA_pal[["SG_COLORS_YELLOW_RED"]], zlim=c(0,1), main="Parabraunerde class (RF)", axes=FALSE, box=FALSE)
points(eberg[eberg$Parabraunerde=="TRUE"&!is.na(eberg$Parabraunerde)&sel.eberg,], pch=19, cex=.5)
points(eberg[eberg$Parabraunerde=="FALSE"&!is.na(eberg$Parabraunerde)&sel.eberg,], pch="+", cex=.8)
#plot(raster(eberg_grid["SSE"]), col=rev(bpy.colors()), main="Shannon Scaled Entropy Index", axes=FALSE, box=FALSE)
#points(eberg[eberg$Parabraunerde=="TRUE"&!is.na(eberg$Parabraunerde)&sel.eberg,], pch=19, cex=.8)
#points(eberg[eberg$Parabraunerde=="FALSE"&!is.na(eberg$Parabraunerde)&sel.eberg,], pch="+", cex=.8)
dev.off()

## predict ALL soil types:
data(eberg)
#eb.s = sample.int(nrow(eberg), 1200)
#eberg = eberg[eb.s,]
coordinates(eberg) <- ~X+Y
proj4string(eberg) <- CRS("+init=epsg:31467")
## Fit model and predict at once:
soiltype <- autopredict(eberg["TAXGRSC"], eberg_grid, auto.plot=FALSE)
pdf("results/eberg/Fig_ebergotzen_TAXGRSC.pdf", width=10, height=6.7)
plot(stack(soiltype$predicted), col=SAGA_pal[["SG_COLORS_YELLOW_RED"]], zlim=c(0,100))
dev.off()

## Plot soil type "G"
r.G = soiltype$predicted["Gley"]
r.G$Gley = ifelse(r.G$Gley>40, 40, r.G$Gley)
plot(raster(r.G), col=SAGA_pal[["SG_COLORS_YELLOW_RED"]], zlim=c(0,40))
points(eberg[eberg$TAXGRSC=="Gley"&!is.na(eberg$TAXGRSC),], pch=19)
points(eberg[!eberg$soiltype=="Gley"&!is.na(eberg$TAXGRSC),], pch="+")
r.P = soiltype$predicted["Parabraunerde"]
r.P$Parabraunerde = ifelse(r.P$Parabraunerde>40, 40, r.P$Parabraunerde)
eberg_grid$SSE_t = entropy_index(soiltype$predicted@data)
ov.eberg <- over(eberg, eberg_grid)
sel.eberg = !is.na(ov.eberg$DEMSRT6)
plot(eberg_grid["SSE_t"], col=rev(bpy.colors()), zlim=c(0,100))
points(eberg[sel.eberg,"TAXGRSC"], pch="+", cex=.6)

shape = "http://maps.google.com/mapfiles/kml/pal2/icon18.png" 
kml(eberg, colour = TAXGRSC, file.name="results/eberg/eberg_TAXGRSC.kml", shape=shape, points_names=eberg$TAXGRSC, LabelScale = .9)
kml(raster(r.G), folder.name="Gley", raster_name="results/eberg/eberg_Gley.png", file.name="results/eberg/eberg_Gley.kml", colour_scale=SAGA_pal[["SG_COLORS_YELLOW_RED"]], zlim=c(0,40), png.type="cairo")
kml(raster(r.P), folder.name="Parabraunerde", raster_name="results/eberg/eberg_Parabraunerde.png", file.name="results/eberg/eberg_Parabraunerde.kml", colour_scale=SAGA_pal[["SG_COLORS_YELLOW_RED"]], zlim=c(0,40), png.type="cairo")
kml(eberg_grid["SSE_t"], folder.name="SSE", raster_name="results/eberg/eberg_SEE.png", file.name="results/eberg/eberg_SEE.kml", colour_scale=rev(bpy.colors()), zlim=c(0,100), png.type="cairo")
## Conclusion: looks like regression-kriging on class probs

## Ebergotzen weighted RF (usually measurement error or sampling probability) ----
## Estimate occurrence probability
data(eberg)
coordinates(eberg) <- ~X+Y
proj4string(eberg) <- CRS("+init=epsg:31467")
ov.eberg <- over(eberg, eberg_grid)
sel.eberg = !is.na(ov.eberg$DEMSRT6)
eberg.xy <- eberg[sel.eberg,"CLYMHT_B"]
sprob <- spsample.prob(eberg.xy, eberg_grid[paste0("PC", 1:10)])
ov.ebergS <- over(eberg.xy, eberg_grid[paste0("PC", 1:10)])
rm.ebergS <- do.call(cbind, list(eberg.xy@data, ov.ebergS, over(eberg.xy, sprob[[1]])))
fs <- as.formula(paste("CLYMHT_B ~ ", paste(paste0("PC", 1:10), collapse="+")))
## the lower the occurrence probability, the higher the weight:
rm.ebergS <- rm.ebergS[complete.cases(rm.ebergS[,all.vars(fs)]),]
m.CLYw <- quantregRanger(fs, rm.ebergS, list(importance='impurity', case.weights=1/rm.ebergS$iprob))
m.CLYw
## R squared (OOB): 0.53
m.CLYn <- quantregRanger(fs, rm.ebergS, list(importance='impurity'))
## Compare predictions:
CLY.rfw <- predict(m.CLYw, eberg_grid@data, quantiles)
eberg_grid$CLY_rfw <- CLY.rfw[,2]
summary(eberg_grid$CLY_rfw)
eberg_grid$CLY_rfw_var <- (CLY.rfw[,3]-CLY.rfw[,1])/2
CLY.rfn <- predict(m.CLYn, eberg_grid@data, quantiles)
eberg_grid$CLY_rfn <- CLY.rfn[,2]
eberg_grid$CLY_rfn_var <- (CLY.rfn[,3]-CLY.rfn[,1])/2
## compare with pure random sampling
rnd <- spsample(eberg_grid, type="random", n=length(sprob[["observations"]]))
sprob2 <- spsample.prob(rnd, eberg_grid[paste0("PC", 1:10)])

pdf(file = "results/eberg/Fig_clay_RF_weighted.pdf", width=11, height=6.5)
par(mfrow=c(2,3), oma=c(0,0,0,0.5), mar=c(0,0,1.5,1))
par(oma=c(0,0,0,0.5), mar=c(0,0,3.5,1))
plot(raster(sprob[[1]]), col=SAGA_pal[["SG_COLORS_YELLOW_RED"]], zlim=c(0,1), main="Occurrence probability", axes=FALSE, box=FALSE)
points(eberg.xy, pch="+", cex=.8)
plot(raster(eberg_grid["CLY_rfn"]), col=leg, main="RF predictions clay", axes=FALSE, box=FALSE, zlim=c(0,55))
plot(raster(eberg_grid["CLY_rfw"]), col=leg, main="RF predictions clay (weigthed)", axes=FALSE, box=FALSE, zlim=c(0,55))
plot(raster(sprob2[[1]]), col=SAGA_pal[["SG_COLORS_YELLOW_RED"]], zlim=c(0,1), main="Occurrence probability", axes=FALSE, box=FALSE)
points(rnd, pch="+", cex=.8)
plot(raster(eberg_grid["CLY_rfn_var"]), col=rev(bpy.colors()), main="Prediction error (RF)", axes=FALSE, box=FALSE, zlim=c(0,28))
plot(raster(eberg_grid["CLY_rfw_var"]), col=rev(bpy.colors()), main="Prediction error (RF weigthed)", axes=FALSE, box=FALSE, zlim=c(0,28))
dev.off()

## Deriving more complex distances ----
#writeGDAL(eberg_grid["DEMSRT6"], "/data/tmp/DEMSRT6.sdat", "SAGA")
r <- raster(eberg_grid["DEMSRT6"])
hd <- transition(r, function(x){x[2] - x[1]}, 8, symm=FALSE)
slope <- geoCorrection(hd)
adj <- adjacent(r, cells=1:ncell(r), pairs=TRUE, directions=8)
speed <- slope
speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] + 0.05))
T <- geoCorrection(speed)
acost <- accCost(T, c(3575290.6,5713305.5))
plot(acost)
writeRaster(acost, "/data/tmp/acost.sdat", "SAGA", overwrite=TRUE)

## Multivariate case ----
## Geochemicals USA (https://mrdata.usgs.gov/geochem/)
geochem = readRDS("data/geochem/geochem.rds")
## negative values are in fact detection limits:
for(i in c("PB_ICP40","CU_ICP40","K_ICP40","MG_ICP40")) { geochem[,i] = ifelse(geochem[,i] < 0, abs(geochem[,i])/2, geochem[,i])  }
coordinates(geochem) = ~coords.x1 + coords.x2
proj4string(geochem) = "+proj=longlat +ellps=clrk66 +towgs84=-9.0,151.0,185.0,0.0,0.0,0.0,0.0 +no_defs"
bubble(geochem[!is.na(geochem$PB_ICP40),"PB_ICP40"])
bubble(geochem[!is.na(geochem$CU_ICP40),"CU_ICP40"])
bubble(geochem[!is.na(geochem$K_ICP40),"K_ICP40"])
bubble(geochem[!is.na(geochem$MG_ICP40),"MG_ICP40"])
geochem$TYPEDESC = as.factor(paste(geochem$TYPEDESC))
summary(geochem$TYPEDESC)
#writeOGR(geochem, "geochem.shp", "geochem", "ESRI Shapefile")
usa5km = readRDS("data/geochem/usa5km.rds")
str(usa5km@data)
geochem = spTransform(geochem, CRS(proj4string(usa5km)))
usa5km.spc = spc(usa5km, ~geomap+globedem+dTRI+nlights03+dairp+sdroads)
ov.geochem = over(x=geochem, y=usa5km.spc@predicted)
t.vars = c("PB_ICP40","CU_ICP40","K_ICP40","MG_ICP40")
df.lst = lapply(t.vars, function(i){cbind(geochem@data[,c(i,"TYPEDESC")], ov.geochem)})
names(df.lst) = t.vars
for(i in t.vars){colnames(df.lst[[i]])[1] = "Y"}
for(i in t.vars){df.lst[[i]]$TYPE = i}
## All variables now with the same name:
rm.geochem = do.call(rbind, df.lst)
#str(rm.geochem)
hist(log1p(rm.geochem$Y))
summary(as.factor(rm.geochem$TYPE))
type.mat = data.frame(model.matrix(~TYPE-1, rm.geochem))
typed.mat = data.frame(model.matrix(~TYPEDESC-1, rm.geochem))
## Final regression matrix
rm.geochem.e = do.call(cbind, list(rm.geochem[,c("Y",paste0("PC",1:21))], type.mat, typed.mat))
fm.g = as.formula(paste0("Y ~ ", paste0(names(rm.geochem.e)[-1], collapse = "+")))
fm.g
m1.geochem <- ranger::ranger(fm.g, rm.geochem.e[complete.cases(rm.geochem.e),], importance = "impurity")
m1.geochem
xl1.g <- as.list(ranger::importance(m1.geochem))
print(t(data.frame(xl1.g[order(unlist(xl1.g), decreasing=TRUE)[1:15]])))
## Predict concentrations
new.usa5km = usa5km.spc@predicted@data
new.usa5km$TYPEDESCSOIL = 0
new.usa5km$TYPEDESCSTRM.SED.DRY = 0
new.usa5km$TYPEDESCSTRM.SED.WET = 1
new.usa5km$TYPEDESCUNKNOWN = 0
for(i in t.vars){
  new.usa5km[,paste0("TYPE",i)] = 1
  for(j in t.vars[!t.vars %in% i]){ new.usa5km[,paste0("TYPE",j)] = 0 }
  x <- predict(m1.geochem, new.usa5km)
  ## fix display options for very skewed variables
  if(any(i %in% c("PB_ICP40","CU_ICP40"))){
    usa5km@data[,paste0(i,"_rf")] = ifelse(x$predictions  < min(geochem@data[,i], na.rm=TRUE), min(geochem@data[,i], na.rm=TRUE), ifelse(x$predictions > quantile(geochem@data[,i], .975, na.rm=TRUE), quantile(geochem@data[,i], .975, na.rm=TRUE), x$predictions)) } else {
    usa5km@data[,paste0(i,"_rf")] = ifelse(x$predictions  < min(geochem@data[,i], na.rm=TRUE), min(geochem@data[,i], na.rm=TRUE), ifelse(x$predictions > max(geochem@data[,i], na.rm=TRUE), max(geochem@data[,i], na.rm=TRUE), x$predictions)) 
  }
}

## Plot predictions next to each other:
pdf(file = "results/geochem/Fig_NGS_elements_RF.pdf", width=7.5, height=7)
par(mfrow=c(2,2), oma=c(0,0,0,0.5), mar=c(0,0,2.5,2))
plot(raster(usa5km["PB_ICP40_rf"]), col=leg, main="Pb (ppm)", axes=FALSE, box=FALSE)
points(geochem[!is.na(geochem$PB_ICP40),], pch="+", cex=.5)
plot(raster(usa5km["CU_ICP40_rf"]), col=leg, main="Cu (ppm)", axes=FALSE, box=FALSE)
points(geochem[!is.na(geochem$CU_ICP40),], pch="+", cex=.5)
plot(raster(usa5km["K_ICP40_rf"]), col=leg, main="K (wt%)", axes=FALSE, box=FALSE)
points(geochem[!is.na(geochem$K_ICP40),], pch="+", cex=.5)
plot(raster(usa5km["MG_ICP40_rf"]), col=leg, main="Mg (wt%)", axes=FALSE, box=FALSE)
points(geochem[!is.na(geochem$MG_ICP40),], pch="+", cex=.5)
dev.off()


## Spatiotemporal prediction using Random Forest ----
## Daily precipitation obtained from https://www.ncdc.noaa.gov/cdo-web/search
## compare with: https://www.r-bloggers.com/part-4a-modelling-predicting-the-amount-of-rain/
## NCDC data access explained at: http://neondataskills.org/R/COOP-precip-data-R
co_prec = readRDS("data/st_prec/boulder_prcp.rds")
str(co_prec)
# 'data.frame':	176467 obs. of  16 variables:
summary(co_prec$PRCP)
## Derive cumulative day:
co_prec$cdate = floor(unclass(as.POSIXct(as.POSIXct(paste(co_prec$DATE), format="%Y-%m-%d")))/86400)
hist(co_prec$cdate)
## Day of the year:
co_prec$doy = as.integer(strftime(as.POSIXct(paste(co_prec$DATE), format="%Y-%m-%d"), format = "%j"))
hist(co_prec$doy)

co_locs.sp = co_prec[!duplicated(co_prec$STATION),c("STATION","LATITUDE","LONGITUDE")]
coordinates(co_locs.sp) = ~ LONGITUDE + LATITUDE
proj4string(co_locs.sp) = CRS("+proj=longlat +datum=WGS84")
## Covariate maps:
co_grids = readRDS("data/st_prec/boulder_grids.rds")
co_grids = as(co_grids, "SpatialPixelsDataFrame")
co_locs.sp = spTransform(co_locs.sp, co_grids@proj4string)
plot(raster(co_grids[1]))
points(co_locs.sp, pch="+")
## overlay:
sel.co <- over(co_locs.sp, co_grids[1])
co_locs.sp <- co_locs.sp[!is.na(sel.co$elev_1km),]
T.lst = paste0("2016-02-0", c(1:6))
for(i in 1:length(T.lst)){
  co_locs.sp$x = plyr::join(co_locs.sp@data, co_prec[co_prec$DATE==T.lst[i], c("STATION","PRCP")])$PRCP
  #bubble(co_locs.sp[!is.na(co_locs.sp$x),"x"])
  writeOGR(co_locs.sp, paste0("results/st_prec/co_prec_", T.lst[i], ".shp"), paste0("co_prec_", T.lst[i]), "ESRI Shapefile")
}
## Geographic distances only:
grid.distP <- GSIF::buffer.dist(co_locs.sp["STATION"], co_grids[1], as.factor(1:nrow(co_locs.sp)))
#plot(stack(grid.distP[1:3]))
dnP <- paste(names(grid.distP), collapse="+")
## spacetime hybrid RF model:
fmP <- as.formula(paste("PRCP ~ cdate + doy + elev_1km + PRISM_prec +", dnP))
#fmP0 <- as.formula(paste("PRCP ~ cdate + doy + elev_1km + PRISM_prec"))
ov.prec <- do.call(cbind, list(co_locs.sp@data, over(co_locs.sp, grid.distP), over(co_locs.sp, co_grids[c("elev_1km","PRISM_prec")])))
rm.prec <- plyr::join(co_prec, ov.prec)
rm.prec <- rm.prec[complete.cases(rm.prec[,c("PRCP","elev_1km","cdate")]),]
## 'data.frame':	157870 obs. of 246 variables
m1.prec <- quantregRanger(fmP, rm.prec[sample.int(size = 1e4, n = nrow(rm.prec)),], list(importance = "impurity", num.trees = 150, mtry = 180))
#m1.prec <- quantregRanger(fmP, rm.prec, list(importance = "impurity", num.trees = 150, mtry = 180))
## mtry needs to be set HIGH --- the higher the better, but this increases computational intensity
## For mtry < 20 R-square is close to 0!!
## TAKES 10 minutes on 8-cores Lenovo Legion
## mtry = 160 takes only 6 minutes
m1.prec
# OOB prediction error (MSE):       0.0052395 
# R squared (OOB):                  0.8511794 
sd(co_prec$PRCP, na.rm = TRUE)
xlP.g <- as.list(m1.prec$variable.importance)
print(t(data.frame(xlP.g[order(unlist(xlP.g), decreasing=TRUE)[1:10]])))
# cdate      2400.80112
# doy        1891.66278
# PRISM_prec   93.48800
# elev_1km     72.80193
# layer.89     17.42405
# layer.221    14.08694
# layer.27     14.07233
# layer.139    13.66933
# layer.219    11.99654
# layer.96     11.34862

## Predict some dates (one week):
T.lst = paste0("2016-02-0", c(1:6))
for(i in 1:length(T.lst)){
  if(!file.exists(paste0("results/st_prec/co_PRCP_", T.lst[i], ".tif"))){
    newdata <- cbind(grid.distP@data, co_grids[c("elev_1km","PRISM_prec")]@data)
    newdata$cdate <- floor(unclass(as.POSIXct(as.POSIXct(T.lst[i], format="%Y-%m-%d")))/86400)
    newdata$doy <- as.integer(strftime(as.POSIXct(T.lst[i], format="%Y-%m-%d"), format = "%j"))
    prec.rfT1 <- predict(m1.prec, newdata, quantiles, all = FALSE) 
    ## TH: to speed up use only subset of observations
    ## TH: takes a lot of RAM (could be further optimized)
    co_grids@data[,paste0("PRCP_", T.lst[i])] = ifelse(prec.rfT1[,2]<0, 0, prec.rfT1[,2])*100
    co_grids@data[,paste0("PRCP_", T.lst[i], "_pe")] = (prec.rfT1[,3]-prec.rfT1[,1])/2*100
    #prec.rfT1 <- predict(m1.prec, newdata)
    #co_grids@data[,paste0("PRCP_", T.lst[i])] = ifelse(prec.rfT1$predictions<0, 0, prec.rfT1$predictions)
    #plot(co_grids[paste0("PRCP_", T.lst[i])])
    #points(co_locs.sp, pch="+")
    writeGDAL(co_grids[paste0("PRCP_", T.lst[i])], fname=paste0("results/st_prec/co_PRCP_", T.lst[i], ".tif"), options="COMPRESS=DEFLATE", type = "Int16", mvFlag = "-32768")
    writeGDAL(co_grids[paste0("PRCP_", T.lst[i], "_pe")], fname=paste0("results/st_prec/co_PRCP_", T.lst[i], "_pe.tif"), options="COMPRESS=DEFLATE", type = "Int16", mvFlag = "-32768")
  }
}
