## Comparison RF vs kriging, see slides at: https://github.com/ISRICWorldSoil/GSIF_tutorials/blob/master/geul/5H_Hengl.pdf
## By: tom.hengl@gmail.com, contributions by: Madlene Nussbaum <madlene.nussbaum@bfh.ch> and Marvin Wright <marv@wrig.de>
## Cite as: Hengl et al., "Random Forest as a Generic Framework for Predictive Modeling of Spatial and Spatio-temporal Variables", https://peerj.com/preprints/26693/
## Licence: GNU GPL

list.of.packages <- c("plyr", "parallel", "randomForest", "quantregForest", "plotKML", "GSIF", "RCurl", "raster", "rgdal", "geoR", "gstat", "scales", "gdistance", "entropy", "lattice", "gridExtra", "intamap", "maxlike", "spatstat", "DescTools", "gdistance")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = TRUE)

setwd("/data/git/GeoMLA/RF_vs_kriging")
load(".RData")
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
library(DescTools)

## RANGER connected packages (best install from github):
#devtools::install_github("imbs-hl/ranger")
library(ranger)
## http://philipppro.github.io/Tuning_random_forest/
#devtools::install_github("mlr-org/mlrMBO")
#devtools::install_github("PhilippPro/tuneRF")
library(tuneRanger)
## Load all functions prepared for this exerice:
source('./R/RFsp_functions.R')

## ** Example Meuse data set ----------------------------------------------
demo(meuse, echo=FALSE)

## Zinc predicted using ordinary kriging ----
zinc.geo <- as.geodata(meuse["zinc"])
ini.v <- c(var(log1p(zinc.geo$data)),500)
zinc.vgm <- likfit(zinc.geo, lambda=0, ini=ini.v, cov.model="exponential")
locs = meuse.grid@coords
zinc.ok <- krige.conv(zinc.geo, locations=locs, krige=krige.control(obj.model=zinc.vgm))
meuse.grid$zinc_ok = zinc.ok$predict
## Prediction error 1 s.d.:
meuse.grid$zinc_ok_range = sqrt(zinc.ok$krige.var)

## Zinc predicted using RF and buffer distances only ----
grid.dist0 <- GSIF::buffer.dist(meuse["zinc"], meuse.grid[1], as.factor(1:nrow(meuse)))
dn0 <- paste(names(grid.dist0), collapse="+")
fm0 <- as.formula(paste("zinc ~ ", dn0))
ov.zinc <- over(meuse["zinc"], grid.dist0)
## Get a good estimate of "mtry" (fine-tuning):
rm.zinc <- cbind(meuse@data["zinc"], ov.zinc)
rt.zinc <- makeRegrTask(data = rm.zinc, target = "zinc", check.data=FALSE)
# Approximated time for tuning: 5M 25S
# Make reproducible tuning
set.seed(1)
t.zinc <- tuneRanger(rt.zinc, num.trees = 150, build.final.model = FALSE, parameters = list(replace = FALSE))
t.zinc
## With seed = 1. 
# Recommended parameter settings: 
#   mtry min.node.size sample.fraction
# 1   76             4         0.70777
# Results: 
#   mse exec.time
# 1 59042.89    0.1022
pars.zinc = list(mtry=t.zinc$recommended.pars$mtry, min.node.size=t.zinc$recommended.pars$min.node.size, sample.fraction=t.zinc$recommended.pars$sample.fraction, num.trees=150, seed=1)
m.zinc <- ranger(fm0, rm.zinc, quantreg=TRUE, mtry=t.zinc$recommended.pars$mtry, min.node.size=t.zinc$recommended.pars$min.node.size, sample.fraction=t.zinc$recommended.pars$sample.fraction, num.trees=150, seed=1)
m.zinc
# R squared (OOB): 0.5126192
## fit model without QRF
m.zinc0 <- ranger(fm0, rm.zinc, mtry=t.zinc$recommended.pars$mtry, min.node.size=t.zinc$recommended.pars$min.node.size, sample.fraction=t.zinc$recommended.pars$sample.fraction, num.trees=150, seed=1)
## predictions
zinc.rfd <- predict(m.zinc, grid.dist0@data, type="quantiles", quantiles=quantiles)$predictions
meuse.grid$zinc_rfd = zinc.rfd[,2]
## This predicts the "median" value; to predict "mean" using ranger without quantreg=TRUE
meuse.grid$zinc_rfd0 = predict(m.zinc0, grid.dist0@data)$predictions
summary(meuse.grid$zinc_rfd0); summary(meuse.grid$zinc_rfd) 
hexbin::hexbinplot(meuse.grid$zinc_rfd0 ~ meuse.grid$zinc_rfd)
## TH: predictions of the "median" values are somewhat smaller than predictions of the "mean" values,
## Prediction error s.d.:
meuse.grid$zinc_rfd_range = (zinc.rfd[,3]-zinc.rfd[,1])/2

## RF and coordinates as covariates only (to be fair, same tuning) ----
fmc <- zinc ~ x + y
rm.zinc.coord <- cbind(meuse@data["zinc"], meuse@coords)
rt.zinc.coord <- makeRegrTask(data = rm.zinc.coord, target = "zinc")
set.seed(1)
t.zinc.coord <- tuneRanger(rt.zinc.coord, num.trees = 150, build.final.model = FALSE, parameters = list(replace = FALSE))
pars.zinc.coord = list(mtry= t.zinc.coord$recommended.pars$mtry, min.node.size=t.zinc.coord$recommended.pars$min.node.size, sample.fraction=t.zinc.coord$recommended.pars$sample.fraction, num.trees=150, seed = 1)
m.zinc.coord <- ranger(fmc, rm.zinc.coord, quantreg=TRUE, mtry=t.zinc.coord$recommended.pars$mtry, min.node.size=t.zinc.coord$recommended.pars$min.node.size, sample.fraction=t.zinc.coord$recommended.pars$sample.fraction, num.trees=150, seed=1)
m.zinc.coord
## R-squared still quite high 0.517
zinc.rfc <- predict(m.zinc.coord, meuse.grid@coords, type="quantiles", quantiles = quantiles)$predictions
meuse.grid$zinc_rfc = zinc.rfc[,2]
meuse.grid$zinc_rfc_range = (zinc.rfc[,3]-zinc.rfc[,1])/2

## Plot predictions UK and RF next to each other:
pdf(file = "results/meuse/Fig_comparison_OK_RF_zinc_meuse.pdf", width=9.5, height=8)
var.max = max(c(meuse.grid$zinc_rfc_range, meuse.grid$zinc_rfd_range, meuse.grid$zinc_ok_range))
var.min = min(c(meuse.grid$zinc_rfc_range, meuse.grid$zinc_rfd_range, meuse.grid$zinc_ok_range))
par(mfrow=c(2,3), oma=c(0,0,0,1), mar=c(0,0,4,3))
plot(log1p(raster(meuse.grid["zinc_ok"])), col=leg, zlim=c(4.7,7.6), main="Ordinary Kriging (OK)", axes=FALSE, box=FALSE, axis.args=axis.ls)
points(meuse, pch="+")
plot(log1p(raster(meuse.grid["zinc_rfd"])), col=leg, zlim=c(4.7,7.6), main="Random Forest (RF), buffers", axes=FALSE, box=FALSE, axis.args=axis.ls)
points(meuse, pch="+")
plot(log1p(raster(meuse.grid["zinc_rfc"])), col=leg, zlim=c(4.7,7.6), main="Random Forest (RF), coordinates", axes=FALSE, box=FALSE, axis.args=axis.ls)
points(meuse, pch="+")
plot(raster(meuse.grid["zinc_ok_range"]), col=rev(bpy.colors()), zlim=c(var.min,var.max), main="OK prediction range", axes=FALSE, box=FALSE)
points(meuse, pch="+")
plot(raster(meuse.grid["zinc_rfd_range"]), col=rev(bpy.colors()), zlim=c(var.min,var.max), main="RF prediction range, buffers", axes=FALSE, box=FALSE)
points(meuse, pch="+")
plot(raster(meuse.grid["zinc_rfc_range"]), col=rev(bpy.colors()), zlim=c(var.min,var.max), main="RF prediction range, coordinates", axes=FALSE, box=FALSE)
points(meuse, pch="+")
dev.off()
## TH: RF creates smoother predictions than OK

## Cross-validation Meuse data set geographical covs only ----
ss <- seq(0,1,0.01) # define probabilites for coverage plot
cv.RF = cv_numeric(varn="zinc", points=meuse, covs=meuse.grid, cpus=1, method="ranger", OK=TRUE, spcT=FALSE, pars.ranger=pars.zinc, predDist = ss)
cv.OK = cv_numeric(varn="zinc", points=meuse, covs=meuse.grid, cpus=1, method="geoR", OK=TRUE, spcT=FALSE)
cv.RF$Summary
cv.OK$Summary
cv.RF$Summary$CCC_est^2; cv.OK$Summary$CCC_est^2
## Compare with the standard error of the mean:
sqrt(var(meuse$zinc)/nrow(meuse))
## Plot residuals vgm:
x = plyr::join_all(list(data.frame(r1=cv.RF$CV_residuals$Observed-cv.RF$CV_residuals$Predicted, id=cv.RF$CV_residuals$SOURCEID), data.frame(r2=cv.OK$CV_residuals$Observed-cv.OK$CV_residuals$Predicted, id=cv.OK$CV_residuals$SOURCEID)))
x = plyr::join(x, data.frame(id=row.names(meuse@data), zinc=meuse$zinc, x=meuse@coords[,1], y=meuse@coords[,2]))
plot_vgm(zinc~1, x, meuse.grid, r1="r1", r2="r2", main="Zinc (Meuse)")

## plot results
lim.zinc = range(meuse$zinc, na.rm = TRUE)
pdf(file = "results/meuse/Fig_correlation_plots_OK_RF_zinc_meuse.pdf", width=9, height=5)
par(oma=c(0,0,0,1), mar=c(0,0,0,2))
pfun.loess <- function(x,y, ...){ 
  panel.xyplot(x,y, ...)
  panel.abline(0,1,lty=2,lw=1,col="grey60")
  panel.loess(x,y,span=0.5, col = "grey30", lwd = 1.7)
}
plt.RF = xyplot(cv.RF[[1]]$Observed~cv.RF[[1]]$Predicted, asp=1, par.settings=list(plot.symbol = list(col=alpha("black", 0.6), fill=alpha("red", 0.6), pch=21, cex=0.9)), scales=list(x=list(log=TRUE, equispaced.log=FALSE), y=list(log=TRUE, equispaced.log=FALSE)), ylab="observed zinc [mg/kg]", xlab="zinc predicted by RF [mg/kg]", panel = pfun.loess, xlim=lim.zinc, ylim=lim.zinc)
plt.OK = xyplot(cv.OK[[1]]$Observed~cv.OK[[1]]$Predicted, asp=1, par.settings=list(plot.symbol = list(col=alpha("black", 0.6), fill=alpha("red", 0.6), pch=21, cex=0.9)), scales=list(x=list(log=TRUE, equispaced.log=FALSE), y=list(log=TRUE, equispaced.log=FALSE)), ylab="observed zinc [mg/kg]", xlab="zinc predicted by geoR [mg/kg]", panel = pfun.loess, xlim=lim.zinc, ylim=lim.zinc)
grid.arrange(plt.OK, plt.RF, ncol=2)
dev.off()

## Quantile coverage plot -------------
cvRF.one.side.coverage <- sapply( ss, FUN = function(xx){ 
  sum( cv.RF$CV_residuals$Observed <= cv.RF$CV_residuals[ grepl(paste0("Pred.Quantile.", xx, "$" ), names(cv.RF$CV_residuals))] ) / nrow( cv.RF$CV_residuals ) } )
cvOK.one.side.coverage <- sapply( ss, FUN = function(xx){ 
  sum( cv.OK$CV_residuals$Observed <= cv.OK$CV_residuals$Predicted + cv.OK$CV_residuals$sdPE * qnorm(xx) ) / nrow( cv.RF$CV_residuals ) })

pdf("results/meuse/Fig_coverage-probabilties_meuse.pdf", width = 9, height = 4.4)
par( oma = c(0,0,0,0), ps = 10, mar=c(2.7, 2.7, 1, 0.8), mfrow = c(1, 2), mgp = c(1.5, 0.3, 0))
## OK
plot(x = ss, y = cvOK.one.side.coverage, type = "l", pch = 20, ylab = "coverage probabilities OK", xlab="nominal probabilities", asp = 1, xaxt = "n", yaxt = "n", lwd = 1.3, col="blue")
abline(0,1, lty = 2, col = "grey60")
axis( 1, tck = 0.02, las = 1, lwd = 0.95, mgp = c(0.5, 0.1, 0))
axis( 2, tck = 0.02, las = 1, lwd = 0.95 )
abline(v = 0.05, lty = "dotted", col = "grey20")
abline(v = 0.95, lty = "dotted", col = "grey20")
## RF
plot(x = ss, y = cvRF.one.side.coverage, type = "l", pch = 20, ylab = "coverage probabilities RFsp", xlab="nominal probabilities", asp = 1, xaxt = "n", yaxt = "n", lwd = 1.3, col="red")
abline(0,1, lty = 2, col = "grey60")
axis( 1, tck = 0.02, las = 1, lwd = 0.95, mgp = c(0.5, 0.1, 0))
axis( 2, tck = 0.02, las = 1, lwd = 0.95 )
abline(v = 0.05, lty = "dotted", col = "grey20")
abline(v = 0.95, lty = "dotted", col = "grey20")
dev.off()

## RF with combined covariates ----
meuse.grid$SW_occurrence = readGDAL("data/meuse/Meuse_GlobalSurfaceWater_occurrence.tif")$band1[meuse.grid@grid.index]
meuse.grid$AHN = readGDAL("data/meuse/ahn.asc")$band1[meuse.grid@grid.index]
meuse.grid$LGN5 = as.factor(readGDAL("data/meuse/lgn5.asc")$band1[meuse.grid@grid.index])
## convert to numeric via PCA
grids.spc = spc(meuse.grid, as.formula("~ SW_occurrence + AHN + ffreq + dist"))
## fit hybrid RF model:
fm1 <- as.formula(paste("zinc ~ ", dn0, " + ", paste(names(grids.spc@predicted), collapse = "+")))
ov.zinc1 <- over(meuse["zinc"], grids.spc@predicted)
rm.zinc1 <- do.call(cbind, list(meuse@data["zinc"], ov.zinc, ov.zinc1))
rt.zinc1 <- makeRegrTask(data = rm.zinc1, target = "zinc")
set.seed(1)
t.zinc1 <- tuneRanger(rt.zinc1, num.trees = 500, build.final.model = FALSE, parameters = list(replace = FALSE))
m1.zinc <- ranger(fm1, rm.zinc1, mtry=t.zinc1$recommended.pars$mtry, min.node.size=t.zinc1$recommended.pars$min.node.size, num.trees=500, importance = "impurity", seed=1)
m1.zinc
## R squared (OOB): 0.6412438
meuse.grid$zinc_rfd1 = predict(m1.zinc, cbind(grid.dist0@data, grids.spc@predicted@data))$predictions
xl <- as.list(ranger::importance(m1.zinc))

pdf(file = "results/meuse/Fig_RF_covs_covar-importance.pdf", width=5, height=5.5)
par(mfrow=c(1,1),oma=c(0.7,2,0,1), mar=c(4,3.5,1,0))
plot(vv <- t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[10:1]])), 1:10, type = "n", ylab = "", yaxt = "n", xlab = "Variable Importance (Node Impurity)")
abline(h = 1:10, lty = "dotted", col = "grey60")
points(vv, 1:10)
axis(2, 1:10, labels = dimnames(vv)[[1]], las = 2)
dev.off()

## RF without geographical coordinates
rm.zinc2 <- cbind(meuse@data["zinc"], ov.zinc1)
rt.zinc2 <- makeRegrTask(data = rm.zinc2, target = "zinc")
set.seed(1)
t.zinc2 <- tuneRanger(rt.zinc2, num.trees = 150, build.final.model = FALSE, parameters = list(replace = FALSE))
m2.zinc <- ranger(paste("zinc ~ ", paste(names(grids.spc@predicted), collapse = "+")), rm.zinc2, mtry=t.zinc2$recommended.pars$mtry, min.node.size=t.zinc2$recommended.pars$min.node.size, num.trees=150, importance = "impurity", seed=1)
m2.zinc
## R squared (OOB): 0.5751939
meuse.grid$zinc_rfd2 = predict(m2.zinc, grids.spc@predicted@data)$predictions

pdf(file = "results/meuse/Fig_RF_covs_bufferdist_zinc_meuse.pdf", width=9, height=5)
par(mfrow=c(1,2), oma=c(0,0,0,1), mar=c(0,0,4,3))
plot(log1p(raster(meuse.grid["zinc_rfd2"])), col=leg, zlim=c(4.7,7.6), main="Random Forest (RF) covs only", axes=FALSE, box=FALSE, axis.args=axis.ls)
points(meuse, pch="+")
plot(log1p(raster(meuse.grid["zinc_rfd1"])), col=leg, zlim=c(4.7,7.6), main="Random Forest (RF) covs + buffer dist.", axes=FALSE, box=FALSE, axis.args=axis.ls)
points(meuse, pch="+")
dev.off()


## ** SIC 1997 data set ---------------------------------------------------
## measurements made in Switzerland on the 8th of May 1986
sic97.sp = readRDS("data/rainfall/sic97.rds")
swiss1km = readRDS("data/rainfall/swiss1km.rds")
ov2 = over(y=swiss1km, x=sic97.sp)
sel.d = which(!is.na(ov2$DEM))
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
#swiss1km$rainfall_UK_range = sd.range(rain.uk$predict, sqrt(rain.uk$krige.var))/2
swiss1km$rainfall_UK_range = sqrt(rain.uk$krige.var)

## SIC 1997 Random Forest example ----
swiss.dist0 <- GSIF::buffer.dist(sic97.sp["rainfall"], swiss1km[1], as.factor(1:nrow(sic97.sp))) ## takes 2-3 mins!
## Deriving buffer distances is computationally very intensive
ov.swiss = over(sic97.sp["rainfall"], swiss.dist0)
sw.dn0 <- paste(names(swiss.dist0), collapse="+")
sw.fm1 <- as.formula(paste("rainfall ~ ", sw.dn0, " + CHELSA_rainfall + DEM"))
ov.rain <- over(sic97.sp["rainfall"], swiss1km[1:2])
sw.rm = do.call(cbind, list(sic97.sp@data["rainfall"], ov.rain, ov.swiss))
## fine-tune RF:
# rt.rain <- makeRegrTask(data = sw.rm[complete.cases(sw.rm[,all.vars(sw.fm1)]),], target = "rainfall")
# estimateTimeTuneRanger(rt.rain, num.threads=88)
## Too time-consuming >> do on number cruncher
# set.seed(1)
# t.rain <- tuneRanger(rt.rain, num.trees = 150, build.final.model = FALSE, num.threads = 88, parameters = list(replace = FALSE))
# t.rain
# pars.rain = list(mtry= t.rain$recommended.pars$mtry, min.node.size=t.rain$recommended.pars$min.node.size, sample.fraction=t.rain$recommended.pars$sample.fraction, num.trees=150, importance = "impurity", seed=1)
pars.rain = list(mtry=27, min.node.size=2, sample.fraction=0.9930754, num.trees=150, importance = "impurity", seed=1)
m1.rain <- ranger(sw.fm1, sw.rm[complete.cases(sw.rm),], mtry=27, min.node.size=2, sample.fraction=0.9930754, num.trees=150, importance = "impurity", seed=1, quantreg=TRUE)
m1.rain
## R squared (OOB): 0.8311812
rain.rfd1 <- predict(m1.rain, cbind(swiss.dist0@data, swiss1km@data), type="quantiles", quantiles = quantiles)$predictions
## now more computational...
swiss1km$rainfall_rfd1 = rain.rfd1[,2]
swiss1km$rainfall_rfd1_var = (rain.rfd1[,3]-rain.rfd1[,1])/2

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
plot(sqrt(raster(swiss1km["rainfall_UK_range"])), col=rev(bpy.colors()), main="Universal kriging (UK) prediction error", axes=FALSE, box=FALSE, zlim=c(0,rainv.max))
points(sic97.sp, pch="+")
plot(raster(swiss1km["rainfall_rfd1_var"]), col=rev(bpy.colors()), main="Random Forest (RF) prediction error", axes=FALSE, box=FALSE, zlim=c(0,rainv.max))
points(sic97.sp, pch="+")
dev.off()

## Cross-validation SIC 1997 ----
## (computationally intensive - takes few minutes because buffer distances have to be re-computed!)
cv.RF2 = cv_numeric(varn="rainfall", points=sic97.sp, covs=swiss1km[c("CHELSA_rainfall","DEM")], cpus=1, method="ranger", spcT=FALSE, pars.ranger = pars.rain[-c(5)])
cv.UK = cv_numeric(varn="rainfall", points=sic97.sp, covs=swiss1km[c("CHELSA_rainfall","DEM")], cpus=1, method="geoR", spcT=FALSE)
cv.RF2$Summary
cv.UK$Summary
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
plt.RF2 = xyplot(cv.RF2[[1]]$Observed~cv.RF2[[1]]$Predicted, asp=1, par.settings=list(plot.symbol = list(col=alpha("black", 0.6), fill=alpha("red", 0.6), pch=21, cex=0.9)), scales=list(x=list(log=TRUE, equispaced.log=FALSE), y=list(log=TRUE, equispaced.log=FALSE)), ylab="observed precipitation [mm]", xlab="precipitation predicted by RF [mm]", panel = pfun.loess, xlim=lim.rain, ylim=lim.rain)
plt.UK = xyplot(cv.UK[[1]]$Observed~cv.UK[[1]]$Predicted, asp=1, par.settings=list(plot.symbol = list(col=alpha("black", 0.6), fill=alpha("red", 0.6), pch=21, cex=0.9)), scales=list(x=list(log=TRUE, equispaced.log=FALSE), y=list(log=TRUE, equispaced.log=FALSE)), ylab="observed precipitation [mm]", xlab="precipitation predicted by geoR [mm]", panel = pfun.loess, xlim=lim.rain, ylim=lim.rain)
grid.arrange(plt.UK, plt.RF2, ncol=2)
dev.off()

## ** Ebergotzen binomial variable --------------------------------------
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
m1.Parabraunerde <- ranger(eb.fm1, rm.eberg2[complete.cases(rm.eberg2),], importance = "impurity", probability = TRUE, seed = 1)
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
par(oma=c(0,0,0,0.5), mar=c(0,0,3.5,1))
plot(raster(eberg_grid["Parabraunerde_TRUE"]), col=SAGA_pal[["SG_COLORS_YELLOW_RED"]], zlim=c(0,1), main="Parabraunerde class (RF)", axes=FALSE, box=FALSE)
points(eberg[eberg$Parabraunerde=="TRUE"&!is.na(eberg$Parabraunerde)&sel.eberg,], pch=19, cex=.5)
points(eberg[eberg$Parabraunerde=="FALSE"&!is.na(eberg$Parabraunerde)&sel.eberg,], pch="+", cex=.8)
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


## ** NRCS: Weighted regression RF -------------------------------------
carson <- read.csv(file="data/NRCS/carson_CLYPPT.csv")
summary(carson$CLYPPT)
carson$DEPTH.f = ifelse(is.na(carson$DEPTH), 20, carson$DEPTH)
carson1km <- readRDS("data/NRCS/carson_covs1km.rds")
coordinates(carson) <- ~X+Y
proj4string(carson) = carson1km@proj4string
rm.carson <- cbind(as.data.frame(carson), over(carson["CLYPPT"], carson1km))
fm.clay <- as.formula(paste("CLYPPT ~ DEPTH.f + ", paste(names(carson1km), collapse = "+")))
rm.carson <- rm.carson[complete.cases(rm.carson[,all.vars(fm.clay)]),]
rm.carson.s <- rm.carson[sample.int(size=2000, nrow(rm.carson)),]
m.clay <- ranger(fm.clay, rm.carson.s, num.trees=150, mtry=25, case.weights=1/(rm.carson.s$CLYPPT.sd^2), quantreg = TRUE)
m.clay
carson1km$DEPTH.f = 30
clay.rfd <- predict(m.clay, carson1km@data, type="quantiles", quantiles=quantiles)$predictions
carson1km$clay_rfd = ifelse(clay.rfd[,2]<10, 10, ifelse(clay.rfd[,2]>35, 35, clay.rfd[,2]))
summary(carson1km$clay_rfd)
carson1km$clay_rfd_var = (clay.rfd[,3]-clay.rfd[,1])/2
m.clay0 <- ranger(fm.clay, rm.carson.s, num.trees=150, mtry=25, seed = 1, quantreg = TRUE)
m.clay0
clay0.rfd <- predict(m.clay0, carson1km@data, type="quantiles", quantiles=quantiles)$predictions
carson1km$clay0_rfd = ifelse(clay0.rfd[,2]<10, 10, ifelse(clay0.rfd[,2]>35, 35, clay0.rfd[,2]))
carson1km$clay0_rfd_var = (clay0.rfd[,3]-clay0.rfd[,1])/2

pdf(file = "results/NRCS/Fig_clay_RF_weighted.pdf", width=10, height=9)
par(mfrow=c(2,2), oma=c(0,0,0,0.5), mar=c(0,0,1.5,1))
par(oma=c(0,0,0,0.5), mar=c(0,0,3.5,1))
plot(raster(carson1km["clay_rfd"]), col=leg, main="RF predictions clay (weigthed)", axes=FALSE, box=FALSE, zlim=c(10,35))
points(rm.carson.s$X, rm.carson.s$Y, pch="+", cex=.8)
plot(raster(carson1km["clay0_rfd"]), col=leg, main="RF predictions clay", axes=FALSE, box=FALSE, zlim=c(10,35))
points(rm.carson.s$X, rm.carson.s$Y, pch="+", cex=.8)
plot(raster(carson1km["clay_rfd_var"]), col=rev(bpy.colors()), main="Prediction range (RF weigthed)", axes=FALSE, box=FALSE, zlim=c(0,28))
points(rm.carson.s$X, rm.carson.s$Y, pch="+", cex=.8)
plot(raster(carson1km["clay0_rfd_var"]), col=rev(bpy.colors()), main="Prediction range (RF)", axes=FALSE, box=FALSE, zlim=c(0,28))
points(rm.carson.s$X, rm.carson.s$Y, pch="+", cex=.8)
dev.off()


## ** Geochemical USA - Multivariate case -----------------------------------
## Geochemicals USA (https://mrdata.usgs.gov/geochem/)
geochem = readRDS("data/geochem/geochem.rds")
## negative values are in fact detection limits:
for(i in c("PB_ICP40","CU_ICP40","K_ICP40","MG_ICP40")) { geochem[,i] = ifelse(geochem[,i] < 0, abs(geochem[,i])/2, geochem[,i])  }
coordinates(geochem) = ~coords.x1 + coords.x2
proj4string(geochem) = "+proj=longlat +ellps=clrk66 +towgs84=-9.0,151.0,185.0,0.0,0.0,0.0,0.0 +no_defs"
geochem$TYPEDESC = as.factor(paste(geochem$TYPEDESC))
summary(geochem$TYPEDESC)
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
type.mat = data.frame(model.matrix(~TYPE-1, rm.geochem))
typed.mat = data.frame(model.matrix(~TYPEDESC-1, rm.geochem))
## Final regression matrix
rm.geochem.e = do.call(cbind, list(rm.geochem[,c("Y",paste0("PC",1:21))], type.mat, typed.mat))
fm.g = as.formula(paste0("Y ~ ", paste0(names(rm.geochem.e)[-1], collapse = "+")))
fm.g
m1.geochem <- ranger::ranger(fm.g, rm.geochem.e[complete.cases(rm.geochem.e),], importance = "impurity", seed = 1)
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


## ** Spatiotemporal prediction using Random Forest ------------------------
## Daily precipitation obtained from https://www.ncdc.noaa.gov/cdo-web/search
## compare with: https://www.r-bloggers.com/part-4a-modelling-predicting-the-amount-of-rain/
## NCDC data access explained at: http://neondataskills.org/R/COOP-precip-data-R
co_prec = readRDS("data/st_prec/boulder_prcp.rds")
str(co_prec)
# 'data.frame':	176467 obs. of  16 variables:
summary(co_prec$PRCP)
## Derive cumulative day:
co_prec$cdate = floor(unclass(as.POSIXct(as.POSIXct(paste(co_prec$DATE), format="%Y-%m-%d")))/86400)
## Day of the year:
co_prec$doy = as.integer(strftime(as.POSIXct(paste(co_prec$DATE), format="%Y-%m-%d"), format = "%j"))

co_locs.sp = co_prec[!duplicated(co_prec$STATION),c("STATION","LATITUDE","LONGITUDE")]
coordinates(co_locs.sp) = ~ LONGITUDE + LATITUDE
proj4string(co_locs.sp) = CRS("+proj=longlat +datum=WGS84")
## Covariate maps:
co_grids = readRDS("data/st_prec/boulder_grids.rds")
co_grids = as(co_grids, "SpatialPixelsDataFrame")
co_locs.sp = spTransform(co_locs.sp, co_grids@proj4string)
## overlay:
sel.co <- over(co_locs.sp, co_grids[1])
co_locs.sp <- co_locs.sp[!is.na(sel.co$elev_1km),]
T.lst = paste0("2016-02-0", c(1:6))
for(i in 1:length(T.lst)){
  co_locs.sp$x = plyr::join(co_locs.sp@data, co_prec[co_prec$DATE==T.lst[i], c("STATION","PRCP")])$PRCP
  writeOGR(co_locs.sp, paste0("results/st_prec/co_prec_", T.lst[i], ".shp"), paste0("co_prec_", T.lst[i]), "ESRI Shapefile")
}
## Geographic distances only:
grid.distP <- GSIF::buffer.dist(co_locs.sp["STATION"], co_grids[1], as.factor(1:nrow(co_locs.sp)))
dnP <- paste(names(grid.distP), collapse="+")
## spacetime 'hybrid' RF model:
fmP <- as.formula(paste("PRCP ~ cdate + doy + elev_1km + PRISM_prec +", dnP))
#fmP0 <- as.formula(paste("PRCP ~ cdate + doy + elev_1km + PRISM_prec"))
ov.prec <- do.call(cbind, list(co_locs.sp@data, over(co_locs.sp, grid.distP), over(co_locs.sp, co_grids[c("elev_1km","PRISM_prec")])))
rm.prec <- plyr::join(co_prec, ov.prec)
rm.prec <- rm.prec[complete.cases(rm.prec[,c("PRCP","elev_1km","cdate")]),]
## 'data.frame':	157870 obs. of 246 variables
# rt.prec <- makeRegrTask(data = rm.prec, target = "PRCP")
# estimateTimeTuneRanger(rt.prec, num.threads = 88)
# # Time consuming >> do on number cruncher
# set.seed(1)
# t.prec <- tuneRanger(rt.prec, num.trees = 150, build.final.model = FALSE)
# pars.prec = list(mtry= t.prec$recommended.pars$mtry, min.node.size=t.prec$recommended.pars$min.node.size, sample.fraction=t.prec$recommended.pars$sample.fraction, num.trees=150, importance = "impurity", seed = 1)
pars.prec = list(mtry=212, min.node.size= 2, sample.fraction=0.9553763, num.trees=150, seed=1)
m1.prec <- ranger(formula = fmP, data = rm.prec, mtry=212, min.node.size= 2, sample.fraction=0.9553763, num.trees=150, seed=1, quantreg=TRUE, importance='impurity')
m1.prec
# Type:                             Regression 
# Number of trees:                  150 
# Sample size:                      157870 
# Number of independent variables:  229 
# Mtry:                             212 
# Target node size:                 2 
# Variable importance mode:         impurity 
# OOB prediction error (MSE):       0.005209039 
# R squared (OOB):                  0.8520446 
sd(co_prec$PRCP, na.rm = TRUE)
# [1] 0.1864649
xlP.g <- as.list(ranger::importance(m1.prec))
print(t(data.frame(xlP.g[order(unlist(xlP.g), decreasing=TRUE)[1:10]])))
# cdate      2333.27929
# doy        1807.15683
# PRISM_prec   86.28802
# elev_1km     68.71017
# layer.89     17.05671
# layer.221    13.69766
# layer.139    12.62467
# layer.27     11.92333
# layer.35     11.63324
# layer.219    11.36218

## Predict some dates (one week):
T.lst = paste0("2016-02-0", c(1:6))
for(i in 1:length(T.lst)){
  if(!file.exists(paste0("results/st_prec/co_PRCP_", T.lst[i], ".tif"))){
    newdata <- cbind(grid.distP@data, co_grids[c("elev_1km","PRISM_prec")]@data)
    newdata$cdate <- floor(unclass(as.POSIXct(as.POSIXct(T.lst[i], format="%Y-%m-%d")))/86400)
    newdata$doy <- as.integer(strftime(as.POSIXct(T.lst[i], format="%Y-%m-%d"), format = "%j"))
    prec.rfT1 <- predict(m1.prec, newdata, type="quantiles", quantiles=quantiles)$predictions
    co_grids@data[,paste0("PRCP_", T.lst[i])] = ifelse(prec.rfT1[,2]<0, 0, prec.rfT1[,2])*100
    co_grids@data[,paste0("PRCP_", T.lst[i], "_pe")] = (prec.rfT1[,3]-prec.rfT1[,1])/2*100
    #plot(co_grids[paste0("PRCP_", T.lst[i], "_pe")])
    #points(co_locs.sp, pch="+")
    writeGDAL(co_grids[paste0("PRCP_", T.lst[i])], fname=paste0("results/st_prec/co_PRCP_", T.lst[i], ".tif"), options="COMPRESS=DEFLATE", type = "Int16", mvFlag = "-32768")
    writeGDAL(co_grids[paste0("PRCP_", T.lst[i], "_pe")], fname=paste0("results/st_prec/co_PRCP_", T.lst[i], "_pe.tif"), options="COMPRESS=DEFLATE", type = "Int16", mvFlag = "-32768")
  }
}

## Cross-validation with Leave-Station-Out ----
## It is very important to take whole stations out!
cv_st = function(rm.prec, co_loc.sp, fmP, idcol="STATION", nfold=5, pars.ranger){
  sel.ul <- dismo::kfold(levels(rm.prec[,idcol]), k=nfold)
  out = list(NULL)
  for(j in 1:nfold){
    rm.prec.s = rm.prec[which(rm.prec[,idcol] %in% levels(rm.prec[,idcol])[sel.ul==j]),]
    rm.prec.t = rm.prec[which(rm.prec[,idcol] %in% levels(rm.prec[,idcol])[!sel.ul==j]),]
    ## subset geographical distances
    dnP.s = paste("layer.", which(co_locs.sp@data[,idcol] %in% levels(rm.prec[,idcol])[!sel.ul==j]), sep="", collapse="+")
    fmP.s = as.formula(paste("PRCP ~ cdate + doy + elev_1km + PRISM_prec +", dnP.s))
    mT.prec <- ranger(formula = fmP.s, data = rm.prec.t, mtry=length(all.vars(fmP.s))-30, min.node.size= 2, sample.fraction=0.9553763, num.trees=150, seed=1)  
    s.prec = predict(mT.prec, rm.prec.s)
    out[[j]] = data.frame(Observed=rm.prec.s$PRCP, Predicted=s.prec$predictions)
  }
  return(out)
}
cv.PRCP = cv_st(rm.prec, co_loc.sp, fmP, idcol="STATION", nfold=5, pars.ranger)
cv.PRCP = do.call(rbind, cv.PRCP)
## RMSE
sqrt(mean((cv.PRCP$Observed - cv.PRCP$Predicted)^2, na.rm = T))
## 0.0696
## Results comparable to krigeST results
## CCC
DescTools::CCC(cv.PRCP$Observed, cv.PRCP$Predicted, ci = "z-transform", conf.level = 0.95, na.rm=TRUE)$rho.c
## 0.925

## ** Spatiotemporeal kriging ----
library(sp)
library(spacetime)
library(gstat)
co_locs.sp@data$PRISM_prec <- over(co_locs.sp, co_grids[c("elev_1km","PRISM_prec")])$PRISM_prec
data <- co_prec[,c("STATION","DATE","PRCP", "cdate", "doy","ELEVATION"), drop=FALSE]
str(data)
# summary(lm(PRCP ~ ELEVATION + cosDoy + cdate, data))
## No correlation = no grounds to use covariates
locs.st = plyr::join(data[c("STATION","DATE","PRCP")], as.data.frame(co_locs.sp))
sel.st = !is.na(locs.st$LONGITUDE) & !is.na(locs.st$PRCP)
co_time <- as.POSIXct(sel.st$DATE)
stsdf <- STIDF(SpatialPoints(locs.st[sel.st, c("LONGITUDE","LATITUDE")]), time=as.POSIXct(locs.st[sel.st,"DATE"]), data=data[sel.st,,drop=FALSE])
## coerce to more compact st class
stsdf <- as(stsdf, "STSDF")
stsdf@sp@proj4string = co_grids@proj4string
str(stsdf)
## variogram modeling:
empStVgm <- gstat::variogramST(PRCP~1, stsdf, tlags = 0:3)
## takes few minutes due to the data size ...
empStVgm$dist <- empStVgm$dist/1e3
empStVgm$avgDist <- empStVgm$avgDist/1e3
smmFit <- fit.StVariogram(empStVgm, vgmST("sumMetric", 
                                    space=vgm(0.015, "Sph", 60, 0.01),
                                    time=vgm(0.035, "Sph", 60, 0.001),
                                    joint=vgm(0.035, "Sph", 30, 0.001),
                                    stAni=1),
                                    lower=c(0,0.01,0,
                                            0,0.01,0,
                                            0,0.01,0,
                                            0.05),
                                    control=list(parscale=c(1,1e3,1,
                                      1,1e3,1,
                                      1,1e3,1,
                                      1)))
## Rescale:
smmFit$space$range <- smmFit$space$range*1e3
smmFit$joint$range <- smmFit$joint$range*1e3
smmFit$stAni <- smmFit$stAni*1e3
## Predict
T.lst = paste0("2016-02-0", c(1:6))
sel.T = which(stsdf@endTime %in% as.POSIXct(T.lst))
predST <- krigeST(PRCP~1, stsdf[,sel.T], STF(co_grids, time = stsdf@time[sel.T]), smmFit, nmax = 15, computeVar = TRUE)
## takes few minutes ...
## Export predictions:
library(rgdal)
for(i in 1:length(T.lst)){
  writeGDAL(predST[,i, "var1.pred"], fname=paste0("RF_vs_kriging/results/st_prec/krigeSt_PRCP_", T.lst[i], ".tif"), options="COMPRESS=DEFLATE")
  writeGDAL(predST[,i, "var1.var"], fname=paste0("RF_vs_kriging/results/st_prec/krigeSt_PRCP_", T.lst[i], "_pe.tif"), options="COMPRESS=DEFLATE")
}

## Plot the two next to each other:
co_tif = raster::stack(c(paste0("RF_vs_kriging/results/st_prec/krigeSt_PRCP_2016-02-0", 2:5, ".tif"), paste0("RF_vs_kriging/results/st_prec/krigeSt_PRCP_2016-02-0", 2:5, "_pe.tif" ), paste0("RF_vs_kriging/results/st_prec/co_PRCP_2016-02-0", 1:4, ".tif" ), paste0("RF_vs_kriging/results/st_prec/co_PRCP_2016-02-0", 1:4, "_pe.tif")))
#plot(co_tif)
co_tif = as(co_tif, "SpatialGridDataFrame")
## scale back to original scale
for(i in grep(glob2rx("krigeSt_PRCP_*_pe"), names(co_tif))){ co_tif@data[,i] = sqrt(co_tif@data[,i]) }
for(i in grep("krigeSt_", names(co_tif))){ co_tif@data[,i] = co_tif@data[,i]*100 }
#str(co_tif@data)
#p.ds = quantile(as.vector(co_tif@data[-grep("_pe", names(co_tif))]), probs = c(0.01,0.99), na.rm = TRUE)
p.ds = range(as.vector(co_tif@data[-grep("_pe", names(co_tif))]))
#p.de = quantile(as.vector(co_tif@data[grep("_pe", names(co_tif))]), probs = c(0.01,0.99), na.rm = TRUE)
p.de = range(as.vector(as.vector(co_tif@data[grep("_pe", names(co_tif))])))
#p.names = c(T.lst[1:4], paste0(T.lst[1:4], " (pe)"), T.lst[1:4], paste0(T.lst[1:4], " (pe)"))
p.names = rep(T.lst[1:4], 4)
colVec <- colorRampPalette(c("lightyellow","lightblue","blue","darkblue"))(100)

pdf(file = "RF_vs_kriging/results/st_prec/Fig_st_prec_predictions.pdf", width=10, height=7)
par(mfrow=c(2,4), oma=c(0,0,0,0), mar=c(0,0,4,0))
for(i in grep(glob2rx("co_PRCP_2016.02.??$"), names(co_tif))){
  plot(raster(co_tif[i]), col=colVec, zlim=p.ds, main=p.names[i], axes=FALSE, box=FALSE, legend=FALSE)
  points(co_locs.sp, pch="+")
}
for(i in grep(glob2rx("krigeSt_PRCP_2016.02.??$"), names(co_tif))){
  plot(raster(co_tif[i]), col=colVec, zlim=p.ds, main=p.names[i], axes=FALSE, box=FALSE, legend=FALSE)
  points(co_locs.sp, pch="+")
}
dev.off()

pdf(file = "RF_vs_kriging/results/st_prec/Fig_st_prec_predictions_pe.pdf", width=10, height=7)
par(mfrow=c(2,4), oma=c(0,0,0,0), mar=c(0,0,4,0))
for(i in grep(glob2rx("co_PRCP_2016.02.??_pe$"), names(co_tif))){
  plot(raster(co_tif[i]), col=rev(bpy.colors()), zlim=p.de, main=p.names[i], axes=FALSE, box=FALSE, legend=FALSE)
  points(co_locs.sp, pch="+")
}
for(i in grep(glob2rx("krigeSt_PRCP_2016.02.??_pe$"), names(co_tif))){
  plot(raster(co_tif[i]), col=rev(bpy.colors()), zlim=p.de, main=p.names[i], axes=FALSE, box=FALSE, legend=FALSE)
  points(co_locs.sp, pch="+")
}
#plot(raster(co_tif[i]), col=leg, zlim=p.de, legend.only=TRUE)
dev.off()
