## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- eval=FALSE, echo=TRUE----------------------------------------------
## devtools::install_github("imbs-hl/ranger")

## ---- echo=TRUE----------------------------------------------------------
library(GSIF)
library(rgdal)
library(raster)
library(ranger)
quantiles = c((1-.682)/2, 0.5, 1-(1-.682)/2)
## color legend:
leg = c("#0000ff", "#0028d7", "#0050af", "#007986", "#00a15e", "#00ca35", 
        "#00f20d", "#1aff00", "#43ff00", "#6bff00", "#94ff00", "#bcff00", 
        "#e5ff00", "#fff200", "#ffca00", "#ffa100", "#ff7900", "#ff5000", "#ff2800", "#ff0000")

## ----meuse---------------------------------------------------------------
demo(meuse, echo=FALSE)

## ---- results=FALSE, tidy=TRUE-------------------------------------------
dir.meuse = "../RF_vs_kriging/data/meuse/"
meuse.grid$SW_occurrence = readGDAL(paste0(dir.meuse, "Meuse_GlobalSurfaceWater_occurrence.tif"))$band1[meuse.grid@grid.index] 
## flooding occurrence
meuse.grid$AHN = readGDAL(paste0(dir.meuse, "ahn.asc"))$band1[meuse.grid@grid.index] 
## AHN.nl precise elevation
meuse.grid$LGN5 = as.factor(readGDAL(paste0(dir.meuse, "lgn5.asc"))$band1[meuse.grid@grid.index]) 
## land use class
## convert to indicators:
meuse.grid@data = cbind(meuse.grid@data, data.frame(model.matrix(~LGN5-1, meuse.grid@data)))

## ----bufferdist----------------------------------------------------------
grid.dist0 <- GSIF::buffer.dist(meuse["zinc"], meuse.grid[1], as.factor(1:nrow(meuse)))

## ------------------------------------------------------------------------
fm1 <- as.formula(paste("zinc ~ ", paste(names(grid.dist0), collapse="+"), 
      " + SW_occurrence + dist + ", paste(paste0("LGN5", levels(meuse.grid$LGN5)), collapse = "+")))
fm1

## ------------------------------------------------------------------------
rm.zinc1 <- do.call(cbind, list(meuse@data["zinc"], 
                             over(meuse["zinc"], meuse.grid), 
                             over(meuse["zinc"], grid.dist0)))

## ------------------------------------------------------------------------
m1.zinc <- ranger(fm1, rm.zinc1, mtry=22, num.trees=500, 
                   importance="impurity", seed=1, quantreg= TRUE)
m1.zinc

## ------------------------------------------------------------------------
xl <- as.list(ranger::importance(m1.zinc))
print(t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:10]])))

## ------------------------------------------------------------------------
pred.zinc.rfq = predict(m1.zinc, 
                        cbind(meuse.grid@data, grid.dist0@data), 
                        type="quantiles", quantiles=quantiles)
str(pred.zinc.rfq)

## ------------------------------------------------------------------------
pred.zinc.rfq$predictions[1,]

## ------------------------------------------------------------------------
meuse.grid$zinc_rfq_U = pred.zinc.rfq$predictions[,3]
meuse.grid$zinc_rfq_L = pred.zinc.rfq$predictions[,1]

## ----rfq-histogram, echo=FALSE, fig.width=6, fig.cap="Histogram of s.d. of the prediction error estimated using QRF."----
meuse.grid$zinc_rfq_r = (meuse.grid$zinc_rfq_U - meuse.grid$zinc_rfq_L)/2
hist(meuse.grid$zinc_rfq_r, main="QRF s.d. of prediction errors", col="grey")

## ------------------------------------------------------------------------
mean(meuse.grid$zinc_rfq_r, na.rm=TRUE); sqrt(m1.zinc$prediction.error)

## ------------------------------------------------------------------------
m2.zinc <- ranger(fm1, rm.zinc1, mtry=22, num.trees=500, seed=1, keep.inbag=TRUE)

## ------------------------------------------------------------------------
pred.zinc.rfj = predict(m2.zinc, cbind(meuse.grid@data, grid.dist0@data), type="se")
str(pred.zinc.rfj)

## ------------------------------------------------------------------------
mean(pred.zinc.rfj$se, na.rm=TRUE); sqrt(m2.zinc$prediction.error)

## ------------------------------------------------------------------------
meuse.grid$zinc_rfj_r = pred.zinc.rfj$se * 
    sqrt(m2.zinc$prediction.error)/mean(pred.zinc.rfj$se, na.rm=TRUE)

## ----jacknife-meuse-maps, echo=FALSE, fig.width=9, fig.cap="Comparison of uncertainty maps based on the QRF vs Jackknife approaches for the Meuse data set."----
r.max = quantile(c(meuse.grid$zinc_rfq_r, meuse.grid$zinc_rfj_r), probs=c(0.025, 0.975), na.rm=TRUE)
meuse.grid$zinc_rfq_r = ifelse(meuse.grid$zinc_rfq_r<r.max[1], r.max[1], ifelse(meuse.grid$zinc_rfq_r>r.max[2], r.max[2], meuse.grid$zinc_rfq_r))
meuse.grid$zinc_rfj_r = ifelse(meuse.grid$zinc_rfj_r<r.max[1], r.max[1], ifelse(meuse.grid$zinc_rfj_r>r.max[2], r.max[2], meuse.grid$zinc_rfj_r))
par(mfrow=c(1,2), oma=c(0,0,0,0.5), mar=c(0,0,1.5,1))
par(oma=c(0,0,0,0.5), mar=c(0,0,3.5,1))
plot(raster(meuse.grid["zinc_rfq_r"]), col=rev(bpy.colors())[1:80], main="Prediction error RF quantreg", axes=FALSE, box=FALSE, zlim=r.max)
points(meuse, pch="+")
plot(raster(meuse.grid["zinc_rfj_r"]), col=rev(bpy.colors())[1:80], main="Prediction error RF Jackknife", axes=FALSE, box=FALSE, zlim=r.max)
points(meuse, pch="+")

## ------------------------------------------------------------------------
sic97.sp = readRDS("../RF_vs_kriging/data/rainfall/sic97.rds")
swiss1km = readRDS("../RF_vs_kriging/data/rainfall/swiss1km.rds")
ov2 = over(y=swiss1km, x=sic97.sp)

## ------------------------------------------------------------------------
swiss.dist0 <- GSIF::buffer.dist(sic97.sp["rainfall"], 
                                 swiss1km[1], as.factor(1:nrow(sic97.sp))) 
## takes 1+ mins!
ov.swiss = over(sic97.sp["rainfall"], swiss.dist0)
sw.dn0 <- paste(names(swiss.dist0), collapse="+")
sw.fm1 <- as.formula(paste("rainfall ~ ", sw.dn0, " + CHELSA_rainfall + DEM"))
sw.fm1
ov.rain <- over(sic97.sp["rainfall"], swiss1km[1:2])
sw.rm = do.call(cbind, list(sic97.sp@data["rainfall"], ov.rain, ov.swiss))

## ------------------------------------------------------------------------
m1.rain <- ranger(sw.fm1, sw.rm[complete.cases(sw.rm),], mtry=27, 
                  min.node.size=2, sample.fraction=0.9930754, 
                  num.trees=150, importance = "impurity", seed=1, quantreg=TRUE)
m1.rain

## ------------------------------------------------------------------------
rain.rfd1 <- predict(m1.rain, cbind(swiss.dist0@data, swiss1km@data), 
                     type="quantiles", quantiles=quantiles)$predictions
## now more computational...
swiss1km$rainfall_rfd1 = rain.rfd1[,2]
## s.d. of the prediction error:
swiss1km$rainfall_rfd1_var = (rain.rfd1[,3]-rain.rfd1[,1])/2
str(swiss1km@data)

## ---- eval=FALSE, include=FALSE------------------------------------------
## library(rgdal)
## writeOGR(sic97.sp["rainfall"], dsn = "../RF_vs_kriging/results/rainfall/rainfall_sic97.shp", "rainfall_sic97", "ESRI Shapefile")
## writeGDAL(swiss1km["rainfall_rfd1"], fname = "../RF_vs_kriging/results/rainfall/pred_rainfall_RFsp.tif", options = c("COMPRESS=DEFLATE"))
## writeGDAL(swiss1km["rainfall_rfd1_var"], fname = "../RF_vs_kriging/results/rainfall/pred.var_rainfall_RFsp.tif", options = c("COMPRESS=DEFLATE"))

## ----qrf-sic97-maps, echo=FALSE, fig.width=9, fig.cap="Predictions and prediction errors for the SIC97 data set."----
par(mfrow=c(1,2), oma=c(0,0,0,0.5), mar=c(0,0,1.5,1))
plot(raster(swiss1km["rainfall_rfd1"]), col=leg, 
     main="Random Forest (RF)", axes=FALSE, box=FALSE)
points(sic97.sp, pch="+")
plot(raster(swiss1km["rainfall_rfd1_var"]), col=rev(bpy.colors()), 
     main="Random Forest (RF) prediction error", axes=FALSE, box=FALSE)
points(sic97.sp, pch="+")

## ------------------------------------------------------------------------
summary(meuse$soil)

## ------------------------------------------------------------------------
meuse@data = cbind(meuse@data, data.frame(model.matrix(~soil-1, meuse@data)))
summary(as.factor(meuse$soil1))

## ------------------------------------------------------------------------
fm.s1 = as.formula(paste("soil1 ~ ", 
                  paste(names(grid.dist0), collapse="+"), 
                  " + SW_occurrence + dist"))
fm.s1
rm.s1 <- do.call(cbind, list(meuse@data["soil1"], 
                             over(meuse["soil1"], meuse.grid), 
                             over(meuse["soil1"], grid.dist0)))

## ------------------------------------------------------------------------
m1.s1 <- ranger(fm.s1, rm.s1, mtry=22, num.trees=500, seed = 1, quantreg=TRUE)
m1.s1

## ------------------------------------------------------------------------
rm.s1$soil1c = as.factor(rm.s1$soil1)
summary(rm.s1$soil1c)

## ------------------------------------------------------------------------
fm.s1c <- as.formula(paste("soil1c ~ ", paste(names(grid.dist0), collapse="+"), 
                           " + SW_occurrence + dist"))
m2.s1 <- ranger(fm.s1c, rm.s1, mtry=22, num.trees=500, 
                seed=1, probability=TRUE, keep.inbag=TRUE)
m2.s1

## ------------------------------------------------------------------------
pred.soil1_rfb = predict(m1.s1, 
                         cbind(meuse.grid@data, grid.dist0@data), 
                         type="quantiles", quantiles=quantiles)
str(pred.soil1_rfb)

## ------------------------------------------------------------------------
meuse.grid$soil1_rfq_U = pred.soil1_rfb$predictions[,3]
meuse.grid$soil1_rfq = pred.soil1_rfb$predictions[,2]
meuse.grid$soil1_rfq_L = pred.soil1_rfb$predictions[,1]
meuse.grid$soil1_rfq_r = (meuse.grid$soil1_rfq_U - meuse.grid$soil1_rfq_L)/2
mean(meuse.grid$soil1_rfq_r, na.rm=TRUE); sqrt(m1.s1$prediction.error)

## ------------------------------------------------------------------------
pred.soil1_rfc = predict(m2.s1, cbind(meuse.grid@data, grid.dist0@data), type="se")
meuse.grid$soil1_rfc = pred.soil1_rfc$predictions[,2]

## ------------------------------------------------------------------------
meuse.grid$soil1_rfc_r = pred.soil1_rfc$se[,2] *
  sqrt(m2.s1$prediction.error)/mean(pred.soil1_rfc$se[,2], na.rm=TRUE)

## ------------------------------------------------------------------------
pred.regr <- predict(m1.s1, cbind(meuse.grid@data, grid.dist0@data), type="response")$predictions
meuse.grid$soil1_rfr <- pred.regr

## ----binomial-meuse-maps, echo=FALSE, fig.width=9, fig.cap="Comparison of uncertainty maps for a binomial variable based on the QRF vs Jackknife approaches for the Meuse data set."----
r.soil1 = quantile(c(meuse.grid$soil1_rfq_r, meuse.grid$soil1_rfc_r), probs=c(0.025, 0.975), na.rm=TRUE)
meuse.grid$soil1_rfq_r = ifelse(meuse.grid$soil1_rfq_r<r.soil1[1], 
                                r.soil1[1], ifelse(meuse.grid$soil1_rfq_r>r.soil1[2],
                                                   r.soil1[2], meuse.grid$soil1_rfq_r))
meuse.grid$soil1_rfc_r = ifelse(meuse.grid$soil1_rfc_r<r.soil1[1], 
                                r.soil1[1], ifelse(meuse.grid$soil1_rfc_r>r.soil1[2],
                                                   r.soil1[2], meuse.grid$soil1_rfc_r))
par(mfrow=c(2,3), oma=c(0,0,0,0.5), mar=c(0,0,1.5,1))
par(oma=c(0,0,0,0.5), mar=c(0,0,3.5,1))
plot(raster(meuse.grid["soil1_rfq"]), col=leg, main="Predictions soil '1' RF quantreg", axes=FALSE, box=FALSE, zlim=c(0,1))
points(meuse, pch="+")
plot(raster(meuse.grid["soil1_rfc"]), col=leg, main="Predictions soil '1' RF probs", axes=FALSE, box=FALSE, zlim=c(0,1))
points(meuse, pch="+")
plot(raster(meuse.grid["soil1_rfr"]), col=leg, main="Predictions soil '1' RF regr", axes=FALSE, box=FALSE, zlim=c(0,1))
points(meuse, pch="+")
plot(raster(meuse.grid["soil1_rfq_r"]), col=rev(bpy.colors())[1:80], main="Prediction error RF quantreg", axes=FALSE, box=FALSE, zlim=r.soil1)
points(meuse, pch="+")
plot(raster(meuse.grid["soil1_rfc_r"]), col=rev(bpy.colors())[1:80], main="Prediction error RF Jackknife", axes=FALSE, box=FALSE, zlim=r.soil1)
points(meuse, pch="+")

## ------------------------------------------------------------------------
summary(meuse$soil)

## ------------------------------------------------------------------------
fm.s = as.formula(paste("soil ~ ", paste(names(grid.dist0), collapse="+"), 
                        " + SW_occurrence + dist"))
rm.s <- do.call(cbind, list(meuse@data["soil"], 
                            over(meuse["soil"], meuse.grid), 
                            over(meuse["soil"], grid.dist0)))
m.s <- ranger(fm.s, rm.s, mtry=22, num.trees=500, seed=1, probability=TRUE, keep.inbag=TRUE)
m.s

## ------------------------------------------------------------------------
m.s0 <- ranger(fm.s, rm.s, mtry=22, num.trees=150, seed=1)
m.s0

## ------------------------------------------------------------------------
pred.soil_rfc = predict(m.s, cbind(meuse.grid@data, grid.dist0@data), type="se")

## ------------------------------------------------------------------------
pred.grids = meuse.grid["soil"]
pred.grids@data = do.call(cbind, list(pred.grids@data, 
                                      data.frame(pred.soil_rfc$predictions),
                                      data.frame(pred.soil_rfc$se)))
names(pred.grids) = c("soil", paste0("pred_soil", 1:3), paste0("se_soil", 1:3))
str(pred.grids@data)

## ----factor-meuse-maps, echo=FALSE, fig.width=9, fig.cap="Predictions of soil types for the meuse data set based on the RFsp: (above) probability for three soil classes, and (below) derived standard errors per class."----
par(mfrow=c(2,3), oma=c(0,0,0,0.5), mar=c(0,0,1.5,1))
par(oma=c(0,0,0,0.5), mar=c(0,0,3.5,1))
plot(raster(pred.grids["pred_soil1"]), col=leg, main="soil type '1' RF probs", 
     axes=FALSE, box=FALSE, zlim=c(0,1))
points(meuse["soil"], pch="+")
plot(raster(pred.grids["pred_soil2"]), col=leg, main="soil type '2' RF probs",
     axes=FALSE, box=FALSE, zlim=c(0,1))
points(meuse["soil"], pch="+")
plot(raster(pred.grids["pred_soil3"]), col=leg, main="soil type '3' RF probs",
     axes=FALSE, box=FALSE, zlim=c(0,1))
points(meuse["soil"], pch="+")
plot(raster(pred.grids["se_soil1"]), col=rev(bpy.colors())[1:80], 
     main="prediction error soil type '1' RF", axes=FALSE, box=FALSE, zlim=c(0,.35))
points(meuse["soil"], pch="+")
plot(raster(pred.grids["se_soil2"]), col=rev(bpy.colors())[1:80], 
     main="prediction error soil type '2' RF", axes=FALSE, box=FALSE, zlim=c(0,.35))
points(meuse["soil"], pch="+")
plot(raster(pred.grids["se_soil3"]), col=rev(bpy.colors())[1:80], 
     main="prediction error soil type '3' RF", axes=FALSE, box=FALSE, zlim=c(0,.35))
points(meuse["soil"], pch="+")

## ------------------------------------------------------------------------
library(mlr)
spatial.taskmeuse = makeRegrTask(data = rm.zinc1[,c("zinc","SW_occurrence","dist","AHN")], target = "zinc", coordinates = data.frame(meuse@coords))
spatial.taskmeuse
learner.rf = makeLearner("regr.ranger")
library("parallelMap")
parallelStartSocket(parallel::detectCores())
resampling = makeResampleDesc("SpRepCV", fold = 5, reps = 5)
cv.meuse = mlr::resample(learner = learner.rf, task = spatial.taskmeuse, resampling = resampling)
## compare with non-spatial CV:
nonspatial.taskmeuse = makeRegrTask(data = rm.zinc1[,c("zinc","SW_occurrence","dist","AHN")], target = "zinc")
resampling0 = makeResampleDesc("RepCV", fold = 5, reps = 5)
cv.meuse0 = mlr::resample(learner = learner.rf, task = nonspatial.taskmeuse, resampling = resampling0)
parallelStop()

