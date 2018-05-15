## Derivation of RF uncertainty (maps) for regression and classification type problems
## Based on: https://github.com/imbs-hl/ranger/issues/136 
## By: tom.hengl@gmail.com and Marvin Wright <marv@wrig.de>

setwd("~/git/GeoMLA/RF_vs_kriging")
#devtools::install_github("imbs-hl/ranger", ref="myquantreg")
library(ranger)
library(rgdal)
library(raster)
library(GSIF)
#library(tuneRF)
source('R/RFsp_functions.R')

## 1 s.d. quantiles
quantiles = c((1-.682)/2, 0.5, 1-(1-.682)/2)
leg = c("#0000ff", "#0028d7", "#0050af", "#007986", "#00a15e", "#00ca35", "#00f20d", "#1aff00", "#43ff00", "#6bff00", "#94ff00", "#bcff00", "#e5ff00", "#fff200", "#ffca00", "#ffa100", "#ff7900", "#ff5000", "#ff2800", "#ff0000")

## Meuse data set:
demo(meuse, echo=FALSE)
meuse.grid$SW_occurrence = readGDAL("data/meuse/Meuse_GlobalSurfaceWater_occurrence.tif")$band1[meuse.grid@grid.index] ## flooding occurrence
meuse.grid$AHN = readGDAL("data/meuse/ahn.asc")$band1[meuse.grid@grid.index] ## AHN.nl precise elevation
meuse.grid$LGN5 = as.factor(readGDAL("data/meuse/lgn5.asc")$band1[meuse.grid@grid.index]) ## land use class
## convert to indicators:
meuse.grid@data = cbind(meuse.grid@data, data.frame(model.matrix(~LGN5-1, meuse.grid@data)))

## Regression problems prediction error ----
grid.dist0 <- GSIF::buffer.dist(meuse["zinc"], meuse.grid[1], as.factor(1:nrow(meuse)))
fm1 <- as.formula(paste("zinc ~ ", paste(names(grid.dist0), collapse="+"), " + SW_occurrence + dist + ", paste(paste0("LGN5", levels(meuse.grid$LGN5)), collapse = "+")))
## regression matrix:
rm.zinc1 <- do.call(cbind, list(meuse@data["zinc"], over(meuse["zinc"], meuse.grid), over(meuse["zinc"], grid.dist0)))
## Prediction error - the quantreg approach:
m1.zinc <- ranger(fm1, rm.zinc1, mtry=22, num.trees=500, importance = "impurity", seed = 1, quantreg = TRUE)
m1.zinc
pred.zinc.rfq = predict(m1.zinc, cbind(meuse.grid@data, grid.dist0@data), type="quantiles", quantiles=quantiles)
meuse.grid$zinc_rfq_U = pred.zinc.rfq$predictions[,3]
meuse.grid$zinc_rfq_L = pred.zinc.rfq$predictions[,1]
## Assuming normal distribution of errors this should match 1 s.d. of the prediction error:
meuse.grid$zinc_rfq_r = (meuse.grid$zinc_rfq_U - meuse.grid$zinc_rfq_L)/2
hist(meuse.grid$zinc_rfq_r)
## compare OOB RMSE and mean s.d. of prediction error:
mean(meuse.grid$zinc_rfq_r, na.rm=TRUE); sqrt(m1.zinc$prediction.error)
## Conclusion: mean prediction error is smaller than the OOB RMSE but 'only' 25%

## Prediction error - the Jackknifing approach approach:
m2.zinc <- ranger(fm1, rm.zinc1, mtry=22, num.trees=500, seed = 1, keep.inbag = TRUE)
pred.zinc.rfj = predict(m2.zinc, cbind(meuse.grid@data, grid.dist0@data), type="se")
hist(pred.zinc.rfj$se)
## compare OOB RMSE and mean s.d. of prediction error:
mean(pred.zinc.rfj$se, na.rm=TRUE); sqrt(m2.zinc$prediction.error)
## Conclusion: mean prediction error is much smaller than OOB RMSE
## TH: Solution is to scale it to RMSE?
meuse.grid$zinc_rfj_r = pred.zinc.rfj$se * sqrt(m2.zinc$prediction.error)/mean(pred.zinc.rfj$se, na.rm=TRUE)

## Plot two maps of errors next to each other:
r.max = quantile(c(meuse.grid$zinc_rfq_r, meuse.grid$zinc_rfj_r), probs=c(0.025, 0.975), na.rm=TRUE)
meuse.grid$zinc_rfq_r = ifelse(meuse.grid$zinc_rfq_r<r.max[1], r.max[1], ifelse(meuse.grid$zinc_rfq_r>r.max[2], r.max[2], meuse.grid$zinc_rfq_r))
meuse.grid$zinc_rfj_r = ifelse(meuse.grid$zinc_rfj_r<r.max[1], r.max[1], ifelse(meuse.grid$zinc_rfj_r>r.max[2], r.max[2], meuse.grid$zinc_rfj_r))

pdf(file = "results/meuse/Fig_comparison_uncertainty_quantreg_vs_jackknife_meuse.pdf", width=9.5, height=5.5)
par(mfrow=c(1,2), oma=c(0,0,0,0.5), mar=c(0,0,1.5,1))
par(oma=c(0,0,0,0.5), mar=c(0,0,3.5,1))
plot(raster(meuse.grid["zinc_rfq_r"]), col=rev(bpy.colors())[1:80], main="Prediction error RF quantreg", axes=FALSE, box=FALSE, zlim=r.max)
points(meuse, pch="+")
plot(raster(meuse.grid["zinc_rfj_r"]), col=rev(bpy.colors())[1:80], main="Prediction error RF Jackknife", axes=FALSE, box=FALSE, zlim=r.max)
points(meuse, pch="+")
dev.off()

## Binomial variable Prediction error ----
meuse@data = cbind(meuse@data, data.frame(model.matrix(~soil-1, meuse@data)))
summary(as.factor(meuse$soil1))
fm.s1 = as.formula(paste("soil1 ~ ", paste(names(grid.dist0), collapse="+"), " + SW_occurrence + dist"))
rm.s1 <- do.call(cbind, list(meuse@data["soil1"], over(meuse["soil1"], meuse.grid), over(meuse["soil1"], grid.dist0)))
## Option 1: treat binomial variable as numeric variable
m1.s1 <- ranger(fm.s1, rm.s1, mtry=22, num.trees=500, seed = 1, quantreg=TRUE)
fm.s1c <- as.formula(paste("soil1c ~ ", paste(names(grid.dist0), collapse="+"), " + SW_occurrence + dist"))
rm.s1$soil1c = as.factor(rm.s1$soil1)
## Option 2: treat binomial variable as factor variable
m2.s1 <- ranger(fm.s1c, rm.s1, mtry=22, num.trees=500, seed = 1, probability=TRUE, keep.inbag = TRUE)
## Conclusion: OOB MSE is about the same as the OOB prediction error
## Derive prediction errors:
pred.soil1_rfb = predict(m1.s1, cbind(meuse.grid@data, grid.dist0@data), type="quantiles", quantiles=quantiles)
meuse.grid$soil1_rfq_U = pred.soil1_rfb$predictions[,3]
meuse.grid$soil1_rfq = pred.soil1_rfb$predictions[,2]
meuse.grid$soil1_rfq_L = pred.soil1_rfb$predictions[,1]
## Assuming normal distribution of errors this should match 1 s.d. of the prediction error:
meuse.grid$soil1_rfq_r = (meuse.grid$soil1_rfq_U - meuse.grid$soil1_rfq_L)/2
mean(meuse.grid$soil1_rfq_r, na.rm=TRUE); sqrt(m1.s1$prediction.error)
## Again, quantreg error smaller than the OOB RMSE
pred.soil1_rfc = predict(m2.s1, cbind(meuse.grid@data, grid.dist0@data), type="se")
meuse.grid$soil1_rfc = pred.soil1_rfc$predictions[,2]
mean(pred.soil1_rfc$se[,2], na.rm=TRUE); sqrt(m2.s1$prediction.error)
## Again much smaller and needs to be scaled:
meuse.grid$soil1_rfc_r = pred.soil1_rfc$se[,2] * sqrt(m2.s1$prediction.error)/mean(pred.soil1_rfc$se[,2], na.rm=TRUE)

r.soil1 = quantile(c(meuse.grid$soil1_rfq_r, meuse.grid$soil1_rfc_r), probs=c(0.025, 0.975), na.rm=TRUE)
meuse.grid$soil1_rfq_r = ifelse(meuse.grid$soil1_rfq_r<r.soil1[1], r.soil1[1], ifelse(meuse.grid$soil1_rfq_r>r.soil1[2], r.soil1[2], meuse.grid$soil1_rfq_r))
meuse.grid$soil1_rfc_r = ifelse(meuse.grid$soil1_rfc_r<r.soil1[1], r.soil1[1], ifelse(meuse.grid$soil1_rfc_r>r.soil1[2], r.soil1[2], meuse.grid$soil1_rfc_r))

# Regression prediction
pred.regr <- predict(m1.s1, cbind(meuse.grid@data, grid.dist0@data), type="response")$predictions
meuse.grid$soil1_rfr <- pred.regr

pdf(file = "results/meuse/Fig_comparison_uncertainty_Binomial_variables_meuse.pdf", width=7, height=8)
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
dev.off()

## Factor variable Prediction error ----
fm.s = as.formula(paste("soil ~ ", paste(names(grid.dist0), collapse="+"), " + SW_occurrence + dist"))
rm.s <- do.call(cbind, list(meuse@data["soil"], over(meuse["soil"], meuse.grid), over(meuse["soil"], grid.dist0)))
m.s <- ranger(fm.s, rm.s, mtry=22, num.trees=500, seed = 1, probability=TRUE, keep.inbag = TRUE)
pred.soil_rfc = predict(m.s, cbind(meuse.grid@data, grid.dist0@data), type="se")
pred.grids = meuse.grid["soil"]
pred.grids@data = do.call(cbind, list(pred.grids@data, data.frame(pred.soil_rfc$predictions), data.frame(pred.soil_rfc$se)))
names(pred.grids) = c("soil", paste0("pred_soil", 1:3), paste0("se_soil", 1:3))
str(pred.grids@data)
## scale errors?

pdf(file = "results/meuse/Fig_comparison_uncertainty_Factor_variables_meuse.pdf", width=9.5, height=7.5)
par(mfrow=c(2,3), oma=c(0,0,0,0.5), mar=c(0,0,1.5,1))
par(oma=c(0,0,0,0.5), mar=c(0,0,3.5,1))
plot(raster(pred.grids["pred_soil1"]), col=leg, main="soil type '1' RF probs", axes=FALSE, box=FALSE, zlim=c(0,1))
points(meuse["soil"], pch="+")
plot(raster(pred.grids["pred_soil2"]), col=leg, main="soil type '2' RF probs", axes=FALSE, box=FALSE, zlim=c(0,1))
points(meuse["soil"], pch="+")
plot(raster(pred.grids["pred_soil3"]), col=leg, main="soil type '3' RF probs", axes=FALSE, box=FALSE, zlim=c(0,1))
points(meuse["soil"], pch="+")
plot(raster(pred.grids["se_soil1"]), col=rev(bpy.colors())[1:80], main="prediction error soil type '1' RF", axes=FALSE, box=FALSE, zlim=c(0,.35))
points(meuse["soil"], pch="+")
plot(raster(pred.grids["se_soil2"]), col=rev(bpy.colors())[1:80], main="prediction error soil type '2' RF", axes=FALSE, box=FALSE, zlim=c(0,.35))
points(meuse["soil"], pch="+")
plot(raster(pred.grids["se_soil3"]), col=rev(bpy.colors())[1:80], main="prediction error soil type '3' RF", axes=FALSE, box=FALSE, zlim=c(0,.35))
points(meuse["soil"], pch="+")
dev.off()
