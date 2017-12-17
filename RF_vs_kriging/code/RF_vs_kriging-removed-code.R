
# Codes removed from Script RF_vs_kriging.R 
# by Madlene Nussbaum <madlene.nussbaum@bfh.ch>  
# Reason: these lines are not directly needed to produce the output for 
## Hengl et al., "Random Forest as a Generic Framework for Predictive Modeling of Spatial and Spatio-temporal Variables", to be submitted to PeerJ Computer Science


# Preamble 

#library(quantregForest)
#library(geostatsp)

load(".RData")



# Meuse ------------------


## Compare GLM vs RF --
m <- glm(zinc~log1p(dist)+ffreq, meuse, family=gaussian(link=log))
plot(m$fitted.values~m$y, asp=1)
abline(0,1)
set.seed(1)
rf <- quantregForest(x=meuse@data[,c("dist","ffreq")], y=meuse$zinc)
plot(rf$predicted~rf$y, asp=1)
abline(0,1)
meuse.grid$glm.zinc <- predict(m, meuse.grid@data, type="response")
meuse.grid$rf.zinc <- predict(rf, meuse.grid@data)[,2]
## Plot predictions next to each other:
meuse.grid$glm.zinc = ifelse(meuse.grid$glm.zinc<expm1(4.8), expm1(4.8), meuse.grid$glm.zinc)
png(file = "results/meuse/Fig_comparison_GLM_RF_zinc_meuse.png", res = 150, width = 1750, height = 1200)
par(mfrow=c(1,2), oma=c(0,0,0,0))
plot(log1p(raster(meuse.grid["glm.zinc"])), col=leg, zlim=c(4.8,7.4), main="GLM")
points(meuse, pch="+")
plot(log1p(raster(meuse.grid["rf.zinc"])), col=leg, zlim=c(4.8,7.4), main="Random Forest")
points(meuse, pch="+")
dev.off()
## TH: Very similar

#plot(variog4(zinc.geo, lambda=0, max.dist=1500, messages=FALSE), lwd=2)


# Prediction error

summary(meuse.grid$zinc_rfd_var)
## Median = 123.6
## Compare with the prediction error derived using the ranger package:
m.zinc.r <- ranger(fm0, rm.zinc, keep.inbag = TRUE)
zinc.rfdR <- predict(m.zinc.r, grid.dist0@data, type="se")
## Prediction error:
meuse.grid$zinc_rfdR_var =  sqrt(zinc.rfdR$se^2 + var(m.zinc.r$predictions-meuse$zinc))
summary(meuse.grid$zinc_rfdR_var)
## Much higher than what we get with the "quantregRanger"


# RF with coordinates only

## Plot
png(file = "results/meuse/Fig_RF_zinc_coordinates_only.png", res = 150, width = 1750, height = 1200)
par(mfrow=c(1,2), oma=c(0,0,0,0))
plot(log1p(raster(meuse.grid["zinc_rfc"])), col=leg, zlim=c(4.8,7.4), main="RF coordinates only", axes=FALSE, box=FALSE, axis.args=axis.ls)
points(meuse, pch="+")
plot(raster(meuse.grid["zinc_rfc_var"]), main="Prediction error", col=rev(bpy.colors()), axes=FALSE, box=FALSE)
points(meuse, pch="+")
dev.off()



## RF with noise data --

xN = cbind(meuse@data["zinc"], ov.zinc)
xN$zinc = runif(length(xN$zinc))
m.zincN <- quantregRanger(fm0, xN)
m.zincN
## Correclty R-square close to 0
zinc.rfdN <- predict(m.zincN, grid.dist0@data, quantiles)
meuse.grid$zinc_rfdN = zinc.rfdN[,2]
meuse.grid$zinc_rfd_varN = (zinc.rfdN[,3]-zinc.rfdN[,1])/2
summary(meuse.grid$zinc_rfd_varN)
## 2nd try:
xN$zinc = runif(length(xN$zinc))
m.zincN2 <- quantregRanger(fm0, xN)
zinc.rfdN2 <- predict(m.zincN2, grid.dist0@data, quantiles)
meuse.grid$zinc_rfdN2 = zinc.rfdN2[,2]
meuse.grid$zinc_rfd_varN2 = (zinc.rfdN2[,3]-zinc.rfdN2[,1])/2
## The expected mean and the se around the mean is:
mean(xN$zinc)
sqrt(var(xN$zinc - 0.5))
## 0.28
mean(meuse.grid$zinc_rfd_varN)
## 0.29
## Matches exactly s.d. from uniform distribution

## Plot 2 realizations:
png(file = "results/meuse/Fig_RF_zinc_pure_noise.png", res = 150, width = 1750, height = 1600)
par(mfrow=c(2,2), oma=c(0,0,0,1), mar=c(0,0,4,3))
plot(raster(meuse.grid["zinc_rfdN"]), col=leg, zlim=c(0,1), main="RF predictions (pure noise) 1", axes=FALSE, box=FALSE)
points(meuse, pch="+")
plot(raster(meuse.grid["zinc_rfd_varN"]), col=rev(bpy.colors()), main="RF prediction error 1", axes=FALSE, box=FALSE, zlim=c(0.1,.43))
points(meuse, pch="+")
plot(raster(meuse.grid["zinc_rfdN2"]), col=leg, zlim=c(0,1), main="RF predictions (pure noise) 2", axes=FALSE, box=FALSE)
points(meuse, pch="+")
plot(raster(meuse.grid["zinc_rfd_varN2"]), col=rev(bpy.colors()), main="RF prediction error 2", axes=FALSE, box=FALSE, zlim=c(0.1,.43))
points(meuse, pch="+")
dev.off()


