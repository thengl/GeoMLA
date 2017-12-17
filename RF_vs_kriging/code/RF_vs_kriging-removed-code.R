
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


## ** SIC 1997 data set ---------------------------------------------------

#plot(stack(swiss1km[1:2]))

#plot(raster(swiss1km["rainfall_UK"]))
#plot(sqrt(raster(swiss1km["rainfall_UK_var"])))


xl1 <- as.list(m1.rain$variable.importance)
# print(t(data.frame(xl1[order(unlist(xl1), decreasing=TRUE)[1:15]])))

# importance is not in the paper, saw it too late, MN
pdf(file = "results/rainfall/Fig_Swiss_rainfall_RF_covs_covar-importance.pdf", width=6, height=4.5)
par(mfrow=c(1,1),oma=c(0.7,2,0,1), mar=c(4,3.5,1,0))
plot(vv <- t(data.frame(xl1[order(unlist(xl1), decreasing=TRUE)[15:1]])), 1:15, type = "n", ylab = "", yaxt = "n", xlab = "Variable Importance (Node Impurity)")
abline(h = 1:15, lty = "dotted", col = "grey60")
points(vv, 1:15)
axis(2, 1:15, labels = dimnames(vv)[[1]], las = 2)
dev.off()


## ** Ebergotzen binomial variable --------------------------------------



#par(mfrow=c(1,2), oma=c(0,0,0,0.5), mar=c(0,0,1.5,1))
#plot(raster(eberg_grid["SSE"]), col=rev(bpy.colors()), main="Shannon Scaled Entropy Index", axes=FALSE, box=FALSE)
#points(eberg[eberg$Parabraunerde=="TRUE"&!is.na(eberg$Parabraunerde)&sel.eberg,], pch=19, cex=.8)
#points(eberg[eberg$Parabraunerde=="FALSE"&!is.na(eberg$Parabraunerde)&sel.eberg,], pch="+", cex=.8)


# This is not in the paper:
  

## ** Intamap example ---------------------------------------------------
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


## ** Geochemical USA - Multivariate case -----------------------------------

bubble(geochem[!is.na(geochem$PB_ICP40),"PB_ICP40"])
bubble(geochem[!is.na(geochem$CU_ICP40),"CU_ICP40"])
bubble(geochem[!is.na(geochem$K_ICP40),"K_ICP40"])
bubble(geochem[!is.na(geochem$MG_ICP40),"MG_ICP40"])

#writeOGR(geochem, "geochem.shp", "geochem", "ESRI Shapefile")

#str(rm.geochem)
hist(log1p(rm.geochem$Y))
summary(as.factor(rm.geochem$TYPE))


## ** Spatiotemporal prediction using Random Forest ------------------------

hist(co_prec$cdate)
hist(co_prec$doy)

plot(raster(co_grids[1]))
points(co_locs.sp, pch="+")


#bubble(co_locs.sp[!is.na(co_locs.sp$x),"x"])

#plot(stack(grid.distP[1:3]))
