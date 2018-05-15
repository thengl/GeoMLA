## Comparison RF vs kriging, see slides at: https://github.com/ISRICWorldSoil/GSIF_tutorials/blob/master/geul/5H_Hengl.pdf
## tom.hengl@gmail.com
## Cite as: Hengl et al., "Random Forest as Best Unbiased Generic Predictor (BUGP) for Spatial and Spatio-temporal Data", to be submitted to PeerJ Computer Science

list.of.packages <- c("plyr", "parallel", "randomForest", "quantregForest", "plotKML", "GSIF", "ranger", "RCurl", "raster", "rgdal", "geoR")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = TRUE)

setwd("~/git/GeoMLA/RF_vs_kriging")
load(".RData")
library(GSIF)
library(rgdal)
library(raster)
library(gstat)
library(randomForest)
library(quantregForest)
library(plotKML)
library(scales)
library(ranger)
library(RCurl)
library(parallel)
library(geoR)
#library(geostatsp)
leg = c("#0000ff", "#0028d7", "#0050af", "#007986", "#00a15e", "#00ca35", "#00f20d", "#1aff00", "#43ff00", "#6bff00", "#94ff00", "#bcff00", "#e5ff00", "#fff200", "#ffca00", "#ffa100", "#ff7900", "#ff5000", "#ff2800", "#ff0000")
source('BGUP_functions.R')

## Load the Meuse data set:
demo(meuse, echo=FALSE)

## Compare GLM vs RF ----
m <- glm(zinc~log1p(dist)+ffreq, meuse, family=gaussian(link=log))
plot(m$fitted.values~m$y, asp=1)
abline(0,1)
rf <- quantregForest(x=meuse@data[,c("dist","ffreq")], y=meuse$zinc)
plot(rf$predicted~rf$y, asp=1)
abline(0,1)
meuse.grid$glm.zinc <- predict(m, meuse.grid@data, type="response")
meuse.grid$rf.zinc <- predict(rf, meuse.grid@data)[,2]

## Plot predictions next to each other:
meuse.grid$glm.zinc = ifelse(meuse.grid$glm.zinc<expm1(4.8), expm1(4.8), meuse.grid$glm.zinc)
png(file = "Fig_comparison_GLM_RF_zinc_meuse.png", res = 150, width = 1750, height = 1200)
par(mfrow=c(1,2), oma=c(0,0,0,0))
plot(log1p(raster(meuse.grid["glm.zinc"])), col=leg, zlim=c(4.8,7.4), main="GLM")
points(meuse, pch="+")
plot(log1p(raster(meuse.grid["rf.zinc"])), col=leg, zlim=c(4.8,7.4), main="Random Forest")
points(meuse, pch="+")
dev.off()
## TH: Very similar

## Zinc predicted using ordinary kriging ----
zinc.geo <- as.geodata(meuse["zinc"])
#plot(variog4(zinc.geo, lambda=0, max.dist=1500, messages=FALSE), lwd=2)
zinc.vgm <- likfit(zinc.geo, lambda=0, messages=FALSE, ini=c(var(log1p(zinc.geo$data)),500), cov.model="exponential")
locs = meuse.grid@coords
zinc.ok <- krige.conv(zinc.geo, locations=locs, krige=krige.control(obj.model=zinc.vgm))
meuse.grid$zinc_ok = zinc.ok$predict
meuse.grid$zinc_ok_var = zinc.ok$krige.var

## Zinc predicted using RF and buffer distances only ----
grid.dist0 <- buffer.dist(meuse["zinc"], meuse.grid[1], as.factor(1:nrow(meuse)))
dn0 <- paste(names(grid.dist0), collapse="+")
fm0 <- as.formula(paste("zinc ~ ", dn0))
ov.zinc <- over(meuse["zinc"], grid.dist0)
m.zinc <- ranger(fm0, cbind(meuse@data["zinc"], ov.zinc), keep.inbag = TRUE)
zinc.rfd <- predict(m.zinc, grid.dist0@data, type = "se")
meuse.grid$zinc_rfd = zinc.rfd$predictions
meuse.grid$zinc_rfd_var = zinc.rfd$se

## Plot predictions next to each other:
#png(file = "Fig_comparison_OK_RF_zinc_meuse.png", res = 150, width = 1750, height = 1200)
var.max = max(c(meuse.grid$zinc_rfd_var, sqrt(meuse.grid$zinc_ok_var)))
axis.ls = list(at=c(4.8,5.7,6.5,7.4), labels=round(expm1(c(4.8,5.7,6.5,7.4))))
pdf(file = "Fig_comparison_OK_RF_zinc_meuse.pdf", width=9, height=9)
par(mfrow=c(2,2), oma=c(0,0,0,1), mar=c(0,0,4,3))
plot(log1p(raster(meuse.grid["zinc_ok"])), col=leg, zlim=c(4.8,7.4), main="Ordinary Kriging (OK)", axes=FALSE, box=FALSE, axis.args=axis.ls)
points(meuse, pch="+")
plot(log1p(raster(meuse.grid["zinc_rfd"])), col=leg, zlim=c(4.8,7.4), main="Random Forest (RF)", axes=FALSE, box=FALSE, axis.args=axis.ls)
points(meuse, pch="+")
plot(sqrt(raster(meuse.grid["zinc_ok_var"])), col=rev(bpy.colors()), zlim=c(0,var.max), main="OK prediction error", axes=FALSE, box=FALSE)
points(meuse, pch="+")
plot(raster(meuse.grid["zinc_rfd_var"]), col=rev(bpy.colors()), zlim=c(0,var.max), main="RF prediction error", axes=FALSE, box=FALSE)
points(meuse, pch="+")
dev.off()
## TH: RF smooths somewhat more than OK

## cross-validation (takes few seconds):
cv.RF = cv_numeric(varn="zinc", points=meuse, covs=meuse.grid, cpus=1, method="ranger", OK=TRUE, spcT=FALSE)
cv.OK = cv_numeric(varn="zinc", points=meuse, covs=meuse.grid, cpus=1, method="geoR", OK=TRUE, spcT=FALSE)
cv.RF$Summary$R.squared; cv.OK$Summary$R.squared
## plot results
library(lattice)
require(gridExtra)
lim.zinc = range(meuse$zinc, na.rm = TRUE)
pdf(file = "Fig_correlation_plots_OK_RF_zinc_meuse.pdf", width=9, height=5)
par(oma=c(0,0,0,1), mar=c(0,0,0,2))
plt.RF = xyplot(cv.RF[[1]]$Predicted~cv.RF[[1]]$Observed, asp=1, par.settings=list(plot.symbol = list(col=alpha("black", 0.6), fill=alpha("red", 0.6), pch=21, cex=0.9)), scales=list(x=list(log=TRUE, equispaced.log=FALSE), y=list(log=TRUE, equispaced.log=FALSE)), xlab="measured", ylab="predicted (machine learning)", panel = pfun.line, xlim=lim.zinc, ylim=lim.zinc)
plt.OK = xyplot(cv.OK[[1]]$Predicted~cv.OK[[1]]$Observed, asp=1, par.settings=list(plot.symbol = list(col=alpha("black", 0.6), fill=alpha("red", 0.6), pch=21, cex=0.9)), scales=list(x=list(log=TRUE, equispaced.log=FALSE), y=list(log=TRUE, equispaced.log=FALSE)), xlab="measured", ylab="predicted (geoR)", panel = pfun.line, xlim=lim.zinc, ylim=lim.zinc)
grid.arrange(plt.OK, plt.RF, ncol=2)
dev.off()

## RF with combined covariates ----
meuse.grid$SW_occurrence = readGDAL("Meuse_GlobalSurfaceWater_occurrence.tif")$band1[meuse.grid@grid.index]
meuse.grid$AHN = readGDAL("ahn.asc")$band1[meuse.grid@grid.index]
meuse.grid$LGN5 = as.factor(readGDAL("lgn5.asc")$band1[meuse.grid@grid.index])
grids.spc = spc(meuse.grid, as.formula("~ SW_occurrence + AHN + ffreq + dist"))
## fit hybrid RF model:
fm1 <- as.formula(paste("zinc ~ ", dn0, " + ", paste(names(grids.spc@predicted), collapse = "+")))
ov.zinc1 <- over(meuse["zinc"], grids.spc@predicted)
m1.zinc <- ranger(fm1, do.call(cbind, list(meuse@data["zinc"], ov.zinc, ov.zinc1)), keep.inbag = TRUE, importance = "impurity")
m1.zinc
zinc.rfd1 <- predict(m1.zinc, cbind(grid.dist0@data, grids.spc@predicted@data), type = "se")
meuse.grid$zinc_rfd1 = zinc.rfd1$predictions
meuse.grid$zinc_rfd1_var = zinc.rfd1$se
xl <- as.list(ranger::importance(m1.zinc))
print(t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:10]])))
m2.zinc <- ranger(paste("zinc ~ ", paste(names(grids.spc@predicted), collapse = "+")), do.call(cbind, list(meuse@data["zinc"], ov.zinc1)))
m2.zinc
meuse.grid$zinc_rfd2 = predict(m2.zinc, grids.spc@predicted@data)$predictions

pdf(file = "Fig_RF_covs_bufferdist_zinc_meuse.pdf", width=9, height=5)
par(mfrow=c(1,2), oma=c(0,0,0,1), mar=c(0,0,4,3))
plot(log1p(raster(meuse.grid["zinc_rfd2"])), col=leg, zlim=c(4.8,7.4), main="Random Forest (RF) covs only", axes=FALSE, box=FALSE, axis.args=axis.ls)
points(meuse, pch="+")
plot(log1p(raster(meuse.grid["zinc_rfd1"])), col=leg, zlim=c(4.8,7.4), main="Random Forest (RF) covs + buffer dist.", axes=FALSE, box=FALSE, axis.args=axis.ls)
points(meuse, pch="+")
dev.off()

## SIC 1997 data set ----
## measurements made in Switzerland on the 8th of May 1986
sic97.sp = readRDS("sic97.rds")
swiss1km = readRDS("swiss1km.rds")
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
plot(raster(swiss1km["rainfall_UK"]))
plot(sqrt(raster(swiss1km["rainfall_UK_var"])))

## SIC 1997 Random Forest example ----
swiss.dist0 <- buffer.dist(sic97.sp["rainfall"], swiss1km[1], as.factor(1:nrow(sic97.sp))) ## takes 2-3 mins
ov.swiss = over(sic97.sp["rainfall"], swiss.dist0)
sw.dn0 <- paste(names(swiss.dist0), collapse="+")
sw.fm1 <- as.formula(paste("rainfall ~ ", sw.dn0, " + CHELSA_rainfall + DEM"))
ov.rain <- over(sic97.sp["rainfall"], swiss1km[1:2])
sw.rm = do.call(cbind, list(sic97.sp@data["rainfall"], ov.rain, ov.swiss))
m1.rain <- ranger::ranger(sw.fm1, sw.rm[complete.cases(sw.rm),], keep.inbag = TRUE, importance = "impurity")
m1.rain
## Predicting errors is computationally very intensive in ranger
rain.rfd1 <- predict(m1.rain, cbind(swiss.dist0@data, swiss1km@data), type = "se")
swiss1km$rainfall_rfd1 = rain.rfd1$predictions
swiss1km$rainfall_rfd1_var = rain.rfd1$se
xl1 <- as.list(ranger::importance(m1.rain))
print(t(data.frame(xl1[order(unlist(xl1), decreasing=TRUE)[1:15]])))

rain.max = max(swiss1km$rainfall_rfd1, na.rm = TRUE)
swiss1km$rainfall_UK = ifelse(swiss1km$rainfall_UK<0, 0, ifelse(swiss1km$rainfall_UK>rain.max, rain.max, swiss1km$rainfall_UK))
## Plot predictions next to each other:
pdf(file = "Fig_Swiss_rainfall_UK_vs_RF.pdf", width=12, height=8)
par(mfrow=c(2,2), oma=c(0,0,0,0.5), mar=c(0,0,1.5,1))
plot(raster(swiss1km["rainfall_UK"]), col=leg, main="Universal kriging (UK)", axes=FALSE, box=FALSE, zlim=c(0, rain.max))
points(sic97.sp, pch="+")
plot(raster(swiss1km["rainfall_rfd1"]), col=leg, main="Random Forest (RF)", axes=FALSE, box=FALSE, zlim=c(0, rain.max))
points(sic97.sp, pch="+")
plot(sqrt(raster(swiss1km["rainfall_UK_var"])), col=rev(bpy.colors()), main="Universal kriging (UK) prediction error", axes=FALSE, box=FALSE, zlim=c(0,140))
points(sic97.sp, pch="+")
plot(raster(swiss1km["rainfall_rfd1_var"]), col=rev(bpy.colors()), main="Random Forest (RF) prediction error", axes=FALSE, box=FALSE, zlim=c(0,140))
points(sic97.sp, pch="+")
dev.off()

## cross-validation (takes few minutes!!):
cv.RF2 = cv_numeric(varn="rainfall", points=sic97.sp, covs=swiss1km[c("CHELSA_rainfall","DEM")], cpus=1, method="ranger", spcT=TRUE)
cv.UK = cv_numeric(varn="rainfall", points=sic97.sp, covs=swiss1km[c("CHELSA_rainfall","DEM")], cpus=1, method="geoR", spcT=FALSE)
cv.RF2$Summary$R.squared; cv.OK$Summary$R.squared
## plot results
library(lattice)
require(gridExtra)
lim.rain = c(20,max(sic97.sp$rainfall, na.rm = TRUE))
pdf(file = "Fig_correlation_plots_OK_RF_rain_SIC97.pdf", width=9, height=5)
par(oma=c(0,0,0,1), mar=c(0,0,0,2))
plt.RF2 = xyplot(cv.RF2[[1]]$Predicted~cv.RF2[[1]]$Observed, asp=1, par.settings=list(plot.symbol = list(col=alpha("black", 0.6), fill=alpha("red", 0.6), pch=21, cex=0.9)), scales=list(x=list(log=TRUE, equispaced.log=FALSE), y=list(log=TRUE, equispaced.log=FALSE)), xlab="measured", ylab="predicted (machine learning)", panel = pfun.line, xlim=lim.rain, ylim=lim.rain)
plt.UK = xyplot(cv.UK[[1]]$Predicted~cv.UK[[1]]$Observed, asp=1, par.settings=list(plot.symbol = list(col=alpha("black", 0.6), fill=alpha("red", 0.6), pch=21, cex=0.9)), scales=list(x=list(log=TRUE, equispaced.log=FALSE), y=list(log=TRUE, equispaced.log=FALSE)), xlab="measured", ylab="predicted (geoR)", panel = pfun.line, xlim=lim.rain, ylim=lim.rain)
grid.arrange(plt.UK, plt.RF2, ncol=2)
dev.off()


## Ebergotzen binomial variable ----
library(plotKML)
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
m1.Parabraunerde <- ranger::ranger(eb.fm1, rm.eberg2[complete.cases(rm.eberg2),], importance = "impurity", probability = TRUE)
m1.Parabraunerde
## Seems to be quite an accurate model
xl1.P <- as.list(ranger::importance(m1.Parabraunerde))
print(t(data.frame(xl1.P[order(unlist(xl1.P), decreasing=TRUE)[1:10]])))
## Predict:
pr.Parabraunerde = predict(m1.Parabraunerde, cbind(eberg.dist0@data, eberg_grid@data[paste0("PC", 1:10)]))
eberg_grid$Parabraunerde_TRUE = pr.Parabraunerde$predictions[,2]
eberg_grid$Parabraunerde_FALSE = pr.Parabraunerde$predictions[,1]
library(entropy)
eberg_grid$SSE = entropy_index(eberg_grid@data[,c("Parabraunerde_TRUE","Parabraunerde_FALSE")])

pdf(file = "Fig_Parabraunerde_RF.pdf", width=7, height=6.5)
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
pdf("Fig_ebergotzen_TAXGRSC.pdf", width=10, height=6.7)
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
kml(eberg, colour = TAXGRSC, file.name="eberg_TAXGRSC.kml", shape=shape, points_names=eberg$TAXGRSC, LabelScale = .9)
kml(raster(r.G), folder.name="Gley", raster_name="eberg_Gley.png", file.name="eberg_Gley.kml", colour_scale=SAGA_pal[["SG_COLORS_YELLOW_RED"]], zlim=c(0,40), png.type="cairo")
kml(raster(r.P), folder.name="Parabraunerde", raster_name="eberg_Parabraunerde.png", file.name="eberg_Parabraunerde.kml", colour_scale=SAGA_pal[["SG_COLORS_YELLOW_RED"]], zlim=c(0,40), png.type="cairo")
kml(eberg_grid["SSE_t"], folder.name="SSE", raster_name="eberg_SEE.png", file.name="eberg_SEE.kml", colour_scale=rev(bpy.colors()), zlim=c(0,100), png.type="cairo")
## Conclusion: looks like regression-kriging on class probs

## Multivariate case ----
## Geochemicals USA (https://mrdata.usgs.gov/geochem/)
geochem = readRDS("geochem.rds")
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
usa5km = readRDS("usa5km.rds")
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
pdf(file = "Fig_NGS_elements_RF.pdf", width=7.5, height=7)
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

nl.rd <- getURL("http://spatialreference.org/ref/sr-org/6781/proj4/")
## Geul data set ----
geul <- read.table("geul.dat", header = TRUE, as.is = TRUE)
geul$pb = as.numeric(geul$pb)
geul = geul[!is.na(geul$pb),]
coordinates(geul) <- ~x+y
proj4string(geul) <- CRS(nl.rd) 
grd25 <- readGDAL("dem25.txt")
grd25 <- as(grd25, "SpatialPixelsDataFrame")
proj4string(grd25) = proj4string(geul) 

## Pb predicted using OK
pb.geo <- as.geodata(geul["pb"])
pb.vgm <- likfit(pb.geo, lambda=0, messages=FALSE, ini=c(var(log1p(pb.geo$data)),500), cov.model="exponential")
locs2 = grd25@coords
pb.ok <- krige.conv(pb.geo, locations=locs2, krige=krige.control(obj.model=pb.vgm))
grd25$pb_ok = pb.ok$predict
## Pb predicted using RF only
ov.geul = over(geul["pb"], grd25)
summary(ov.geul$band1)
geul.s = geul[!is.na(ov.geul$band1),"pb"]
grid.dist1 <- buffer.dist(geul.s, grd25[1], as.factor(1:nrow(geul.s)))
dn1 <- paste(names(grid.dist1), collapse="+")
fm1 <- as.formula(paste("pb ~", dn1))
m1 <- fit.gstatModel(geul.s, fm1, grid.dist1, method="ranger", rvgm=NULL)
rk.m1 <- predict(m1, grid.dist1)

## Plot predictions next to each other:
grd25$pb_ok = ifelse(grd25$pb_ok<expm1(4.2), expm1(4.2), grd25$pb_ok)
png(file = "Fig_comparison_OK_RF_Pb_Geul.png", res = 150, width = 1750, height = 1200)
par(mfrow=c(1,2), oma=c(0,0,0,0))
plot(log1p(raster(grd25["pb_ok"])), col=leg, zlim=c(4.2,6.6), main="geoR (krige.conv)")
points(geul.s, pch="+")
plot(log1p(raster(rk.m1@predicted[2])), col=leg, zlim=c(4.2,6.6), main="Random Forest")
points(geul.s, pch="+")
dev.off()

## RF with both buffer dist and covariates ----
grd25$swi <- readGDAL("swi.sdat")$band1[grd25@grid.index]
grd25$dis <- readGDAL("riverdist.txt")$band1[grd25@grid.index]
plot(stack(grd25))
grd25T <- grd25[c("band1","swi","dis")]
grd25T@data <- cbind(grd25T@data, grid.dist1@data)
## Run principal component analysis:
grd25.spc <- spc(grd25T, as.formula(paste("~", paste(names(grd25T), collapse = "+"))))
plot(stack(grd25.spc@predicted[1:6]))
fm2 <- as.formula(paste("pb ~", paste(names(grd25.spc@predicted), collapse = "+")))
m2 <- fit.gstatModel(geul.s, fm2, grd25.spc@predicted, method="quantregForest", rvgm=NULL)
plot(m2)
dev.off()
rk.m2 <- predict(m2, grd25.spc@predicted)
#plot(rk.m2, col=leg)
varImpPlot(m1@regModel)

rk.m2@predicted$pb = ifelse(rk.m2@predicted$pb<expm1(4.2), expm1(4.2), rk.m2@predicted$pb)
png(file = "Fig_comparison_RF_covariates_Pb_Geul.png", res = 150, width = 1750, height = 1200)
par(mfrow=c(1,2), oma=c(0,0,0,0))
plot(log1p(raster(rk.m2@predicted[2])), col=leg, zlim=c(4.2,6.6), main="Random Forest + covs")
points(geul.s, pch="+")
plot(log1p(raster(rk.m1@predicted[2])), col=leg, zlim=c(4.2,6.6), main="Random Forest")
points(geul.s, pch="+")
dev.off()
