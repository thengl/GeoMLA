## Extra examples

library(RCurl)
library(rgdal)
nl.rd <- getURL("http://spatialreference.org/ref/sr-org/6781/proj4/")
## Geul data set ----
geul <- read.table("data/geul/geul.dat", header = TRUE, as.is = TRUE)
geul$pb = as.numeric(geul$pb)
geul = geul[!is.na(geul$pb),]
coordinates(geul) <- ~x+y
proj4string(geul) <- CRS(nl.rd) 
grd25 <- readGDAL("data/geul/dem25.txt")
grd25 <- as(grd25, "SpatialPixelsDataFrame")
proj4string(grd25) = proj4string(geul) 

## Pb predicted using OK
library(geoR)
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
png(file = "results/geul/Fig_comparison_OK_RF_Pb_Geul.png", res = 150, width = 1750, height = 1200)
par(mfrow=c(1,2), oma=c(0,0,0,0))
plot(log1p(raster(grd25["pb_ok"])), col=leg, zlim=c(4.2,6.6), main="geoR (krige.conv)")
points(geul.s, pch="+")
plot(log1p(raster(rk.m1@predicted[2])), col=leg, zlim=c(4.2,6.6), main="Random Forest")
points(geul.s, pch="+")
dev.off()

## RF with both buffer dist and covariates ----
grd25$swi <- readGDAL("data/geul/swi.sdat")$band1[grd25@grid.index]
grd25$dis <- readGDAL("data/geul/riverdist.txt")$band1[grd25@grid.index]
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
# varImpPlot(m1@regModel)

rk.m2@predicted$pb = ifelse(rk.m2@predicted$pb<expm1(4.2), expm1(4.2), rk.m2@predicted$pb)
png(file = "results/geul/Fig_comparison_RF_covariates_Pb_Geul.png", res = 150, width = 1750, height = 1200)
par(mfrow=c(1,2), oma=c(0,0,0,0))
plot(log1p(raster(rk.m2@predicted[2])), col=leg, zlim=c(4.2,6.6), main="Random Forest + covs")
points(geul.s, pch="+")
plot(log1p(raster(rk.m1@predicted[2])), col=leg, zlim=c(4.2,6.6), main="Random Forest")
points(geul.s, pch="+")
dev.off()
