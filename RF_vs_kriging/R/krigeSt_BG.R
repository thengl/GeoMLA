##
library(sp)
library(spacetime)
library(gstat)

co_prec =  readRDS("./RF_vs_kriging/data/st_prec/boulder_prcp.rds")

co_prec$cdate = floor(unclass(as.POSIXct(as.POSIXct(paste(co_prec$ DATE), format="%Y-%m-%d")))/86400)
co_prec$doy = as.integer(strftime(as.POSIXct(paste(co_prec$DATE), format="%Y-%m-%d"), format = "%j"))

co_locs.sp = co_prec[,c("STATION","LATITUDE","LONGITUDE")]
coordinates(co_locs.sp) = ~ LONGITUDE + LATITUDE
proj4string(co_locs.sp) = CRS("+proj=longlat +datum=WGS84")
co_grids = readRDS("./RF_vs_kriging/data/st_prec/boulder_grids.rds")
co_grids = as(co_grids, "SpatialPixelsDataFrame")
co_locs.sp = spTransform(co_locs.sp, co_grids@proj4string)

co_time <- as.POSIXct(co_prec$DATE)

stsdf <- STIDF(as(co_locs.sp, "SpatialPoints"), co_time, co_prec)
stsdf <- as(stsdf, "STSDF")

empStVgm <- gstat::variogramST(PRCP~1, stsdf, tlags = 0:3)
plot(empStVgm, wireframe=TRUE, scales=list(arrows=F))

## fit of st vgm

# rescale to ease optim
empStVgm$dist <- empStVgm$dist/1e3
empStVgm$avgDist <- empStVgm$avgDist/1e3

metFit <- fit.StVariogram(empStVgm, vgmST("metric", joint=vgm(0.035, "Sph", 30, 0.005), stAni=1))

smmFit <- fit.StVariogram(empStVgm, vgmST("sumMetric", 
                                          space=vgm(0.035, "Sph", 30, 0.05),
                                          time=vgm(0.035, "Sph", 60, 0.001),
                                          joint=vgm(0.035, "Sph", 30, 0.001),
                                          nugget=0.005, stAni=1))

plot(empStVgm, list(metFit, smmFit), all=T, wireframe=F)

plot(empStVgm, wireframe=T)
plot(empStVgm, smmFit, wireframe=T)

