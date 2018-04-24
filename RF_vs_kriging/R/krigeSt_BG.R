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

data <- co_prec[,"PRCP", drop=FALSE]

stsdf <- STIDF(as(co_locs.sp, "SpatialPoints")[!is.na(data)], co_time[!is.na(data)], data[!is.na(data),,drop=F])
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
                                          stAni=1),
                          lower=c(0,0.01,0,
                                  0,0.01,0,
                                  0,0.01,0,
                                  0.05),
                          control=list(parscale=c(1,1e3,1,
                                                  1,1e3,1,
                                                  1,1e3,1,
                                                  1)))

attr(smmFit, "MSE")
plot(empStVgm, smmFit, wireframe=TRUE)

plot(empStVgm, list(metFit, smmFit), all=T, wireframe=F)
plot(empStVgm, list(metFit, smmFit), all=T, wireframe=T)
plot(empStVgm, list(metFit, smmFit), diff=T, wireframe=T)

plot(empStVgm, wireframe=T)
plot(empStVgm, smmFit, wireframe=T)

## interpolation
# re-scale
smmFit$space$range <- smmFit$space$range*1e3
smmFit$joint$range <- smmFit$joint$range*1e3
smmFit$stAni <- smmFit$stAni*1e3
  
predST <- krigeST(PRCP~1, stsdf[,1:10], STF(co_grids, time = stsdf@time[1:10]), smmFit)

stplot(predST)
stplot(stsdf[,1:11,"PRCP"], number=10)

# LOOCV:
stationId <- 7
boolSel <- stsdf@index[,1] == stationId

stopifnot(sum(boolSel)>0)

pred_loocv <- krigeST(PRCP~1, stsdf[!boolSel,,drop=F], stsdf[boolSel,,drop=F], smmFit, nmax=10)

pred_loocv@data$var1.pred[pred_loocv@data$var1.pred < 0] <- 0

library(lattice)
xyplot(var1.pred~PRCP, pred_loocv@data)



## some rather crazy ideas using a marginal tansform
hist(stsdf@data$PRCP)

stsdf@data$rnkPRCP <- (rank(stsdf@data$PRCP, ties.method = "min")-1)/length(stsdf@data$PRCP)
stsdf@data$rnkPRCP[stsdf@data$rnkPRCP == 0] <- 0.5
stsdf@data$rnkPRCP <- qnorm(stsdf@data$rnkPRCP)
hist(stsdf@data$rnkPRCP)

stplot(stsdf[,1:10,"rnkPRCP"])

empStVgmRnk <- variogramST(rnkPRCP~1, stsdf, tlags = 0:3, width=30000/10)
plot(empStVgmRnk, wireframe=TRUE, scales=list(arrows=F))


# rescale to ease optim
empStVgmRnk$dist <- empStVgmRnk$dist/1e3
empStVgmRnk$avgDist <- empStVgmRnk$avgDist/1e3

metFit <- fit.StVariogram(empStVgmRnk, 
                          vgmST("metric", joint=vgm(0.35, "Sph", 30, 0.05), stAni=1))
attr(metFit, "MSE")

smmFit <- fit.StVariogram(empStVgmRnk, 
                          vgmST("sumMetric", 
                                space=vgm(0.15, "Sph", 30, 0.01),
                                time=vgm(0.15, "Sph", 60, 0.01),
                                joint=vgm(0.15, "Sph", 30, 0.01),
                                stAni=1),
                          lower=c(0,0.01,0,
                                  0,0.01,0,
                                  0,0.01,0,
                                  0.05),
                          control=list(parscale=c(1,1e3,1,
                                                  1,1e3,1,
                                                  1,1e3,1,
                                                  1)))
attr(smmFit, "MSE")
plot(empStVgmRnk, smmFit, wireframe=TRUE)

plot(empStVgmRnk, list(metFit, smmFit), all=T, wireframe=T, diff=F, scales=list(arrows=F))
plot(empStVgmRnk, list(metFit, smmFit), all=T, wireframe=T, diff=T, scales=list(arrows=F))

## interpolation
# re-scale
smmFit$space$range <- smmFit$space$range*1e3
smmFit$joint$range <- smmFit$joint$range*1e3
smmFit$stAni <- smmFit$stAni*1e3

predST <- krigeST(rnkPRCP~1, stsdf[,1:10], STF(co_grids, time = stsdf@time[1:10]), smmFit, nmax=10)

stplot(predST, sp.layout=list(spatial.points=stsdf@sp))
stplot(stsdf[,1:11,"rnkPRCP"], number=10)

hist(pnorm(predST@data$var1.pred))

hist(predST@data$var1.pred)


# LOOCV:
stationId <- 1
boolSel <- stsdf@index[,1] == stationId

stopifnot(sum(boolSel)>0)

pred_loocv <- krigeST(PRCP~1, stsdf[!boolSel,,drop=F], stsdf[boolSel,,drop=F], smmFit, nmax=5)

library(lattice)
xyplot(var1.pred~PRCP, pred_loocv@data)

