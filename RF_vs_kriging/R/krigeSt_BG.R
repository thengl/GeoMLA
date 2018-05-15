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

co_locs.sp@data$PRISM_prec <- over(co_locs.sp, co_grids[c("elev_1km","PRISM_prec")])$PRISM_prec

co_time <- as.POSIXct(co_prec$DATE)

data <- co_prec[,c("PRCP", "cdate", "doy","ELEVATION"), drop=FALSE]
data$cosDoy <- cos((data$doy)*pi/365)

# summary(lm(PRCP ~ ELEVATION + cosDoy + cdate, data))
# round(cor(data[sample(176467,1e5),], use = "p")[,4],3)

stsdf <- STIDF(as(co_locs.sp, "SpatialPointsDataFrame")[!is.na(data$PRCP),], co_time[!is.na(data$PRCP)], data[!is.na(data$PRCP),,drop=F])
stsdf <- as(stsdf, "STSDF")

# df <- as.data.frame(stsdf)
# summary(lm(PRCP ~ cdate + doy + PRISM_prec, df[]))
# 
# stsdf@data$PRISM_prec <- df$PRISM_prec
# stsdf@data$diff_PRCP <- stsdf@data$PRISM_prec - stsdf@data$PRCP
#   
# spplot(co_grids, "PRISM_prec")
# 
# round(cor(df[,-c(3:7)], use = "p", method = "kendall")[,4],2)


## no covariates - no zero correction
empStVgm <- gstat::variogramST(PRCP~1, stsdf, tlags = 0:3)
plot(empStVgm, wireframe=TRUE, scales=list(arrows=F))
plot(empStVgm, wireframe=FALSE, scales=list(arrows=F))

# ## no covariates - no zero correction
# empStVgmNoZero <- gstat::variogramST(diff_PRCP~1, stsdf[stsdf@data$PRCP > 0,,drop=F], tlags = 0:3, boundaries=c(0,2000,5000,10000,15000, 20e3, 30e3,40e3), na.omit = TRUE)
# empStVgmNoZero
# # <- empStVgmNoZero[!is.nan(empStVgmNoZero$avgDist),]
# empStVgmNoZero$avgDist
# 
# plot(empStVgmNoZero, wireframe=T, scales=list(arrows=F))
# plot(empStVgmNoZero, wireframe=F, scales=list(arrows=F))


## fit of st vgm

# rescale to ease optim
empStVgm$dist <- empStVgm$dist/1e3
empStVgm$avgDist <- empStVgm$avgDist/1e3

metFit <- fit.StVariogram(empStVgm, vgmST("metric", joint=vgm(0.035, "Sph", 120, 0.005), stAni=1.25))

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

attr(smmFit, "MSE")*1e6
plot(empStVgm, smmFit, wireframe=TRUE)
plot(empStVgm, list(metFit, smmFit), all=T, wireframe=F)
plot(empStVgm, list(metFit, smmFit), diff=T, wireframe=T)

plot(empStVgm, wireframe=T)
plot(empStVgm, smmFit, wireframe=T)


png("RF_vs_kriging/results/st_prec/ST-vgm.png", width = 2600, height = 1000, res = 200)
plot(empStVgm, list(metFit, smmFit), all=T, wireframe=T, scales=list(arrows=FALSE))
dev.off()

## interpolation
# re-scale
smmFit$space$range <- smmFit$space$range*1e3
smmFit$joint$range <- smmFit$joint$range*1e3
smmFit$stAni <- smmFit$stAni*1e3
  
predST <- krigeST(PRCP~1, stsdf[,818:833], STF(co_grids, time = stsdf@time[823:828]), 
                  smmFit, nmax = 10, computeVar = TRUE)

predST@data$var1.sd <- sqrt(predST@data$var1.var) 

library(rgdal)
T.lst = paste0("2016-02-0", c(1:6))
for (i in 1:6) {
  writeGDAL(predST[,i, "var1.pred"], 
            fname=paste0("RF_vs_kriging/results/st_prec/krigeSt_PRCP_", T.lst[i], ".tif"),
            options="COMPRESS=DEFLATE", type = "Int16", mvFlag = "-32768")
  writeGDAL(predST[,i, "var1.sd"], 
            fname=paste0("RF_vs_kriging/results/st_prec/krigeSt_PRCP_", T.lst[i], "_pe.tif"),
            options="COMPRESS=DEFLATE", type = "Int16", mvFlag = "-32768")
}


png("RF_vs_kriging/results/st_prec/precip_kriged.png", width = 2000, height = 2000, res = 200)
stplot(predST[,,"var1.pred"], main="Kriged Precipitations")
dev.off()

png("RF_vs_kriging/results/st_prec/precip_kriged_sd.png", width = 2000, height = 2000, res = 200)
stplot(predST[,,"var1.sd"], main="Kriging standard deviation")
dev.off()


colVec <- colorRampPalette(c("white","lightblue","blue","darkblue"))(100)
pointCol <- function(x) colVec[findInterval(x, seq(0,1, length.out = 100))]

for (i in 1:6) {
  ptSub <- stsdf[,(823:828)[i],"PRCP"]
  png(paste0("RF_vs_kriging/results/st_prec/precip_kriged_", T.lst[i], ".png"),
      width = 1000, height = 1000, res = 200)
  print(spplot(predST[,i], "var1.pred", col.regions=colVec, 
               at=seq(0,1, length.out = 100),
               main=paste("Kriged precipitation at:", T.lst[i]),
               sp.layout=list(list(sp.points=ptSub, pch=21, col="darkgrey", cex=1.4),
                              list(sp.points=ptSub, pch=16, col=pointCol(ptSub@data$PRCP), cex=1.2 ) ),
               colorkey=list(colVec, at=seq(0,1,length.out = 100))))
  dev.off()
}


# LOOCV:
stsdf@data$loocv_pred <- NA
nSp <- length(stsdf@sp$STATION)

for (stationId in 221:nSp) {
  cat("Started to predict station",stationId,"just now. \n") 
  timing <- Sys.time()
  pred <- krigeST(PRCP~1, stsdf[(1:nSp)[-stationId],,drop=F], stsdf[stationId,,drop=F], smmFit, nmax=10)@data$var1.pred
  
  pred[pred < 0] <- 0
  
  # pred_loocv[[stationId]] <- pred
  predVec <- stsdf@data$loocv_pred
  predVec[stsdf@index[,1]==stationId] <- pred
  stsdf@data$loocv_pred <- predVec
  
  cat("The prediction took:",Sys.time()-timing,"\n") 
}

library(lattice)
plot(stsdf[1,,drop=F]@data$PRCP, stsdf[1,,drop=F]@data$loocv_pred)
abline(0,1)

# MAE:
mean(abs(stsdf@data$PRCP - stsdf@data$loocv_pred), na.rm = T) 

# RMSE:
sqrt(mean((stsdf@data$PRCP - stsdf@data$loocv_pred)^2, na.rm = T)) # 0.069

# ME/bias
mean(stsdf@data$PRCP - stsdf@data$loocv_pred, na.rm = T) # 0.00005

# COR
cor(stsdf@data$PRCP, stsdf@data$loocv_pred, use = "p") # 0.930

# MAE:  0.022
# RMSE: 0.069
# ME:  <0.000
# COR:  0.930

hist(stsdf@data$PRCP)
hist(stsdf@data$loocv_pred)
