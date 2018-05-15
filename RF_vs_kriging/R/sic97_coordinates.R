## http://www.meteoschweiz.admin.ch/home/wetter/messwerte/messwerte-an-stationen.html?
laea.prj = "+proj=laea +lat_0=46.95240555555556 +lon_0=7.439583333333333 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
ch.prj = "+proj=somerc +lat_0=46.95240555555556 +lon_0=7.439583333333333 +k_0=1 +x_0=600000 +y_0=200000 +ellps=bessel +towgs84=674.4,15.1,405.3,0,0,0,0 +units=m +no_defs"

ch.st = read.csv("data/rainfall/CH_stations.csv")
coordinates(ch.st) = ~ X+Y
proj4string(ch.st) = ch.prj
shape = "http://maps.google.com/mapfiles/kml/pal2/icon18.png" 
kml(ch.st, colour=Altitude, shape=shape)
XY = spTransform(ch.st, CRS(laea.prj))
ch.st$X0 = XY@coords[,1]
ch.st$Y0 = XY@coords[,2]

summary(sic_full@coords)
#proj4string(sic_full) = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
proj4string(sic_full) = "+proj=laea +lat_0=50 +lon_0=11 +x_0=5000000 +y_0=3200000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
kml(sic_full, colour=rainfall)

data(SIC)
#proj4string(swissAltitude)
#plot(swissAltitude, main="elevation")
#points(swissRain)
#plot(swissBorder, add=TRUE)
#plot(swissLandType)
s.x = sic.all$coords[,"V2"]*1000 - 17791.29 + 2672591
s.y = sic.all$coords[,"V3"]*1000 - 13224.66 + 1200225 
sic.sp = SpatialPointsDataFrame(coords=data.frame(s.x, s.y), data = data.frame(rainfall=sic.all$data))
proj4string(sic.sp) = CRSargs(CRS("+init=epsg:2056"))

data("swissRain")
laea.prj = "+proj=laea +lat_0=46.95240555555556 +lon_0=7.439583333333333 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" 
swissRain = read.table("sic_full.dat", sep=',', col.names=c('ID','x','y','rainfall'))
## Real pity that the coordinate system for the SIC 1997 exercise is still unknown and need to be guessed...
swissRain$X = swissRain$x + 54925
swissRain$Y = swissRain$y - 13264
sic97.sp = SpatialPointsDataFrame(coords=swissRain[,c("X","Y")], data = swissRain["rainfall"])
proj4string(sic97.sp) = laea.prj
#kml(sic97.sp, colour=rainfall, shape=shape)
ch.dem = reproject(swissAltitude, laea.prj)
#plot(ch.dem)
ch.pol = spTransform(swissBorder, CRS(laea.prj))
system(paste0('gdalwarp /mnt/cartman/MERIT/MERIT_100m.vrt swissDEM.tif -t_srs \"', laea.prj,'\" -te -133321.8 -153663.2 259108.2 123836.8 -tr 1000 1000 -co \"COMPRESS=DEFLATE\"'))
system(paste0('gdalwarp /mnt/cartman/CHELSA/CHELSA_prec_5_V1.2_land.tif swissChelsa_May.tif -t_srs \"', laea.prj,'\" -te -133321.8 -153663.2 259108.2 123836.8 -tr 1000 1000 -co \"COMPRESS=DEFLATE\"'))
ch1km = readGDAL("data/rainfall/swissChelsa_May.tif")
ch1km$DEM = readGDAL("data/rainfall/swissDEM.tif")$band1
ch.bor = rasterize(ch.pol, raster(ch1km))
ch1km$border = as(ch.bor, "SpatialGridDataFrame")$NAME_ENGLISH
ch1km = as(ch1km, "SpatialPixelsDataFrame")
ch1km = ch1km[!is.na(ch1km$border),]
str(ch1km@data)
names(ch1km)[1] = "CHELSA_rainfall"
saveRDS(ch1km, "data/rainfall/swiss1km.rds")
saveRDS(sic97.sp, "data/rainfall/sic97.rds")
