## Preparation of the geochem data set:
## http://spatial-analyst.net/book/USgrids5km

library(rgdal)
library(utils)

download.file("http://spatial-analyst.net/book/system/files/usgrids5km.zip", "/home/tom/Downloads/usgrids5km.zip")
unzip("/home/tom/Downloads/usgrids5km.zip")
sapply(c("geomap.asc","globedem.asc","dTRI.asc","nlights03.asc","dairp.asc","sdroads.asc"), function(i){system(paste0('gdalwarp /home/tom/Downloads/usgrids5km/', i,' ', gsub(".asc", ".tif", i),' -te 360000 1555000 985000 2210000 -r \"near\" -co \"COMPRESS=DEFLATE\"'))})
usa5km = stack(gsub(".asc", ".tif", c("geomap.asc","globedem.asc","dTRI.asc","nlights03.asc","dairp.asc","sdroads.asc")))
usa5km = as(usa5km, "SpatialGridDataFrame")
proj4string(usa5km) = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
usa5km$geomap = as.factor(usa5km$geomap)
usa5km = as(usa5km, "SpatialPixelsDataFrame")
usa5km = usa5km[!is.na(usa5km$geomap),]
saveRDS(usa5km, "usa5km.rds")

download.file("https://mrdata.usgs.gov/geochem/ngs.zip", "/home/tom/Downloads/ngs.zip")
unzip("/home/tom/Downloads/ngs.zip")
geochem = readOGR("/home/tom/Downloads/ngs/ngs.shp", "ngs")
ov = over(x=spTransform(geochem, CRS(proj4string(usa5km))), y=usa5km)
geochem = geochem[!is.na(ov$geomap),]
geochem = as.data.frame(geochem[,c("LABNO","PB_ICP40","CU_ICP40","CD_ICP40","AS_ICP40","MG_ICP40","K_ICP40","TYPEDESC","CATEGORY")])
saveRDS(geochem, "geochem.rds")
