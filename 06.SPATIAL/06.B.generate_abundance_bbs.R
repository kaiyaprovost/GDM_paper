library(tidyr)
library(gstat)
library(raster)
setwd("~/BBS_data/")
bbsdata = "Both_fulldata_onlybig10.csv"
bbs = read.csv(bbsdata)
names(bbs)
agg = aggregate(bbs$SPECIESTOTAL~bbs$LATITUDE+bbs$LONGITUDE+bbs$SPECIES,FUN=function(x) {mean(x,na.rm=T)}) 
names(agg) = c("LAT","LON","SPP","NUM")
completeagg = (as.data.frame(complete(agg,SPP,nesting(LAT,LON),fill=list(NUM=0))))
head(completeagg)
nrow(completeagg)
agg = completeagg
agg$NUM[agg$NUM==0] = 0.0001
ramp = colorRamp(c("blue","yellow", "red"),alpha=T)
colors <- ramp(log10(agg$NUM)/log10(1000))
agg$COL = colors/255
agg$NUM[agg$NUM==0.0001] = 0
agg$COL[is.na(agg$COL)] = 0
head(agg)
for (spp in unique(agg$SPP)) {
  png(paste(spp,"_abundance.png",sep=""))
  layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
  print(spp)
  temp = agg[agg$SPP==spp,]
  plot(agg$LON,agg$LAT,col="lightgrey",pch=0,main=spp,cex=0.5,
       ylab="lat",xlab="long")
  points(temp$LON,temp$LAT,col=rgb(temp$COL[,1],temp$COL[,2],temp$COL[,3],temp$COL[,4]),pch=15)
  legend_image <- as.raster(matrix(rgb(ramp((20:1)/20),maxColorValue = 255), ncol=1))
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Log Mean Abundance')
  rasterImage(legend_image, 0, 0, 1,1)
  text(x=1.5, y = seq(0,1,l=4), labels = c(1,10,100,1000))
  dev.off()
}
## set up the raster stuff for interpolating 
ras = raster("AddedPredThin_BestModel_Amphispiza bilineata_addedLayers_worldclim.asc_-6.46640157699585_-69.5545196533203.rescaled.asc")
values(ras)[values(ras) >= 0] = 1
x_range <- as.numeric(c(extent(ras)[1], extent(ras)[2]))  # min/max longitude of the interpolation area
y_range <- as.numeric(c(extent(ras)[3], extent(ras)[4]))  # min/max latitude of the interpolation area
# create an empty grid of values ranging from the xmin-xmax, ymin-ymax
## convert the grd by an order of mag
grd <- expand.grid(x = seq(from = x_range[1],
                           to = x_range[2], 
                           by = res(ras)[1]*10),
                   y = seq(from = y_range[1],                                           
                           to = y_range[2], 
                           by = res(ras)[2]*10))  # expand points to grid
class(grd)
# Convert grd object to a matrix and then turn into a spatial
# points object
coordinates(grd) <- ~x + y
# turn into a spatial pixels object
gridded(grd) <- TRUE
for (spp in unique(agg$SPP)) { 
  print(spp)
  temp = agg[agg$SPP==spp,]
  coordinates(temp) <- ~LON + LAT
  idw_pow1 <- idw(formula = NUM ~ 1,
                  locations = temp,
                  newdata = grd,
                  idp = 5)
  idw = (raster(idw_pow1, layer=1, values=TRUE))
  writeRaster(idw,paste(spp,"_Idw_abundance.asc",sep=""),format="ascii",overwrite=T)
  idw2 = resample(idw,ras)
  extent(idw2)
  extent(ras)
  cropped = idw2*ras
  writeRaster(cropped,paste(spp,"_Idw_abundance_clipped.asc",sep=""),format="ascii",overwrite=T)
}
files = list.files(pattern=".asc$")
for (file in files) {
  print(file)
  spp = strsplit(file,"_")[[1]][1]
  ras = raster(file)
  png(paste(file,".png",sep=""))
  plot(ras,main=spp)
  dev.off()
}
