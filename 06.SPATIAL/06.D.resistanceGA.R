library(ResistanceGA)
library(GA)
library(gdistance)
library(scales)
rm(list = ls())
aggregationFactor = 1
dosmallraster = T
doabund = T
outputdirectory = "~/Distances/"
setwd(outputdirectory)
if (dosmallraster == T) {
  samples = read.csv("LATLONGS_BY_SPECIES.txt",sep="\t")
  allspp = as.character(sort(unique(samples$spp)))
  whereraster="~/3.enms/"
  for (j in 1:length(allspp)) {
    print(j)
    spp = allspp[j]
    print(spp)
    pathlist = list.files(path=whereraster,
                          pattern=glob2rx(paste("*",tolower(spp),"*rescaled*.asc",sep="")),
                          recursive = T)
    pathlist = pathlist[lapply(pathlist,function(x) length(grep("USA",x,value=FALSE))) == 0]
    pathlist = pathlist[lapply(pathlist,function(x) length(grep("zoomed",x,value=FALSE))) == 0]
    sample.locales = (samples[samples$spp==spp,c(3,2)])
    sample.locales = unique(sample.locales)
    sample.locales = sample.locales[complete.cases(sample.locales),]
    maxlat = max(sample.locales$lat,na.rm = T)+1
    maxlong = max(sample.locales$long,na.rm = T)+1
    minlat = min(sample.locales$lat,na.rm = T)-1
    minlong = min(sample.locales$long,na.rm = T)-1
    print(paste(minlong,maxlong,minlat,maxlat))
    print("setting extent")
    ext = extent(minlong,maxlong,minlat,maxlat)
    print(ext)
    sample.locales = SpatialPoints(sample.locales)
    for (i in 1:length(pathlist)) {
      print(paste(i,"of",length(pathlist)))
      old <- Sys.time()
      print(old)
      rasterpath = paste(whereraster,pathlist[i],sep="")
      print(rasterpath)
      splitpath = strsplit(rasterpath,"/")[[1]]
      suffix = splitpath[length(splitpath)]
      continuous = raster(rasterpath)
      ## crop raster
      print("cropping")
      continuous = crop(continuous,ext)
      png(paste(outputdirectory,suffix,".AGGREGATE-",aggregationFactor,".png",sep=""))
      par(mfrow=c(1,2))
      plot(continuous)
      if (aggregationFactor > 1) {
        print(paste("aggregating by factor of",aggregationFactor))
        continuous <- aggregate(continuous, fact=aggregationFactor, fun=mean)
      } else {
        print("skipping aggregation")
      }
      plot(continuous)
      points(sample.locales)
      dev.off()
      temp = sample.locales
      temp = SpatialPoints(temp)
      temp = temp[complete.cases(extract(continuous,temp)),]
      individuals = as.numeric(rownames(temp@coords))
      included = samples[individuals,]
      latlong = paste(included$lat,included$long)
      rawvalues = as.data.frame(extract(continuous,temp))
      rownames(rawvalues) = latlong
      rawdist = as.data.frame(as.matrix(dist(rawvalues)))
      rawfile = paste(outputdirectory,suffix,"AGGFACT-",aggregationFactor,".RAW-EUCLIDIAN.txt",sep="")
      write.csv(rawdist,rawfile)
      outfile = paste(outputdirectory,suffix,"AGGFACT-",aggregationFactor,".MORPH-AND-GENE-DISTANCES.txt",sep="")
      gdist.inputs <- gdist.prep(length(temp),
                                 samples = temp,
                                 method = 'costDistance') 
      print("gdist response")
      gdist.response <- Run_gdistance(gdist.inputs = gdist.inputs,
                                      r = continuous,
                                      scl=F) 
      distdf = as.matrix(gdist.response)
      rownames(distdf) = latlong
      colnames(distdf) = latlong
      write.csv(distdf,outfile)
      new <- Sys.time() - old # calculate difference
      print(paste("Time elapsed:",new))
    }
  }
}
if (doabund == T) {
  samples = read.csv("LATLONGS_ONLY.txt",sep="\t")
  whereraster="~/abundance_asciis/"
  pathlist = list.files(path=whereraster,
                        pattern=glob2rx(paste("*.asc",sep="")),
                        recursive = T)
  sample.locales = (samples[,c("long","lat")])
  sample.locales = unique(sample.locales)
  sample.locales = sample.locales[complete.cases(sample.locales),]
  maxlat = max(sample.locales$lat,na.rm = T)+1
  maxlong = max(sample.locales$long,na.rm = T)+1
  minlat = min(sample.locales$lat,na.rm = T)-1
  minlong = min(sample.locales$long,na.rm = T)-1
  sample.locales = unique(sample.locales)
  sample.locales = SpatialPoints(sample.locales)
  print("setting extent")
  ext = extent(minlong,maxlong,minlat,maxlat)
  print(ext)
  for (i in 1:length(pathlist)) {
    print(paste(i,"of",length(pathlist)))
    old <- Sys.time()
    print(old)
    rasterpath = paste(whereraster,pathlist[i],sep="")
    print(rasterpath)
    splitpath = strsplit(rasterpath,"/")[[1]]
    suffix = splitpath[length(splitpath)]
    continuous = raster(rasterpath)
    ## crop raster
    print("cropping")
    continuous = crop(continuous,ext)
    png(paste(outputdirectory,suffix,".AGGREGATE-",aggregationFactor,".png",sep=""))
    par(mfrow=c(1,2))
    plot(continuous[[1]])
    ## need to rescale the individual layers 
    print("rescaling")
    todo = continuous
    vals = values(todo)
    vals2 = rescale(vals)
    values(todo) = vals2
    continuous = todo
    if (aggregationFactor > 1) {
      print(paste("aggregating by factor of",aggregationFactor))
      continuous <- aggregate(continuous, fact=aggregationFactor, fun=mean)
    } else {
      print("skipping aggregation")
    }
    plot(continuous[[1]])
    points(sample.locales)
    dev.off()
    ## need to transform the abundances by Inverse Monomolecular (7)
    r.tran <- Resistance.tran(transformation = 7,
                              shape = 3.5,
                              max = 150,
                              r = continuous)
    values(r.tran) = rescale(values(r.tran))
    png(paste(outputdirectory,suffix,"AGGFACT-",aggregationFactor,".FULL_RESISTANCE_SURFACE.png",sep=""))
    plot(r.tran)
    dev.off()
    writeRaster(r.tran,paste(outputdirectory,suffix,"AGGFACT-",aggregationFactor,".FULL_RESISTANCE_SURFACE.asc",sep=""),
                format="ascii")
    temp = sample.locales
    temp = SpatialPoints(temp)
    temp = temp[complete.cases(extract(continuous,temp)),]
    individuals = as.numeric(rownames(temp@coords))
    included = samples[individuals,]
    latlong = paste(included$lat,included$long)
    ## get euclidian distances
    rawvalues = as.data.frame(extract(continuous,temp))
    rownames(rawvalues) = latlong
    rawdist = as.data.frame(as.matrix(dist(rawvalues)))
    rawfile = paste(outputdirectory,suffix,"AGGFACT-",aggregationFactor,".RAW-EUCLIDIAN.txt",sep="")
    write.csv(rawdist,rawfile)
    outfile = paste(outputdirectory,suffix,"AGGFACT-",aggregationFactor,".MORPH-AND-GENE-DISTANCES.txt",sep="")
    gdist.inputs <- gdist.prep(length(temp),
                               samples = temp,
                               method = 'costDistance') 
    print("gdist response")
    gdist.response <- Run_gdistance(gdist.inputs = gdist.inputs,
                                    r=r.tran,
                                    scl=F) 
    distdf = as.matrix(gdist.response)
    rownames(distdf) = latlong
    colnames(distdf) = latlong
    write.csv(distdf,outfile)
    new <- Sys.time() - old # calculate difference
    print(paste("Time elapsed:",new))
  }
}
