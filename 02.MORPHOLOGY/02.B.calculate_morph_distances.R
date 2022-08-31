library(gdm)
library(sgdm)
library(raster)
library(rgdal)
library(vegan)
library(ecodist)
doPCA = T
morph = "KLP_Master_Spreadsheet_Morphology_AGGREGATED_19_July_2019.csv"
morphdf = read.csv(morph)
morphdf = morphdf[,2:40]
morphlocs = morphdf[,c("LAT","LONG")]
morphdata = morphdf[,c(3:9)]
morphspp = morphdf$SPP
setwd("")
if (doPCA==T) {
  ## pca morph 
  morphpca = "KLP_Master_Spreadsheet_Morphology_AGGREGATED_19_July_2019_PCASCORES_FULL.csv"
  pcadf = read.csv(morphpca)
  pcalocs = pcadf[,c("LAT","LONG")]
  pcadata = pcadf[,c(1:3)]
  pcaspp = pcadf$SPP
}
## pure environments
envfile = "ENMS_multilayer_onlyworldclim.tif"
env = stack(envfile)
plot(env[[1]])
points(morphlocs$LONG,morphlocs$LAT)
envdata = as.data.frame(extract(env,morphlocs[,2:1]))
envdata$SPP = morphspp
envdata$CATALOG.NUMBER = morphdf$CATALOG.NUMBER
inputdata = merge(morphdf,envdata)
## try a mantel test between matrices
plotPairwiseLM = function(spp="BELLII",
                          pngsuffix="morph_geog_dist.png",
                          data=merged,
                          xcol="GEOGDIST",
                          ycol="MORPHDIST",
                          xlab=xcol,
                          ylab=ycol,
                          lmtype="linear") {
  png(paste(spp,pngsuffix,sep="_"))
  data_x = data[,xcol]
  data_y = data[,ycol]
  final = cbind(data_x,data_y)
  tokeep = which(complete.cases(final))
  data_x = data_x[tokeep]
  data_y = data_y[tokeep]
  if(length(unique(data_x)) > 1 & length(unique(data_y)) > 1) {
    if(lmtype=="linear") {
      mod=lm(data_y~data_x,na.action=na.exclude)
      rsq=summary(mod)$adj.r.squared
      pval=summary(mod)$coefficients[2,4]
      newx <- seq(min(data_x), max(data_x), length.out=100)
      preds <- predict(mod, newdata=data.frame(data_x=newx),
                       interval = 'confidence')
      plot(data_x,data_y,
           col=rgb(0,0,0,1),
           xlab=xlab,
           ylab=ylab,
           main=paste(spp,"\nR^2: ",round(rsq,2),"\np: ",formatC(pval, format = "e", digits = 2),sep=""))
      polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = rgb(1,0,0,0.3), border = NA)
      lines(newx, preds[ ,3], lty = 'dashed', col = 'red')
      lines(newx, preds[ ,2], lty = 'dashed', col = 'red')
      abline(mod,col="red")
    } else if (lmtype=="logistic") {
      data_y = scales::rescale(data_y)
      data_y[data_y>=0.5] = 1
      data_y[data_y<0.5] = 0
      mod=glm(data_y~data_x,na.action=na.exclude,family="binomial")
      aic=AIC(mod)
      pval=summary(mod)$coefficients[2,4]
      newx <- seq(min(data_x), max(data_x), length.out=100)
      plot(data_x,data_y,
           col=rgb(0,0,0,1),
           xlab=xlab,
           ylab=ylab,
           main=paste(spp,"\nAIC: ",round(aic,2),"\np: ",formatC(pval, format = "e", digits = 2),sep=""))
      curve(predict(mod,data.frame(data_x=x),type="resp"),add=TRUE,col="red")
      abline(h=0.5,col="grey",lty="dashed")
    }
  } else {
    pval="N/A"
    rsq="N/A"
    plot(data_x,data_y,
         col=rgb(0,0,0,1),
         xlab=xlab,
         ylab=ylab,
         main=paste(spp,"\nR^2: ",rsq,"\np: ",pval,sep=""))
  }
  dev.off()
}
for (spp in unique(morphspp)) {
  print(spp)
  if (spp != "") {
    locs = morphlocs[morphspp==spp,]
    df = morphdf[morphspp==spp,]
    if (doMorph==T) {
      data = morphdf[morphspp==spp,c(3:9)]
    } else {
      data = df
    }
    envspp = envdata[envdata$SPP==spp,]
    if (unique(morphdf$GENETIC.SPLIT[morphspp==spp]) == 0) {
      structure = morphdf[morphspp==spp,c("GENETIC.SIDE","GENETIC.SPLIT","WHICH.SIDE.OF.CFB")]
      structure$GENETIC.SIDE[structure$WHICH.SIDE.OF.CFB=="SONORAN"] = 0
      structure$GENETIC.SIDE[structure$WHICH.SIDE.OF.CFB=="CHIHUAHUAN"] = 1
      structure = structure[,c("GENETIC.SIDE","GENETIC.SPLIT")]
      structure$GENETIC.SIDE = scales::rescale(structure$GENETIC.SIDE)
    } else {
      structure = morphdf[morphspp==spp,c("GENETIC.SIDE","GENETIC.SPLIT")]
      structure$GENETIC.SIDE = scales::rescale(structure$GENETIC.SIDE)
    }
    strdist = as.data.frame(as.matrix(dist(structure, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)))
    colnames(strdist) = df$CATALOG.NUMBER
    rownames(strdist) = df$CATALOG.NUMBER
    pairwisestr = data.frame( t(combn(names(strdist),2)), dist=t(strdist)[lower.tri(strdist)] )
    names(pairwisestr) = c("FIRST","SECOND","STRDIST")
    geogdist = as.data.frame(as.matrix(dist(locs, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)))
    colnames(geogdist) = df$CATALOG.NUMBER
    rownames(geogdist) = df$CATALOG.NUMBER
    pairwisegeog = data.frame( t(combn(names(geogdist),2)), dist=t(geogdist)[lower.tri(geogdist)] )
    names(pairwisegeog) = c("FIRST","SECOND","GEOGDIST")
    if (doMorph == T) {
      morphdist = as.data.frame(as.matrix(dist(data, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)))
      colnames(morphdist) = df$CATALOG.NUMBER
      rownames(morphdist) = df$CATALOG.NUMBER
      pairwisemorph = data.frame( t(combn(names(morphdist),2)), dist=t(morphdist)[lower.tri(morphdist)] )
      names(pairwisemorph) = c("FIRST","SECOND","MORPHDIST")
    }
    envdist = as.data.frame(as.matrix(dist(envspp, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)))
    colnames(envdist) = df$CATALOG.NUMBER
    rownames(envdist) = df$CATALOG.NUMBER
    pairwiseenv = data.frame( t(combn(names(envdist),2)), dist=t(envdist)[lower.tri(envdist)] )
    names(pairwiseenv) = c("FIRST","SECOND","ENVDIST")
    write.csv(geogdist,paste(spp,"_distancematrix_geog.csv",sep=""))
    write.csv(envdist,paste(spp,"_distancematrix_env.csv",sep=""))
    write.csv(strdist,paste(spp,"_distancematrix_str.csv",sep=""))
    if (doMorph == T) {
      write.csv(morphdist,paste(spp,"_distancematrix_morph.csv",sep=""))
      merged = merge(pairwisegeog,pairwisemorph)
      merged = merge(merged,pairwiseenv)
      merged = merge(merged,pairwisestr)
      write.csv(merged,paste(spp,"_distances.csv",sep=""))
      plotPairwiseLM(spp=spp,pngsuffix="morph_geog_dist.png",
                     data=merged,
                     xcol="GEOGDIST",
                     ycol="MORPHDIST",
                     xlab="geographic distance",
                     ylab="morphological distance")
      plotPairwiseLM(spp=spp,pngsuffix="morph_env_dist.png",
                     data=merged,
                     xcol="ENVDIST",
                     ycol="MORPHDIST",
                     xlab="environmental distance",
                     ylab="morphological distance")
      plotPairwiseLM(spp=spp,pngsuffix="env_geog_dist.png",
                     data=merged,
                     xcol="GEOGDIST",
                     ycol="ENVDIST",
                     xlab="geographic distance",
                     ylab="environmenal distance")
      plotPairwiseLM(spp=spp,pngsuffix="morph_str_dist.png",
                     data=merged,
                     ycol="STRDIST",
                     xcol="MORPHDIST",
                     ylab="structure distance",
                     xlab="morphological distance",
                     lmtype="logistic")
      plotPairwiseLM(spp=spp,pngsuffix="str_env_dist.png",
                     data=merged,
                     xcol="ENVDIST",
                     ycol="STRDIST",
                     xlab="environmental distance",
                     ylab="structure distance",
                     lmtype="logistic")
      plotPairwiseLM(spp=spp,pngsuffix="str_geog_dist.png",
                     data=merged,
                     xcol="GEOGDIST",
                     ycol="STRDIST",
                     xlab="geographic distance",
                     ylab="structure distance",
                     lmtype="logistic")
      forpca = df[,c("BILL.HEIGHT","BILL.LENGTH",
                     "BILL.WIDTH","TARSUS.LENGTH",
                     "WING.LENGTH.TO.PRIMARIES",
                     "WING.LENGTH.TO.SECONDARIES","SPP","CATALOG.NUMBER")]
      geopca = df[complete.cases(forpca),c("LAT","LONG","CATALOG.NUMBER")]
      forpca = forpca[complete.cases(forpca),]
      morphpca = prcomp(forpca[,1:6])
      pcadata = as.data.frame(morphpca$x)
      pcadata$CATALOG.NUMBER = forpca$CATALOG.NUMBER
      pcamerge = merge(geopca,pcadata)
      pca_morphdist = as.data.frame(as.matrix(dist(pcamerge[,4:6], method = "euclidean", diag = FALSE, upper = FALSE, p = 2)))
      colnames(pca_morphdist) = pcamerge$CATALOG.NUMBER
      rownames(pca_morphdist) = pcamerge$CATALOG.NUMBER
      pairwise_pcamorph = data.frame( t(combn(names(pca_morphdist),2)), dist=t(pca_morphdist)[lower.tri(pca_morphdist)] )
      names(pairwise_pcamorph) = c("FIRST","SECOND","MORPHDIST")
      pca_geodist = as.data.frame(as.matrix(dist(pcamerge[,2:3], method = "euclidean", diag = FALSE, upper = FALSE, p = 2)))
      colnames(pca_geodist) = pcamerge$CATALOG.NUMBER
      rownames(pca_geodist) = pcamerge$CATALOG.NUMBER
      pairwise_pcageo = data.frame( t(combn(names(pca_geodist),2)), dist=t(pca_geodist)[lower.tri(pca_geodist)] )
      names(pairwise_pcageo) = c("FIRST","SECOND","GEOGDIST")
      merged_pcadist = merge(pairwise_pcamorph,pairwise_pcageo)
      png(paste(spp,"pcamorph_geog_dist.png",sep="_"))
      data_x = merged_pcadist$GEOGDIST
      data_y = merged_pcadist$MORPHDIST
      mod=lm(data_y~data_x)
      rsq=summary(mod)$adj.r.squared
      pval=summary(mod)$coefficients[2,4]
      newx <- seq(min(data_x), max(data_x), length.out=100)
      preds <- predict(mod, newdata=data.frame(data_x=newx),
                       interval = 'confidence')
      plot(data_x,data_y,
           col=rgb(0,0,0,1),
           xlab="geographic distance",
           ylab="morphological (PCA) distance",
           main=paste(spp,"\nR^2: ",round(rsq,2),"\np: ",formatC(pval, format = "e", digits = 2),sep=""))
      polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = rgb(1,0,0,0.3), border = NA)
      lines(newx, preds[ ,3], lty = 'dashed', col = 'red')
      lines(newx, preds[ ,2], lty = 'dashed', col = 'red')
      abline(mod,col="red")
      dev.off()
    }
  }
}
if (doPCA==T) {
  for (spp in unique(pcaspp)) {
    print(spp)
    locs = pcalocs[pcaspp==spp,]
    df = pcadf[pcaspp==spp,]
    pc1 = df[,"PC1"]
    pc2 = df[,"PC2"]
    pc3 = df[,"PC3"]
    pc1dist = as.data.frame(as.matrix(dist(pc1, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)))
    colnames(pc1dist) = df$CATALOG.NUMBER
    rownames(pc1dist) = df$CATALOG.NUMBER
    pairwisepc1 = data.frame( t(combn(names(pc1dist),2)), dist=t(pc1dist)[lower.tri(pc1dist)] )
    names(pairwisepc1) = c("FIRST","SECOND","PC1DIST")
    pc2dist = as.data.frame(as.matrix(dist(pc2, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)))
    colnames(pc2dist) = df$CATALOG.NUMBER
    rownames(pc2dist) = df$CATALOG.NUMBER
    pairwisepc2 = data.frame( t(combn(names(pc2dist),2)), dist=t(pc2dist)[lower.tri(pc2dist)] )
    names(pairwisepc2) = c("FIRST","SECOND","PC2DIST")
    pc3dist = as.data.frame(as.matrix(dist(pc3, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)))
    colnames(pc3dist) = df$CATALOG.NUMBER
    rownames(pc3dist) = df$CATALOG.NUMBER
    pairwisepc3 = data.frame( t(combn(names(pc3dist),2)), dist=t(pc3dist)[lower.tri(pc3dist)] )
    names(pairwisepc3) = c("FIRST","SECOND","PC3DIST")
    pairwisepca = merge(pairwisepc1,pairwisepc2)
    pairwisepca = merge(pairwisepca,pairwisepc3)
    write.csv(pc1dist,paste(spp,"_distancematrix_pc1morph.csv",sep=""))
    write.csv(pc2dist,paste(spp,"_distancematrix_pc2morph.csv",sep=""))
    write.csv(pc3dist,paste(spp,"_distancematrix_pc3morph.csv",sep=""))
    write.csv(pairwisepca,paste(spp,"_distances_pairwisepca.csv",sep=""))
  }
}
