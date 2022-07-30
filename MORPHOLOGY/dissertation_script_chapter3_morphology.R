dynamic_require <- function(package) {
  if (eval(parse(text = paste("require(", package, ")"))))
    return(TRUE)
  install.packages(package)
  return(eval(parse(text = paste(
    "require(", package, ")"
  ))))
}
packages = c("PopGenome", "moments", "R.utils", "gtools","progress", "RColorBrewer",
             "corrplot", "devtools", "ggbiplot", "dabestr", "factoextra", "cluster",
             "tidyverse","caret","MASS","plotly","lme4")
for (p in packages) {
  dynamic_require(p)
}
#install_github("vqv/ggbiplot",force=T)

generateData=F
plotResid=F

## perform aggregation
if(generateData==T){
print("AGGREGATE"); {
  path = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/"
  setwd(path)
  #filename = "/Users/kprovost/Dropbox (AMNH)/Dissertation/Stephanie/morphology/KLP_Master_Spreadsheet_Morphology_19june2018.csv"
  filename = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/measurements/KLP_Master_Spreadsheet_Morphology_18_July_2019.csv"
  
  data <- read.delim(filename,sep=",",header=T,stringsAsFactors = F,fill=T)
  
  #plot(data[,c(11:17)])
  data$BILL.LENGTH = as.numeric(data$BILL.LENGTH)
  data$BILL.WIDTH = as.numeric(data$BILL.WIDTH)
  data$BILL.HEIGHT = as.numeric(data$BILL.HEIGHT)
  data$TARSUS.LENGTH = as.numeric(data$TARSUS.LENGTH)
  data$WING.LENGTH.TO.PRIMARIES = as.numeric(data$WING.LENGTH.TO.PRIMARIES)
  data$WING.LENGTH.TO.SECONDARIES = as.numeric(data$WING.LENGTH.TO.SECONDARIES)
  data$TAIL = as.numeric(data$TAIL)
  data$LAT = as.numeric(data$LAT)
  data$LONG = as.numeric(data$LONG)
  
  data$KIPPSINDEX = as.numeric(data$KIPPSINDEX)
  data$BEAKBASEAREA = as.numeric(data$BEAKBASEAREA)
  data$BEAKVOL = as.numeric(data$BEAKVOL)
  data$BEAKLATERALSURFACE = as.numeric(data$BEAKLATERALSURFACE)
  data$BEAKTOTALSURFACE = as.numeric(data$BEAKTOTALSURFACE)
  data$BEAKAREAVOLUME = as.numeric(data$BEAKAREAVOLUME)
  
  #data$GENETIC.SIDE[data$GENETIC.SIDE==2] = 0
  #data$GENETIC.SIDE[data$GENETIC.SIDE==3] = -1
  data$GENETIC.SIDE = as.numeric(data$GENETIC.SIDE)
  
  #data$GENETIC.SPLIT[data$GENETIC.SPLIT=="?"] = 0.5
  data$GENETIC.SPLIT = as.numeric(data$GENETIC.SPLIT)
  
  data$SEQUENCED. = as.numeric(data$SEQUENCED.)
  
  
  
  #plot(data[data$SPP=="FUSCA",c(11:17)])
  #data[data$SPP=="FUSCA",]
  #correlation = cor(data[,c(11:17)],use="pairwise.complete.obs")
  #corrplot(correlation,method="number")
  
  ## notes on correcting for body mass:
  ## To correct for allometric scaling of beak size, we included
  ## body mass as a covariate in each analysis. This approach is
  ## preferred among contemporary phylogenetic comparative studies 
  ## because the use of residuals may cause collinearity issues
  # (Freckleton 2009; Symonds and Tattersall 2010; Baab et al. 2014;
  # Benson-Amram et al. 2015).
  
  ## Morphological variables and
  # body mass were log transformed before analysis. Since
  # nonnormality was not a problem with latitude, altitude,
  # and Tmin, these were untransformed, except in the case of
  # the group analyses of gulls and terns, where Nudds and
  # Oswald’s log-transformed measure of latitude was used.
  # To control for body mass, we calculated multiple regres-
  # sions in COMPARE, predicting bill length (or leg-element
  # length), with both latitude (or altitude or minimum tem-
  # perature) and body mass as predictors.
  
  #nobad = morph[morph$CONDITION=="",]
  #noblank1 = nobad[nobad$BEAKAREAVOLUME!="",]
  #noblank2 = noblank1[!(is.na(noblank1$BEAKAREAVOLUME)),]
  #noblank = noblank2[!(is.na(noblank2$KIPPSINDEX)),]
  #males = noblank[noblank$SEX=="MALE",]
  
  morph_0 = data
  #morph_0[is.na(morph_0)] = -1
  
  ## doesn't work with so many NA values -- replace with zero?
  ## PROBLEM: TAIL, LONG
  agg = aggregate(cbind(BILL.HEIGHT,BILL.LENGTH,BILL.WIDTH,TARSUS.LENGTH,WING.LENGTH.TO.PRIMARIES,
                        WING.LENGTH.TO.SECONDARIES,as.numeric(TAIL),LAT,as.numeric(LONG),KIPPSINDEX,STRUCTURE,BEAKBASEAREA,
                        BEAKVOL,BEAKLATERALSURFACE,BEAKTOTALSURFACE,BEAKAREAVOLUME,GENETIC.SIDE,GENETIC.SPLIT,SEQUENCED.) ~ CATALOG.NUMBER, data = morph_0, FUN=mean,na.action=na.pass)
  
  agg2 = aggregate(cbind(SPECIES,WHICH.SIDE.OF.CFB,COUNTRY,STATE,COUNTY,LOCALITY,SEX,AGE,CONDITION,BAD.MEASUREMENT,MEASURER,DATE,NOTES,
                         GENUS,SPP,SUBSPP,GEOREF.BY,LOCALITY2)~CATALOG.NUMBER,data=morph_0,
                   FUN=function(x) {paste(unique(x),collapse="; ")},na.action=na.pass)
  agg3 = aggregate(cbind(MEASUREMENT)~CATALOG.NUMBER,data=morph_0,FUN=base::max,na.action=na.pass)
  
  
  mergeAgg = merge(agg3,agg,by="CATALOG.NUMBER",all=T)
  mergeAgg = merge(mergeAgg,agg2,by="CATALOG.NUMBER",all=T)
  
  kipps = function(primary,secondary) {
    #kipps index calculated (primary_length-secondary_length/primary_length*100)
    return(((primary-secondary)/primary)*100)
  }
  beakbase = function(height,width){
    #beak base area calculated as area base cone (pi*beak_height*beak_width)
    return(pi*height*width)
  }
  beakvolume = function(height,width,length){
    #beak volume calculated as beak vol cone (1/3*beak_base*beak_length)
    base = beakbase(height,width)
    return((1/3)*base*length)
  }
  beakappxlatsurf = function(height,width,length){
    #beak lateral surface area calculated as appx lateral surf area cone ((beak_width+beak_height)/4)*beak_length*pi
    return(((width+height)/4)*length*pi)
  }
  beaktotsurf = function(height,width,length){
    #beak total surface is beak_lateral_surface + beak_base_area
    base = beakbase(height,width)
    lat = beakappxlatsurf(height,width,length)
    return(base+lat)
  }
  beaklatSAVR = function(height,width,length){
    #beak area volume is surface area / beak volume
    vol = beakvolume(height,width,length)
    lat = beakappxlatsurf(height,width,length)
    return(lat/vol)
  }
  beaktotSAVR = function(height,width,length){
    #beak area volume is surface area / beak volume
    vol = beakvolume(height,width,length)
    tot = beaktotsurf(height,width,length)
    return(tot/vol)
  }
  beakgeom = function(height,width,length){
    base = beakbase(height,width)
    vol = beakvolume(height,width,length)
    lat = beakappxlatsurf(height,width,length)
    tot = beaktotsurf(height,width,length)
    latvolratio = beaklatSAVR(height,width,length)
    totvolratio = beaktotSAVR(height,width,length)
    x = as.data.frame(cbind(base,vol,lat,tot,latvolratio,totvolratio))
    return(x)
  }
  
  head(mergeAgg)
  #mergeAgg[mergeAgg==-1] = NA
  
  kipp = kipps(primary=mergeAgg$WING.LENGTH.TO.PRIMARIES,secondary=mergeAgg$WING.LENGTH.TO.SECONDARIES)
  beaks = beakgeom(height=mergeAgg$BILL.HEIGHT,width=mergeAgg$BILL.WIDTH,length=mergeAgg$BILL.LENGTH)
  calculated = cbind(kipp,beaks,mergeAgg$CATALOG.NUMBER)
  names(calculated)
  names(calculated) = c("KIPPSINDEX_CALC","BEAKBASEAREA_CALC","BEAKVOL_CALC",
                        "BEAKLATERALSURFACE_CALC","BEAKTOTALSURFACE_CALC",
                        "BEAKAREAVOLUME_CALC","BEAKTOTAREAVOLUME_CALC","CATALOG.NUMBER")
  names(calculated)
  mergeAgg = merge(mergeAgg,calculated,by="CATALOG.NUMBER",all=T)
  
  path = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/measurements/aggregated datasets/"
  setwd(path)
  write.csv(mergeAgg,file=paste("KLP_Master_Spreadsheet_Morphology_AGGREGATED_",format(Sys.time(),"%d_%B_%Y"),".csv",sep=""),
            na="")
  
  morph = read.csv("KLP_Master_Spreadsheet_Morphology_AGGREGATED_19_July_2019.csv",
                     sep=",",header=T)
}; #dev.off()

  
## control for body size -- with and without log

print("BODY SIZE"); {
  
  names(morph)[10] = "TAIL"
  names(morph)[12] = "LONG"
  
  
  log = log(morph[,c(4:10,13,15:19,41:47)])
  names(log) = sapply(names(log),simplify=T,USE.NAMES = F, FUN=function(x) {paste("LOG",x,sep="_")})
  log = cbind(morph[,c(-4:-10,-13,-15:-19,-41:-47)],log)
  
  path = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/measurements/aggregated datasets/"
  setwd(path)
  write.csv(log,file=paste("KLP_Master_Spreadsheet_Morphology_AGGREGATED_LOG_",format(Sys.time(),"%d_%B_%Y"),".csv",sep=""),
            na="")
  
  model = lm(morph$BILL.HEIGHT~morph$TARSUS.LENGTH)
  plot(morph$TARSUS.LENGTH,morph$BILL.HEIGHT)
  abline(model)
  
  #model = lm(morph$BILL.HEIGHT~morph$TARSUS.LENGTH+morph$SPP)
  #plot(morph$TARSUS.LENGTH,morph$BILL.HEIGHT)
  #abline(model)
  
  #model = lm(log$LOG_BILL.HEIGHT~log$LOG_TARSUS.LENGTH)
  #plot(log$LOG_TARSUS.LENGTH,log$LOG_BILL.HEIGHT)
  #abline(model)
  
  #model = lm(log$LOG_BILL.HEIGHT~log$LOG_TARSUS.LENGTH+log$SPP)
  #plot(log$LOG_TARSUS.LENGTH,log$LOG_BILL.HEIGHT)
  #abline(model)
  
  ## get residuals after log
  
  model$residuals
  
  res_table = log[,2:27]
  rownames(res_table) = log$X
  for (i in 28:length(names(log))) {
    name = names(log)[i]
    print(name)
    
    temp=log[is.finite(rowSums(log[,c(i,31)])),]
    
    model = lm(temp[,i]~temp$LOG_TARSUS.LENGTH)
    res = residuals(model)
    res = unname(res)
    x = merge(res_table,res,by="row.names",all.x=TRUE)
    names(x)[i] = paste("RES_",name,sep="")
    rownames(x) = x$Row.names
    res_table = x[2:length(colnames(x))]
  }
  
  path = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/measurements/aggregated datasets/"
  setwd(path)
  write.csv(res_table,
            file=paste("KLP_Master_Spreadsheet_Morphology_AGGREGATED_LOG_residuals_",
                       format(Sys.time(),"%d_%B_%Y"),".csv",sep=""),
            na="")
  
  ## check to see if residuals are correlated
  #corrplot::corrplot(cor(res_table[,28:40],use="pairwise.complete.obs"),
  #                   method="ellipse",type="upper",diag=F)
  
  
  
  ## calculate residuals before log and compare beak/wing
  mor1 = morph[,c(4:10,13,15:19,41:47)]
  mor1 = cbind(morph[,c(-4:-10,-13,-15:-19,-41:-47)],mor1)
  
  
  nolog_res = mor1[,2:27]
  rownames(nolog_res) = mor1$X
  for (i in 28:length(names(mor1))) {
    name = names(mor1)[i]
    print(name)
    ## tarsus length is 31
    torun = mor1[,c(name,"TARSUS.LENGTH")]
    model = lm(torun[,name]~torun[,"TARSUS.LENGTH"],na.action=na.exclude)
    res = residuals(model)
    res = unname(res)
    x = merge(nolog_res,res,by="row.names",all.x=TRUE)
    names(x)[i] = paste("RES_",name,sep="")
    rownames(x) = x$Row.names
    nolog_res = x[2:length(colnames(x))]
  }
  
  path = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/measurements/aggregated datasets/"
  setwd(path)
  write.csv(nolog_res,
            file=paste("KLP_Master_Spreadsheet_Morphology_AGGREGATED_residuals_",
                       format(Sys.time(),"%d_%B_%Y"),".csv",sep=""),
            na="")
}; dev.off()
}
  
## MAKE DATA SUMMARY
print("SUMMARIZE"); {
  
  datafile="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/measurements/aggregated datasets/KLP_Master_Spreadsheet_Morphology_AGGREGATED_19_July_2019.csv"
  df=read.table(datafile,header=T,stringsAsFactors = T,sep="\t",blank.lines.skip = F)
  data=df[df$CATALOG.NUMBER!="",]
  data=data[complete.cases(data$CATALOG.NUMBER),]
  data=data[df$MEASUREMENT!="",]  
  data=data[complete.cases(data$MEASUREMENT),]
  ## get rid of ones without beak measurements
  data=data[complete.cases(data$BILL.HEIGHT),]
  ## get rid of the ones without tarsus measurements
  data=data[complete.cases(data$TARSUS.LENGTH),]
  ## there is now 1 without secondaries and 27 without tails
  ## this dataset already has females and juvs removed
  nobad = data[data$CONDITION=="",]
  
  path = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/"
  setwd(path)
  
  {
  #corrplot::corrplot(cor(nolog_res[,24:36],use="pairwise.complete.obs"),
  #                   method="ellipse",type="upper",diag=F)
  
  #nobad = nolog_res[nolog_res$CONDITION=="",]
  #noblank1 = nobad[nobad$RES_BEAKAREAVOLUME_CALC!="",]
  #noblank2 = noblank1[!(is.na(noblank1$RES_BEAKAREAVOLUME_CALC)),]
  #noblank = noblank2[!(is.na(noblank2$RES_KIPPSINDEX_CALC)),]
  #males = noblank[noblank$SEX=="MALE",]
  #males = nobad[nobad$SEX=="MALE" | nobad$SEX=="?MALE",]
  
  #males$SPP<-factor(males$SPP, levels=c("CARDINALIS","SINUATUS","FUSCA",
  #                                      "CURVIROSTRE","CRISSALE/DORSALE",
  #                                      "BILINEATA","BRUNNEICAPILLUS","NITENS",
  #                                      "FLAVICEPS","BELLII","MELANURA"))
  }
  
  males$SPP = factor(males$SPP, levels=sort(c("SINUATUS",#"CARDINALIS",
                                         "FUSCA",
                                         "BILINEATA","NITENS","FLAVICEPS",
                                         "MELANURA","CRISSALE",
                                         "CURVIROSTRE","BRUNNEICAPILLUS","BELLII")))
  
  names=sort(c("SIN",#"CAR",
          "FUS","BIL","NIT","FLA","MEL","CRI","CUR","BRU","BEL"))
  #col=c("white",#"cyan",
  #      "cyan","white","white","red","grey","white","cyan","red","cyan")
  col=rev(viridis::plasma(12))
  
  par(mfrow=c(1,1),lwd=1,font.main=4,mar=c(4,4,0,0))
  png("Analyses and Images/BOXPLOTS/aggregated_beakAreaVolume.png")
  boxplot(as.numeric(males$BEAKAREAVOLUME_CALC)~males$SPP,horizontal=F,
          names=names,las=2,
          xlab="Spp",ylab="Beak Surface Area to Volume Ratio",
          col=col)
  dev.off()
  
  ## DOESNT WORK
  {
  # permuteDeserts = function(SPECIES="MELANURA",N=10000,data=males,colNumber=30){
  #   
  #   ## if no measurements for a species will break
  #   
  #   d = data[data$SPP==SPECIES,]
  #   if (sum(is.na(d[,colNumber])) != 0) {
  #     w = which(is.na(d[,colNumber]))
  #     d = d[-w,]
  #   }
  #   
  #   if(nrow(d) != 0) {
  #     
  #     val = d[,colNumber]
  #     
  #     #val = val[!is.na(val)]
  #     name = colnames(d)[colNumber]
  #     des = d$WHICH.SIDE.OF.CFB
  #     
  #     if (length(unique(des)) == 1) {
  #       print("NOT DONE, TOO FEW TO COMPARE")
  #       return(NA)
  #     } else {
  #       
  #       ## will not work if only a single individual is present? 
  #       
  #       y = as.data.frame(cbind(as.numeric(val),as.character(des),z),stringsAsFactors = F)
  #       if (length(y$VALUE[y$ORIGINAL=="SONORAN"])  > 1 && length(y$VALUE[y$ORIGINAL=="CHIHUAHUAN"]) > 1)  {
  #         z = replicate(N,sample(des),simplify=T)
  #         
  #         colnames(y) = c("VALUE","ORIGINAL",paste("REP",1:N,sep="_"))
  #         y$VALUE = as.numeric(y$VALUE)
  #         
  #         print("org")
  #         
  #         
  #         
  #         org_ttest = t.test(y$VALUE[y$ORIGINAL=="SONORAN"],y$VALUE[y$ORIGINAL=="CHIHUAHUAN"])$p.value
  #         new_ttests = c()
  #         print("new")
  #         for (i in 1:N+2) {
  #           #print(colnames(y[i]))
  #           des1 = y[,c(1,i)]; colnames(des1) = c("VALUE","PERM")
  #           new_ttest = t.test(des1$VALUE[des1$PERM=="SONORAN"],des1$VALUE[des1$PERM=="CHIHUAHUAN"])$p.value
  #           new_ttests = c(new_ttests,new_ttest)
  #         }
  #         png(paste("Analyses and Images/permutedPValues_",name,"_",SPECIES,".png",sep=""))
  #         hist(new_ttests,main=paste(name,SPECIES,sep="\n"),xlab="p-value",ylab="freq",breaks=20)
  #         abline(v=org_ttest,col="red",lwd=2,lty=2)
  #         dev.off()
  #         print(summary(new_ttests))
  #         return(new_ttests)
  #       }
  #       else {
  #         print("ONLY SINGLETON, CANNOT DO")
  #         return(NA)
  #         
  #       }
  #     }
  #   }
  #   else {
  #     return(NA)
  #   }
  # }
  
  ## DOESNT WORK
  # for (colNumber in 27:length(names(males))) {
  #   for (SPECIES in unique(males$SPP)) {
  #     print(c(names(males)[colNumber],SPECIES))
  #     par(ask=F)
  #     
  #     ## think if no measurements for a species will break
  #     
  #     
  #     spp1 = males[males$SPP== SPECIES,] 
  #     spp_s = spp1[spp1$WHICH.SIDE.OF.CFB=="SONORAN",colNumber]
  #     spp_c = spp1[spp1$WHICH.SIDE.OF.CFB=="CHIHUAHUAN",colNumber]
  #     
  #     s = sum(males$SPP[males$WHICH.SIDE.OF.CFB=="SONORAN"] == SPECIES)
  #     c = sum(males$SPP[males$WHICH.SIDE.OF.CFB=="CHIHUAHUAN"] == SPECIES)
  #     print(s)
  #     print(c)
  #     if (c != 0 && s != 0) {
  #       permuteDeserts(SPECIES=SPECIES,N=100,data=males,colNumber=colNumber)
  #     }
  #   }
  # }
  }
    
  chi = males[males$WHICH.SIDE.OF.CFB=="CHIHUAHUAN",]
  son = males[males$WHICH.SIDE.OF.CFB=="SONORAN",]
  unk = males[males$WHICH.SIDE.OF.CFB=="UNCLEAR",]
  no_unk = males[males$WHICH.SIDE.OF.CFB!="UNCLEAR",]
  
  males$sppside=paste(substr(males$SPP,1,3),substr(males$WHICH.SIDE.OF.CFB,1,3))

  palette(c("cyan","green","grey"))
  par(mfrow=c(1,3),lwd=1,font.main=4,mar=c(4,4,0,0))
  png("Analyses and Images/BOXPLOTS/aggregated_beakLateralSurfArea_overlap_byDesert.png")
  boxplot(males$BEAKLATERALSURFACE_CALC~males$sppside,las=2,
          xlab="",col=c(1,2,3,1,2,1,2,1,2,1,2,1,2,3,1,2,3,1,2,1,2,1,2),
          ylab="Lateral Surface Area")
  abline(v=c(3.5,5.5,7.5,9.5,11.5,14.5,17.5,19.5,21.5),col="lightgrey",lty=2)
  dev.off()
  
  
  
  if(plotResid==T){
  par(mfrow=c(1,1),bg="black",col.axis="white",col.lab="white",
      col.main="white",col.sub="white",
      lwd=1,font.main=4,mar=c(2,2,2,1))
  png("Analyses and Images/BOXPLOTS/aggregated_beakAreaVolume_RESIDUALS.png")
  boxplot(as.numeric(males$RES_BEAKAREAVOLUME_CALC)~males$SPP,horizontal=F,
          names=names,las=2,
          xlab="Spp",ylab="Corr Lat Surf Area (mm2) / Corr Beak Vol (mm3)",
          col=col)
  dev.off()
  
  png("Analyses and Images/BOXPLOTS/aggregated_beakLateralSurfArea_RESIDUALS.png")
  boxplot(as.numeric(males$RES_BEAKLATERALSURFACE_CALC)~males$SPP,horizontal=F,
          names=names,las=2,xlab="Spp",ylab="Corrected Lat Surface Area (mm2)",
          col=col)
  dev.off()
  
  png("Analyses and Images/BOXPLOTS/aggregated_beakVolume_RESIDUALS.png")
  boxplot(as.numeric(males$RES_BEAKVOL_CALC)~males$SPP,horizontal=F,
          names=names,las=2,xlab="Spp",ylab="Corrected Log Beak Vol (mm3)",
          col=col)
  dev.off()
  
  summary(males$RES_BEAKAREAVOLUME_CALC)
  meanSAV = aggregate(males$RES_BEAKAREAVOLUME_CALC, by=list(males$SPP), FUN=mean)
  medSAV =  aggregate(males$RES_BEAKAREAVOLUME_CALC, by=list(males$SPP), FUN=median)
  
  par(mfrow=c(1,1),bg="black",col.axis="white",col.lab="white",
      col.main="white",col.sub="white",col="white",
      lwd=1,font.main=4,mar=c(5,4,2,1))
  png("Analyses and Images/BOXPLOTS/aggregated_beakLateralSurfArea_overlap_byDesert_RESIDUALS.png")
  boxplot(son$RES_BEAKLATERALSURFACE_CALC
          ~son$SPP,border=rgb(0,0.5,0,1),ylim=c(-100,200),add=F,col=rgb(0,1,0,1),
          las=2,boxwex=0.33,at=(seq(1,11)-2/3))
  boxplot(chi$RES_BEAKLATERALSURFACE_CALC
          ~chi$SPP,border=rgb(0,0.5,1,1),ylim=c(-100,200),add=T,col=rgb(0,1,1,1),
          las=2,boxwex=0.33,at=(seq(1,11)-1/3))
  boxplot(chi$RES_BEAKLATERALSURFACE_CALC
          ~chi$SPP,border=rgb(0,0.5,1,0),ylim=c(-100,200),add=T,col=rgb(0,1,1,0),
          las=2,boxwex=0.33,at=seq(0.5,10.5))
  dev.off()
  
  png("Analyses and Images/BOXPLOTS/aggregated_beakVolume_overlap_byDesert_RESIDUALS.png")
  boxplot(son$RES_BEAKVOL_CALC
          ~son$SPP,border=rgb(0,0.5,0,1),ylim=c(-800,2100),add=F,col=rgb(0,1,0,1),
          las=2,boxwex=0.33,names=rep("",11),at=(seq(1,11)-2/3))
  boxplot(chi$RES_BEAKVOL_CALC
          ~chi$SPP,border=rgb(0,0.5,1,1),ylim=c(-800,2100),add=T,col=rgb(0,1,1,1),
          las=2,boxwex=0.33,names=rep("",11),at=(seq(1,11)-1/3))
  boxplot(chi$RES_BEAKVOL_CALC
          ~chi$SPP,border=rgb(0,0.5,1,0),ylim=c(-800,2100),add=T,col=rgb(0,1,1,0),
          las=2,boxwex=0.33,at=seq(0.5,10.5))
  dev.off()
  
  png("Analyses and Images/BOXPLOTS/aggregated_beakAreaVolume_overlap_byDesert_RESIDUALS.png")
  boxplot(son$RES_BEAKAREAVOLUME_CALC
          ~son$SPP,border=rgb(0,0.5,0,1),ylim=c(-0.2,0.4),add=F,col=rgb(0,1,0,1),
          las=2,boxwex=0.33,names=rep("",11),at=(seq(1,11)-2/3))
  boxplot(chi$RES_BEAKAREAVOLUME_CALC
          ~chi$SPP,border=rgb(0,0.5,1,1),ylim=c(-0.2,0.4),add=T,col=rgb(0,1,1,1),
          las=2,boxwex=0.33,names=rep("",11),at=(seq(1,11)-1/3))
  boxplot(chi$RES_BEAKAREAVOLUME_CALC
          ~chi$SPP,border=rgb(0,0.5,1,0),ylim=c(-0.2,0.4),add=T,col=rgb(0,1,1,0),
          las=2,boxwex=0.33,at=seq(0.5,10.5))
  dev.off()
  
  png("Analyses and Images/BOXPLOTS/aggregated_kipps_overlap_byDesert_RESIDUALS.png")
  boxplot(son$RES_KIPPSINDEX_CALC
          ~son$SPP,border=rgb(0,0.5,0,1),ylim=c(-7,11),add=F,col=rgb(0,1,0,1),
          las=2,boxwex=0.33,names=rep("",11),at=(seq(1,11)-2/3))
  boxplot(chi$RES_KIPPSINDEX_CALC
          ~chi$SPP,border=rgb(0,0.5,1,1),ylim=c(-7,11),add=T,col=rgb(0,1,1,1),
          las=2,boxwex=0.33,names=rep("",11),at=(seq(1,11)-1/3))
  boxplot(chi$RES_KIPPSINDEX_CALC
          ~chi$SPP,border=rgb(0,0.5,1,0),ylim=c(-7,11),add=T,col=rgb(0,1,1,0),
          las=2,boxwex=0.33,at=seq(0.5,10.5))
  dev.off()
  
  png("Analyses and Images/BOXPLOTS/aggregated_kipps_overlap_byDesert_RESIDUALS.png")
  boxplot(son$RES_KIPPSINDEX_CALC
          ~son$SPP,border=rgb(0,0.5,0,1),ylim=c(-7,11),add=F,col=rgb(0,1,0,1),
          las=2,boxwex=0.33,names=rep("",11),at=(seq(1,11)-2/3))
  boxplot(chi$RES_KIPPSINDEX_CALC
          ~chi$SPP,border=rgb(0,0.5,1,1),ylim=c(-7,11),add=T,col=rgb(0,1,1,1),
          las=2,boxwex=0.33,names=rep("",11),at=(seq(1,11)-1/3))
  boxplot(chi$RES_KIPPSINDEX_CALC
          ~chi$SPP,border=rgb(0,0.5,1,0),ylim=c(-7,11),add=T,col=rgb(0,1,1,0),
          las=2,boxwex=0.33,at=seq(0.5,10.5))
  dev.off()
  
  for (i in c(27:length(names(son)))) {
    n = names(son)[i]
    min = floor(min(son[,i],chi[,i],na.rm=T)*4)/4
    max = ceiling(max(son[,i],chi[,i],na.rm=T)*4)/4
    png(paste("Analyses and Images/BOXPLOTS/aggregated_",n,"_byDesert_RESIDUALS.png",sep=""),
        pointsize=18)
    par(mfrow=c(1,1),bg="black",col.axis="white",col.lab="white",
        col.main="white",col.sub="white",col="white",
        lwd=1,font.main=4,mar=c(5,4,2,1))
    boxplot(son[,i]~son$SPP,border=rgb(0,0.5,0,1),ylim=c(min,max),add=F,col=rgb(0,1,0,1),
            las=2,boxwex=0.33,names=rep("",11),at=(seq(1,11)-2/3),main=n)
    boxplot(chi[,i]~chi$SPP,border=rgb(0,0.5,1,1),ylim=c(min,max),add=T,col=rgb(0,1,1,1),
            las=2,boxwex=0.33,names=rep("",11),at=(seq(1,11)-1/3))
    boxplot(chi[,i]~chi$SPP,border=rgb(0,0.5,1,0),ylim=c(min,max),add=T,col=rgb(0,1,1,0),
            las=2,boxwex=0.33,at=seq(0.5,10.5))
    dev.off()
  }
  }
  
  if(generateData==T){
  for (spp in unique(males$SPP)) {
    for (des in unique(males$WHICH.SIDE.OF.CFB)) {
      print(paste(spp,des))
      tab = males[males$SPP==spp & males$WHICH.SIDE.OF.CFB==des,c(27:length(names(males)))]
      m = as.data.frame(sapply(tab,mean))
      #s = sapply(tab,sd)
      #names(spp) = "SPP"
      #names(des) = "WHICH.SIDE.OF.CFB"
      sppdes = paste(substring(spp,1,3),substring(des,1,3),sep="_")
      names(m) = sppdes
      tm = t(m)
      write(paste("SPP_DES",names(tab),collapse=",",sep=","),"measurements/aggregated datasets/KLP_Master_Spreadsheet_Morphology_AGGREGATED_RESIDUALS_MEAN.csv",
            append=T,sep=",")
      write.table(tm,"measurements/aggregated datasets/KLP_Master_Spreadsheet_Morphology_AGGREGATED_RESIDUALS_MEAN.csv",
                  append=T,sep=",",col.names=F)
    }
  }
  
  for (spp in unique(males$SPP)) {
    print(spp)
    write(spp,"Analyses and Images/BOXPLOTS/output_RES.txt",append=T,sep="\t",ncolumns=3)
    tabS = males[males$SPP==spp & males$WHICH.SIDE.OF.CFB=="SONORAN",c(27:length(names(males)))]
    tabC = males[males$SPP==spp & males$WHICH.SIDE.OF.CFB=="CHIHUAHUAN",c(27:length(names(males)))]
    for (colnum in 1:length(colnames(tabS))) {
      write(colnames(tabS)[colnum],"Analyses and Images/BOXPLOTS/output_RES.txt",append=T,sep="\t",ncolumns=3)
      print(colnum)
      if ((nrow(tabC) != 0) && (nrow(tabS) != 0) && sum(!is.na(tabC[,colnum])) > 1 && sum(!is.na(tabS[,colnum])) > 1) {
        test = t.test(tabS[,colnum],tabC[,colnum])
        write(test$p.value,"Analyses and Images/BOXPLOTS/output_RES.txt",append=T,sep="\n",ncolumns=3)
      } else {
        write("#####","Analyses and Images/BOXPLOTS/output_RES.txt",append=T,sep="\n",ncolumns=3)
      }
    }
  }
  }
  
  ## need to exit script and manually reformat this? 
  
  #testtable = read.table("sonVsChit-test_RES.txt",header=T,sep="\t")
  
}; dev.off()

## import data
print("IMPORT"); {
  ## for background
  #Env = raster::stack(list.files(
  #  path='/Users/kprovost/Dropbox (AMNH)/Classes/Spatial Bioinformatics/spatial_bioinformatics-master/ENM/wc2-5/',
  #  pattern="\\.bil$",
  #  full.names=T))
  Env=raster("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/ENMS_multilayer.tif")
  #ext = raster::extent(c(-125,-60,10,50)) ## make sure this will play nice with your points
  ext = raster::extent(c(-118,-96,21,37)) ## make sure this will play nice with your points
  Env = raster::crop(Env, ext)
  bg = Env[[1]] ## just for plotting 
  
  path = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/measurements/"
  # fileprefix1 = "KLP_Master_Spreadsheet_Morphology_"
  # fileprefix2 = "aggregated datasets/KLP_Master_Spreadsheet_Morphology_"
  # date1 = "19_July_2019"
  # date2 = "19_July_2019"
  # 
  # ## agg refers to aggregated (as in, individual readings are merged)
  # ## res refers to residuals in measurements after accounting for tarsus length
  # notagg = read.csv(paste(path,fileprefix1,date2,".csv",sep=""))
  # agg = read.csv(paste(path,fileprefix2,"AGGREGATED_",date2,".csv",sep=""))
  # agg_res = read.csv(paste(path,fileprefix2,"AGGREGATED_residuals_",date2,".csv",sep=""))
  # agg_log_res = read.csv(paste(path,fileprefix2,"AGGREGATED_LOG_residuals_",date2,".csv",sep=""))
  # names(agg)[c(10,12)] = c("TAIL","LONG")
  # 
  # notagg_good = notagg[notagg$CONDITION=="",]
  # notagg_good = notagg_good[!(notagg_good$AGE %in% c("JUV","?JUV")),]
  # notagg_good_males = notagg_good[notagg_good$SEX %in% c("MALE","?MALE"),]
  # 
  # 
  # agg_good = agg[agg$CONDITION=="",]
  # agg_good = agg_good[!(agg_good$AGE %in% c("JUV","?JUV")),]
  # agg_good_males = agg_good[agg_good$SEX %in% c("MALE","?MALE"),]
  # 
  # agg_res_good = agg_res[agg_res$CONDITION=="",]
  # agg_res_good = agg_res_good[agg_res_good$AGE!="JUV" && agg_res_good$AGE!="?JUV",]
  # agg_res_good_males = agg_res_good[agg_res_good$SEX=="MALE",]
  # 
  # agg_fus = agg_good_males[agg_good_males$SPP=="FUSCA",]
  # agg_res_fus = agg_res_good_males[agg_res_good_males$SPP=="FUSCA",]
  
  datafile="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/measurements/aggregated datasets/KLP_Master_Spreadsheet_Morphology_AGGREGATED_19_July_2019.csv"
  df=read.table(datafile,header=T,stringsAsFactors = T,sep="\t",blank.lines.skip = F)
  data=df[df$CATALOG.NUMBER!="",]
  data=data[complete.cases(data$CATALOG.NUMBER),]
  data=data[df$MEASUREMENT!="",]  
  data=data[complete.cases(data$MEASUREMENT),]
  ## get rid of ones without beak measurements
  data=data[complete.cases(data$BILL.HEIGHT),]
  ## get rid of the ones without tarsus measurements
  data=data[complete.cases(data$TARSUS.LENGTH),]
  ## there is now 1 without secondaries and 27 without tails
  ## this dataset already has females and juvs removed
  nobad = data[data$CONDITION=="",]
  males=nobad
  
  pal = colorRampPalette(RColorBrewer::brewer.pal(9,"YlOrBr"))
  col <- pal(10)[as.numeric(cut(agg_fus$TARSUS.LENGTH,breaks = 10))]
  col[is.na(col)] = "black"
  #raster::plot(bg, col="grey",colNA="darkgrey",legend=F,main="TARSUS")
  #points(agg$LONG,agg$LAT)
  #points(agg_fus$LONG,agg_fus$LAT,col=col,pch=2,bg=col)
  
  #col <- pal(10)[as.numeric(cut(agg_fus$BEAKTOTAREAVOLUME_CALC,breaks = 10))]
  #col[is.na(col)] = "black"
  #raster::plot(bg, col="grey",colNA="darkgrey",legend=F,main="Beak Surface Area to Volume")
  #points(agg$LONG,agg$LAT)
  #points(agg_fus$LONG,agg_fus$LAT,col=col,pch=2,bg=col)
  
  #col <- pal(10)[as.numeric(cut(agg_res_fus$RES_BEAKTOTAREAVOLUME_CALC,breaks = 10))]
  #col[is.na(col)] = "black"
  #raster::plot(bg, col="grey",colNA="darkgrey",legend=F,main="Beak Surface Area to Volume (residuals)")
  #points(agg$LONG,agg$LAT)
  #points(agg_res_fus$LONG,agg_res_fus$LAT,col=col,pch=2,bg=col)
}; #dev.off()

## export imagery
print("EXPORT"); {
  agg=males
  pdf("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/SAMPLING/Tarsus_and_beakSAVR.pdf")
  par(mfrow=c(4,3),mar=c(1,1,1,1),ask=F)
  for(species in unique(agg$SPP)){
    print(species)
    
    if (species != "") {
      
      agg_spp = agg[agg$SPP==species,]
      agg_res_spp = agg_res_good_males[agg_res_good_males$SPP==species,]
      
      par(mfrow=c(1,2))
      
      pal = colorRampPalette(RColorBrewer::brewer.pal(9,"Spectral"))
      col <- pal(10)[as.numeric(cut(agg_spp$TARSUS.LENGTH,breaks = 10))]
      col[is.na(col)] = "black"
      raster::plot(bg, col="grey",colNA="darkgrey",legend=F,main=paste(species,"Tarsus length"))
      #points(agg$LONG,agg$LAT)
      points(agg_spp$LONG,agg_spp$LAT,col=col,pch=2,bg=col)
      
      col <- pal(10)[as.numeric(cut(agg_res_spp$RES_BEAKTOTAREAVOLUME_CALC,breaks = 10))]
      col[is.na(col)] = "black"
      raster::plot(bg, col="grey",colNA="darkgrey",legend=F,main=paste(species,"Beak Surface Area to Volume (residuals)"))
      #points(agg$LONG,agg$LAT)
      points(agg_res_spp$LONG,agg_res_spp$LAT,col=col,pch=2,bg=col)
    }
  }
  dev.off()
  
  pdf("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/SAMPLING/Which_side_of_CFB_by_species.pdf")
  par(mfrow=c(3,4),mar=c(1,1,3,1),ask=F)
  for(species in unique(agg$SPP)){
    print(species)
    if (species != "") {
    agg_spp = agg[agg$SPP==species,]
    agg_spp$GENETIC.SIDE[agg_spp$GENETIC.SIDE==-1] = 2
    #agg_res_spp = agg_res_good_males[agg_res_good_males$SPP==species,]
    #agg_res_spp$GENETIC.SIDE[agg_res_spp$GENETIC.SIDE==-1] = 2
    
    
    if ((unique(agg_spp$GENETIC.SPLIT)!=1)) {
      print("not split or unknown")
      tocheck = as.numeric(agg_spp$WHICH.SIDE.OF.CFB)
    } else {
      print("split")
      tocheck = as.numeric(agg_spp$GENETIC.SIDE)
    }
    
    print(tocheck)
    
    
    #par(mfrow=c(1,2))
    
    ## add an
    
    pal = colorRampPalette(RColorBrewer::brewer.pal(9,"Spectral"))
    col <- pal(3)[tocheck]
    col[is.na(col)] = "black"
    raster::plot(bg, col="grey",colNA="darkgrey",legend=F,main=paste(species))
    #points(agg$LONG,agg$LAT)
    points(agg_spp$LONG,agg_spp$LAT,col=col,pch=2,bg=col)
    }
  }
  dev.off()
}; #dev.off()

## quick t-tests
if(generateData==T){
print("T TESTS"); {
  ## check if the value of corrected surface area to volume is different between deserts
  ## except crissals 
  
  towrite = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/MODELS/t.tests.csv"
  
  for(species in unique(agg_res_good_males$SPP)){
    
    print(species)
    write(species,towrite,append=T)
    agg_res_spp = agg_res_good_males[agg_res_good_males$SPP==species,]
    meanVal = mean(agg_res_spp$RES_BEAKLATERALSURFACE_CALC,na.rm=T)
    print(paste("Species mean value RES_BEAKLATERALSURFACE_CALC:",meanVal))
    write(paste("Species mean value RES_BEAKLATERALSURFACE_CALC:",meanVal),towrite,append=T)
    
    
    #if (species != "CRISSALE/DORSALE"){
    son = agg_res_spp$RES_BEAKLATERALSURFACE_CALC[agg_res_spp$WHICH.SIDE.OF.CFB=="SONORAN"]
    chi = agg_res_spp$RES_BEAKLATERALSURFACE_CALC[agg_res_spp$WHICH.SIDE.OF.CFB=="CHIHUAHUAN"]
    print(paste("N Son:",length(son),"N Chi:",length(chi)))
    write(paste("N Son:",length(son),"N Chi:",length(chi)),towrite,append=T)
    print(paste("Pvalue: ",t.test(son,chi)$p.value))
    write(paste("Pvalue: ",t.test(son,chi)$p.value),towrite,append=T)
    
    
    if (t.test(son,chi)$p.value < 0.05) {
      print("SIG")
      write("SIG",towrite,append=T)
      
    }
    
    #}
  }
  
  
  ## for fun compare within genera
  t.test(agg_good_males$TARSUS.LENGTH[agg_res_good_males$SPP=="CRISSALE"],
         agg_good_males$TARSUS.LENGTH[agg_res_good_males$SPP=="CURVIROSTRE"])
  #t.test(agg_good_males$TARSUS.LENGTH[agg_res_good_males$SPP=="CARDINALIS"],
  #       agg_good_males$TARSUS.LENGTH[agg_res_good_males$SPP=="SINUATUS"])
  
  
  names(agg_res_good_males[28:47])
  for(i in c(28:33,35:47)){
    
    
    print("xxxxxxxxxxxxxxxxxxxxxxxxxx")
    write("xxxxxxxxxxxxxxxxxxxxxxxxxx",towrite,append=T)
    print(names(agg_res_good_males[i]))
    write(names(agg_res_good_males[i]),towrite,append=T)
    print("TOXOSTOMA")
    write("TOXOSTOMA",towrite,append=T)
    print(t.test(agg_res_good_males[agg_res_good_males$SPP=="CRISSALE",i],
                 agg_res_good_males[agg_res_good_males$SPP=="CURVIROSTRE",i])$p.value)
    write(t.test(agg_res_good_males[agg_res_good_males$SPP=="CRISSALE",i],
                 agg_res_good_males[agg_res_good_males$SPP=="CURVIROSTRE",i])$p.value,towrite,append=T)
    #print("CARDINALIS")
    #print(t.test(agg_res_good_males[agg_res_good_males$SPP=="CARDINALIS",i],
    #             agg_res_good_males[agg_res_good_males$SPP=="SINUATUS",i])$p.value)
  }
}; dev.off()


## run a simple GLMM -- is longitude predictable?  

print("GLMM"); {
  
  torun = agg[,c(4:10,14,24,37,41:47)]
  #agg_good_males = agg_good_males[!is.na(agg_good_males$BEAKTOTALSURFACE_CALC),]
  model = glm(STRUCTURE ~ BILL.HEIGHT + BILL.LENGTH + BILL.WIDTH + WING.LENGTH.TO.PRIMARIES +
                WING.LENGTH.TO.SECONDARIES + KIPPSINDEX_CALC + BEAKBASEAREA_CALC + BEAKVOL_CALC +
                BEAKLATERALSURFACE_CALC + BEAKAREAVOLUME_CALC 
              #+ SPP
              ,data=torun,family = "binomial")
  summary(model) ## long vs rest, w/ spp: beak vol sig, beak lat surf near sig, beaktotsurf NA, bil sig
  car::Anova(model)
  ## significant ones long vs rest w/ spp: beak vol sig, spp near sig 
  
  library(AICcmodavg)
  
  AICc(model)
  
  model1 = glm(STRUCTURE ~ BILL.HEIGHT,data=torun,family = "binomial")
  model2 = glm(STRUCTURE ~ BILL.LENGTH,data=torun,family = "binomial")
  model3 = glm(STRUCTURE ~ BILL.WIDTH,data=torun,family = "binomial")
  model4 = glm(STRUCTURE ~ TARSUS.LENGTH,data=torun,family = "binomial")
  model5 = glm(STRUCTURE ~ WING.LENGTH.TO.PRIMARIES,data=torun,family = "binomial")
  model6 = glm(STRUCTURE ~ WING.LENGTH.TO.SECONDARIES,data=torun,family = "binomial")
  model7 = glm(STRUCTURE ~ TAIL,data=torun,family = "binomial")
  model8 = glm(STRUCTURE ~ KIPPSINDEX_CALC,data=torun,family = "binomial")
  model9 = glm(STRUCTURE ~ BEAKBASEAREA_CALC,data=torun,family = "binomial")
  model10 = glm(STRUCTURE ~ BEAKVOL_CALC,data=torun,family = "binomial")
  model11 = glm(STRUCTURE ~ BEAKLATERALSURFACE_CALC,data=torun,family = "binomial")
  model12 = glm(STRUCTURE ~ BEAKTOTALSURFACE_CALC,data=torun,family = "binomial")
  model13 = glm(STRUCTURE ~ BEAKAREAVOLUME_CALC,data=torun,family = "binomial")
  model14 = glm(STRUCTURE ~ BEAKTOTAREAVOLUME_CALC,data=torun,family = "binomial")
  
  modellist = list(model1,model2,model3,model4,
                   model5,model6,model7,model8,
                   model9,model10,model11,model12,
                   model13,model14)
  for (i in 1:14) {
    model = modellist[[i]]
    print(paste(AICcmodavg::AICc(model),i,sep=" "))
  }
  
  ## now model 2 is the best -- TAIL plus BILL LENGTH 
  ## now model 6 -- BILL.LENGTH + TAIL + WING.LENGTH.TO.SECONDARIES
  ## stabilized on 6 -- add species 
  
  
  model6 = glm(STRUCTURE ~ BILL.LENGTH + TAIL + WING.LENGTH.TO.SECONDARIES,data=torun,family = "binomial")
  model6a = glm(STRUCTURE ~ BILL.LENGTH + WING.LENGTH.TO.SECONDARIES,data=torun,family = "binomial")
  
  ## model 7 is the best so far -- TAIL only 
  variable = torun$TAIL
  names(variable) = "TAIL"
  plot(torun$TAIL,torun$GENETIC.SPLIT,xlab="tailsize",ylab="structure",pch=22,
       bg=torun$GENETIC.SPLIT) # plot with body size on x-axis and survival (0 or 1) on y-axis
  points(torun$TAIL, predict(model6a, newdata=torun, type="response"),pch=22,
         bg=torun$GENETIC.SPLIT)
  
  model_nospp = glm(LONG ~ BILL.HEIGHT + BILL.LENGTH + BILL.WIDTH + TARSUS.LENGTH + WING.LENGTH.TO.PRIMARIES +
                      WING.LENGTH.TO.SECONDARIES + TAIL + KIPPSINDEX_CALC + BEAKBASEAREA_CALC + BEAKVOL_CALC +
                      BEAKLATERALSURFACE_CALC + BEAKTOTALSURFACE_CALC + BEAKAREAVOLUME_CALC + BEAKTOTAREAVOLUME_CALC 
                    ,data=agg_good_males)
  summary(model_nospp) ## bill height/length/width sig, beak base, beak vol, beat kat surf sig, beak tot area NA, 
  car::Anova(model_nospp) ## bil he/le/wi sig, beak vol sig
  
  
  
  
}; dev.off()

## do univariate models with species 
print("UNIVAR"); {
  model_onlyspp = glm(LONG ~ SPP,data=agg_good_males)
  model_BH = glm(LONG ~ BILL.HEIGHT+ SPP,data=agg_good_males)
  model_BL = glm(LONG ~ BILL.LENGTH+ SPP,data=agg_good_males)
  model_BW = glm(LONG ~ BILL.WIDTH+ SPP,data=agg_good_males)
  model_TL = glm(LONG ~ TARSUS.LENGTH+ SPP,data=agg_good_males)
  model_WP = glm(LONG ~ WING.LENGTH.TO.PRIMARIES+ SPP,data=agg_good_males)
  model_WS = glm(LONG ~ WING.LENGTH.TO.SECONDARIES+ SPP,data=agg_good_males)
  model_TA = glm(LONG ~ TAIL+ SPP,data=agg_good_males)
  model_KI = glm(LONG ~ KIPPSINDEX_CALC+ SPP,data=agg_good_males)
  model_BBA = glm(LONG ~ BEAKBASEAREA_CALC+ SPP,data=agg_good_males)
  model_BV = glm(LONG ~ BEAKVOL_CALC+ SPP,data=agg_good_males)
  model_BLS = glm(LONG ~ BEAKLATERALSURFACE_CALC+ SPP,data=agg_good_males)
  model_BTS = glm(LONG ~ BEAKTOTALSURFACE_CALC+ SPP,data=agg_good_males)
  model_BAV = glm(LONG ~ BEAKAREAVOLUME_CALC+ SPP,data=agg_good_males)
  model_BTAV = glm(LONG ~ BEAKTOTAREAVOLUME_CALC + SPP,data=agg_good_males)
  
  summary(model_onlyspp) ## fus sig, cris near sig
  car::Anova(model_onlyspp) ## spp sig
  
  summary(model_BH) ## cris near sig 
  car::Anova(model_BH) ## spp sig
  
  summary(model_BL) ## fus near sig 
  car::Anova(model_BL) ## spp sig 
  
  summary(model_BW) ## width sig, card cris curv fus nit sin sig 
  car::Anova(model_BW) ## width and spp sig 
  
  summary(model_TL) ## tarsis sig, cur sin sig 
  car::Anova(model_TL) ## tarsus spp sig 
  
  model_TL_rev = glm(TARSUS.LENGTH ~ LONG + I(LONG^2) + SPP,data=agg_good_males)
  summary(model_TL_rev) ## no square basically everything is signif except bil, mel, nit
  car::Anova(model_TL_rev) ## long and spp sig 
  ## if you add the square nit is sig, mel is close, bil still not 
}; dev.off()

## quadratic plots
print("QUADRATIC"); {
  generateQuadraticPlot = function(spp="SINUATUS",xvariable="LONG",
                                   yvariable="TARSUS.LENGTH",maxquadratic=2,data=agg_good_males,
                                   colornumber=11) {
    
    smalldata=data[,c(xvariable,yvariable,"SPP")]
    sppdata = data[data$SPP==spp,c(xvariable,yvariable)]
    names(sppdata) = c("X","Y")
    names(smalldata) = c("X","Y","SPP")
    
    model = glm(Y ~ poly(X,maxquadratic) + SPP,data=smalldata)
    #print(summary(model))
    sppmodel = glm(Y ~ poly(X,maxquadratic),data=sppdata)
    #print(summary(sppmodel))
    pvalues = as.numeric(coef(summary(sppmodel))[,4][2:(maxquadratic+1)])
    sig = pvalues <= 0.05
    
    plot(smalldata$X,smalldata$Y,col="white",
         pch=as.numeric(data$SPP))
    points(sppdata$X,sppdata$Y,col=colornumber,pch=colornumber)
    topredict = as.data.frame(seq(min(smalldata$X),max(smalldata$X),0.01))
    names(topredict) = "X"
    fitted = predict.lm(object=sppmodel,newdata=topredict)
    points(topredict$X,fitted,type="l")
    legend("topleft",legend=c(spp,formatC(pvalues, format = "e", digits = 2)),bty="n",
           text.col = c("black",as.numeric(sig)+1))
  }
  
  pdfstring = paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/MODELS/Results_of_models.pdf",sep="")
  pdf(pdfstring,height=3,width=8)
  for (variablenum in c(4:9,13,15:19,41:47)) {
    ## tail is skipped
    var = names(agg_good_males)[variablenum]
    print(var)
    
    par(mfrow=c(2,6),mar=c(1,1,1,1),ask=F)
    pal = c("red","orange","green","blue","purple","cyan","black","grey","magenta","pink","darkgreen")
    palette(pal)
    options(digits = 3, scipen = -2)
    plot(0,type='n',axes=FALSE,ann=FALSE)
    legend("center",legend=c("Longitude","vs",var))
    for (i in 1:length(unique(agg_good_males$SPP))) {
      spp = (as.character(sort(unique(agg_good_males$SPP))[i]))
      
      generateQuadraticPlot(spp=spp,xvariable="LONG",yvariable=var,
                            maxquadratic = 2,data=agg_good_males,colornumber=i)
      
    }
    
  }
  dev.off()
  
  pdfstring = paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/MODELS/Results_of_models_residuals.pdf",sep="")
  pdf(pdfstring,height=3,width=8)
  for (variablenum in c(28:33,35:47)) {
    ## tail is skipped
    var = names(agg_res_good_males)[variablenum]
    print(var)
    
    par(mfrow=c(2,6),mar=c(1,1,1,1),ask=F)
    pal = c("red","orange","green","blue","purple","cyan","black","grey","magenta","pink","darkgreen")
    palette(pal)
    options(digits = 3, scipen = -2)
    plot(0,type='n',axes=FALSE,ann=FALSE)
    legend("center",legend=c("Longitude","vs",var))
    for (i in 1:length(unique(agg_res_good_males$SPP))) {
      spp = (as.character(sort(unique(agg_res_good_males$SPP))[i]))
      
      generateQuadraticPlot(spp=spp,xvariable="LONG",yvariable=var,
                            maxquadratic = 2,data=agg_res_good_males,colornumber=i)
      
    }
    
  }
  dev.off()
  
  # pdfstring = paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Results_of_models_notaggregated.pdf",sep="")
  # pdf(pdfstring,height=3,width=8)
  # for (variablenum in c(11:16,30,32:36)) {
  #   ## tail is skipped
  #   var = names(notagg_good_males)[variablenum]
  #   print(var)
  #   
  #   par(mfrow=c(2,6),mar=c(1,1,1,1),ask=F)
  #   pal = c("red","orange","green","blue","purple","cyan","black","grey","magenta","pink","darkgreen")
  #   palette(pal)
  #   options(digits = 3, scipen = -2)
  #   plot(0,type='n',axes=FALSE,ann=FALSE)
  #   legend("center",legend=c("Longitude","vs",var))
  #   for (i in 1:length(unique(notagg_good_males$SPP))) {
  #     spp = (as.character(sort(unique(notagg_good_males$SPP))[i]))
  #     
  #     generateQuadraticPlot(spp=spp,xvariable="LONG",yvariable=var,
  #                           maxquadratic = 2,data=notagg_good_males,colornumber=i)
  #     
  #   }
  #   
  # }
  # dev.off()
}; dev.off()


## get error between measurements

if(generateData==T)
print("ERROR"); {
  names(notagg)
  outfile = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/ERROR/Measurements_standard_error_within_individuals.temp"
  table = data.frame(CATALOG.NUMBER=character(),
                     COLUMN=character(), 
                     AVERAGE=numeric(), 
                     STDEV=numeric(),
                     MEASUREMENTS=numeric(),
                     MAX=numeric(),
                     MIN=numeric(),
                     DIFMAXMIN=numeric(),
                     stringsAsFactors=F) 
  for(ind in unique(notagg$CATALOG.NUMBER)){
    x = notagg[notagg$CATALOG.NUMBER==ind,]
    #print(x)
    ## cols are 11 to 17
    if(nrow(x)!=1){
      for(colnum in 11:17){
        if(!(is.na(mean(x[,colnum])))){
          row = cbind(as.character(ind),
                      names(x[colnum]),
                      as.numeric(mean(x[,colnum],na.rm=T)),
                      as.numeric(sd(x[,colnum],na.rm=T)),
                      as.numeric(nrow(x)),
                      max(x[,colnum],na.rm=T),
                      min(x[,colnum],na.rm=T),
                      abs(max(x[,colnum],na.rm=T) - min(x[,colnum],na.rm=T)))
          colnames(row) = colnames(table)
          write.table(row,outfile,append=T)
          table = rbind(table,row)
        }
      }
      
      #print(paste(as.character(ind),names(x[colnum]),"avg",mean(x[,colnum]),"stdev",sd(x[,colnum]),"MEASUREMENTS",nrow(x),sep=" "))
    }
  }
  
  
  ## errors 
  
  table2 = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/ERROR/Measurements_standard_error_within_individuals.txt",
                    sep="\t")
  
  #table$STDEV = as.numeric(table$STDEV)
  #table$AVERAGE = as.numeric(table$AVERAGE)
  #table$DIFMAXMIN = as.numeric(table$DIFMAXMIN)
  #table$MAX=as.numeric(table$MAX)
  #table$MIN=as.numeric(table$MIN)
  hist(table2$DIFMAXMIN)
  zzz = (table2[which(table2$DIFMAXMIN>=1),1:2])
  zzz[zzz$COLUMN=="TAIL",]
  
  boxplot(table2$DIFMAXMIN~table2$COLUMN)
  
  
  hist(table2$STDEV,xlim=c(0,2))
  yyy = table2[which(table2$STDEV>=0.5),1:2]
}; dev.off()
  
#### do a clustering analysis 
print("CLUSTER"); {
  #tocluster = agg_good_males[,c(11:16,30,32:36)] ## not tail, too many missing
  tocluster = agg[,c(4:12,24,26,37,41:47)]
  tocluster = tocluster[complete.cases(tocluster),]
  #clustered = kmeans(tocluster[,c(11:16,30,32:36)],centers=11)
  clustered = kmeans(tocluster[,c(1:7,13:19)],centers=10)
  #fpc::plotcluster(tocluster, clustered$cluster)
  
  postcluster = cbind(tocluster,clustered$cluster )
  names(postcluster) = c(names(tocluster),"assignment")
  boxplot(postcluster$assignment~postcluster$SPP)
  
  ## try within deserts 
  
  pdf("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/KMEANS/cluster_assignment_k3.pdf")
  palette("default")
  for (spp in sort(unique(tocluster$SPP))){
    print(spp)
    temp = tocluster[tocluster$SPP==spp,]
    temp_clus = kmeans(temp[,c(1:7,13:19)],centers=3)
    posttemp = cbind(temp,temp_clus$cluster)
    names(posttemp) = c(names(temp),"assignment")
    par(mfrow=c(2,2))
    #print(table(posttemp$STATE,posttemp$assignment))
    plot(posttemp$BILL.HEIGHT,posttemp$BILL.WIDTH,col=posttemp$assignment,main=spp,
         pch=posttemp$assignment)
    m = table(posttemp$assignment,as.character(posttemp$STATE))
    tab = prop.table(m, margin=2)
    barplot(tab,las=2,col=c("black","red","green"))
    
    
    m2 = table(posttemp$assignment,as.character(posttemp$WHICH.SIDE.OF.CFB))
    tab2 = prop.table(m2, margin=2)
    barplot(tab2,las=2,col=c("black","red","green"))
    
    
    raster::plot(bg, col="grey",colNA="darkgrey",legend=F,main=paste(spp))
    points(as.numeric(as.character(posttemp$LONG)),as.numeric(as.character(posttemp$LAT)),col=posttemp$assignment,main=spp,
           pch=posttemp$assignment)
    #readline(prompt="Press any key to continue.")
  }
  dev.off()
  
}; dev.off()

#### clustering with silouette
print("CLUSTER2"); {
  
  subsetted = agg[complete.cases(agg),]
  subsetted = subsetted[,c(4:10,41:47,37)]
  
  avg_sil <- function(k) {
    km.res <- kmeans(df, centers = k, nstart = 25)
    ss <- silhouette(km.res$cluster, dist(df))
    mean(ss[, 3])
  }
  
  par(mfrow=c(2,5))
  for (spp in unique(subsetted$SPP)) {
    print(spp)
    df = subsetted[subsetted$SPP==spp,-15]
    
    # k.values <- 2:min(nrow(df)-1,10)
    # 
    # # extract avg silhouette for 2-15 clusters
    # avg_sil_values <- map_dbl(k.values, avg_sil)
    # 
    # plot(k.values, avg_sil_values,
    #      type = "b", pch = 19, frame = FALSE,
    #      xlab = "Number of clusters K",
    #      ylab = "Average Silhouettes",
    #      main=spp)
    # 
    #readline(prompt="Press [enter] to continue")
    
    X = fviz_nbclust(df, kmeans, method = "silhouette",
                     k.max=(min(nrow(df)-1,10)))
    #text(x=4,y=0.2,labels=spp)
    
    #readline(prompt=spp)
    
  }
  
  bestk=c(2,2,2,3,2,3,2,2,3,2)
}; dev.off()
}

#### run dabest analysis 
if(generateData==T){
print("DABEST"); {
  
  ## CUSTOM PLOT DABEST FOR WHITE LINES
  {
    get_tick_labels <- function(plot.obj, axes) {
      # This works for ggplot2 v3.0.0; not fully tested with other versions yet.
      plot.obj.build <- ggplot2::ggplot_build(plot.obj)
      
      if (axes == "x") {
        ticks <- plot.obj.build$layout$panel_params[[1]]$x.labels
      } else if (axes == "y") {
        ticks <- plot.obj.build$layout$panel_params[[1]]$y.labels
      } else {
        stop('`axes` must be either "x" or "y"')
      }
      
      return(ticks)
    }
    
    
    
    get_tick_positions <- function(plot.obj, axes) {
      # This works for ggplot2 v3.0.0; not fully tested with other versions yet.
      plot.obj.build <- ggplot2::ggplot_build(plot.obj)
      
      if (axes == "x") {
        pos <- plot.obj.build$layout$panel_params[[1]]$y.minor_source
      } else if (axes == "y") {
        pos <- plot.obj.build$layout$panel_params[[1]]$y.major_source
      } else {
        stop('`axes` must be either "x" or "y"')
      }
      
      return(pos)
    }
    
    
    
    max_nchar_ticks <- function(tick_list) {
      tick_nchars <- sapply(tick_list, nchar)
      return(max(tick_nchars))
    }
    
    plot.custom.dabest <- function(x, ...,
                                   color.column        = NULL,
                                   palette             = "Set1",
                                   float.contrast      = TRUE,
                                   slopegraph          = TRUE,
                                   group.summaries     = "mean_sd",
                                   
                                   rawplot.type        = c("swarmplot", "sinaplot"),
                                   rawplot.ylim        = NULL,
                                   rawplot.ylabel      = NULL,
                                   rawplot.markersize  = 2,
                                   rawplot.groupwidth  = 0.3,
                                   
                                   effsize.ylim        = NULL,
                                   effsize.ylabel       = NULL,
                                   effsize.markersize  = 4,
                                   
                                   theme               = ggplot2::theme_classic(),
                                   tick.fontsize       = 11,
                                   axes.title.fontsize = 14,
                                   
                                   swarmplot.params    = NULL,
                                   sinaplot.params     = NULL,
                                   slopegraph.params   = NULL )
    {
      
      #### Check object class ####
      if (class(x)[1] != "dabest") {
        stop(paste(
          "The object you are plotting is not a `dabest` class object. ",
          "Please check again! ")
        )
      } else {
        dabest.object <- x
      }
      
      
      #### Extract variables ####
      
      # Create handles for easy access to the items in `dabest.object`.
      raw.data           <-  dabest.object$data
      boot.result        <-  dabest.object$result
      idx                <-  dabest.object$idx
      id.col             <-  dabest.object$id.column
      summary            <-  dabest.object$summary
      
      plot.groups.sizes  <-  unlist(lapply(idx, length))
      all.groups         <-  unlist(idx)
      
      # The variables below should are quosures!
      x_enquo            <-  dabest.object$x
      y_enquo            <-  dabest.object$y
      x_quoname          <-  rlang::quo_name(x_enquo)
      y_quoname          <-  rlang::quo_name(y_enquo)
      # `func` is not a quosure but a string.
      func               <-  boot.result$func[1]
      is.paired          <-  boot.result$paired[1]
      
      
      #### Decide if floating or slopegraph. ####
      # float.contrast and slopegraph
      if (isFALSE(is.paired)) slopegraph <- FALSE
      
      if (max(plot.groups.sizes) > 2) {
        float.contrast <- FALSE
        slopegraph     <- FALSE
      }
      
      if (length(all.groups) > 2) {
        float.contrast <- FALSE
      }
      
      
      
      #### Select data for plotting. ####
      if (length(all.groups)     == 2 &&
          plot.groups.sizes[[1]] == 2) {
        # Not multiplot. Add it to an empty list.
        for.plot <- raw.data
        
      } else {
        # Reorder the plot data according to idx.
        for.plot <- list()
        for (subplot_groups in idx) {
          subplot  <- dplyr::filter(raw.data, !!x_enquo %in% subplot_groups)
          for.plot[[length(for.plot) + 1]] <- subplot
        }
        
        for.plot <- dplyr::bind_rows(for.plot)
      }
      
      for.plot[[x_quoname]] <-
        for.plot[[x_quoname]] %>%
        factor(all.groups, ordered = TRUE)
      
      
      #### Compute the Ns. ####
      Ns <-
        for.plot %>%
        dplyr::group_by(!!x_enquo) %>%
        dplyr::count()
      
      Ns$swarmticklabs <-
        do.call(paste, c(Ns[c(x_quoname, "n")],
                         sep = "\nN = "))
      
      
      #### Compute stats for Tufte lines. ####
      for.tufte.lines <-
        for.plot %>%
        dplyr::group_by(!!x_enquo) %>%
        dplyr::summarize(mean           = mean(!!y_enquo),
                         median         = median(!!y_enquo),
                         sd             = sd(!!y_enquo),
                         low.quartile   = stats::quantile(!!y_enquo)[2],
                         upper.quartile = stats::quantile(!!y_enquo)[4]) %>%
        dplyr::mutate(low.sd = mean - sd, upper.sd = mean + sd)
      
      
      #### Parse keywords. ####
      # color.column
      color.col_enquo      <-  rlang::enquo(color.column)
      swarm.dodge        <-  0
      
      if (rlang::quo_is_null(color.col_enquo)) {
        color.aes          <-  ggplot2::aes(col = !!x_enquo)
        # swarm.dodge        <-  0
      } else {
        color.col_quoname  <-  rlang::quo_name(color.col_enquo)
        color.aes          <-  ggplot2::aes(col = !!color.col_enquo)
        # swarm.dodge        <-  0.1
      }
      
      
      # rawplot.type
      # If rawplot is not specified, defaults to 'swarmplot'.
      if (length(rawplot.type) > 1) {
        rawplot.type <- rawplot.type[1]
      }
      
      
      #### swarmplot/sinaplot params. ####
      if (isFALSE(slopegraph)) {
        
        if (rawplot.type == 'swarmplot') {
          if (is.null(swarmplot.params)) {
            swarmplot.params <- list(size        = rawplot.markersize,
                                     width       = rawplot.groupwidth,
                                     dodge.width = swarm.dodge,
                                     mapping     = color.aes,
                                     alpha       = 0.95)
          } else if (class(swarmplot.params) != "list") {
            stop("`swarmplot.params` is not a list.")
          } else swarmplot.params[['mapping']] = color.aes
          
        } else if (rawplot.type == 'sinaplot') {
          swarm.width = 0.3
          if (is.null(sinaplot.params)) {
            sinaplot.params <- list(size = rawplot.markersize,
                                    maxwidth = swarm.width,
                                    mapping = color.aes)
          } else if (class(sinaplot.params) != "list") {
            stop("`sinaplot.params` is not a list.")
          } else sinaplot.params[['mapping']] = color.aes
          
        } else stop(paste(rawplot.type, "is not a recognized plot type. ",
                          "Accepted plot types: 'swarmplot' and 'sinaplot'."))
      } else {
        rawplot.type <- "slopegraph"
      }
      
      # y-axes labels.
      if (is.null(rawplot.ylabel)) {
        rawplot.ylabel <- stringr::str_interp("${y_quoname}\n")
      } else {
        rawplot.ylabel <- stringr::str_interp("${rawplot.ylabel}\n")
      }
      
      if (is.null(effsize.ylabel)) {
        if (isTRUE(is.paired)) {
          effsize.ylabel <-
            stringr::str_interp("Paired ${func} difference\n")
        } else {
          effsize.ylabel <-
            stringr::str_interp("Unpaired ${func} difference\n")
        }
      } else {
        effsize.ylabel <- stringr::str_interp("${effsize.ylabel}\n")
      }
      
      
      
      
      #### Create themes. ####
      horizontal.line.width = 0.4
      
      non.floating.theme <-
        theme +
        ggplot2::theme(
          axis.text            =  ggplot2::element_text(size = tick.fontsize),
          axis.title           =  ggplot2::element_text(size = axes.title.fontsize),
          axis.ticks.length    =  ggplot2::unit(7, "points"),
          axis.ticks.x.bottom  =  ggplot2::element_blank(),
          axis.title.x.bottom  =  ggplot2::element_blank()
        )
      
      floating.theme <-
        non.floating.theme +
        ggplot2::theme(
          axis.title.x.bottom  =  ggplot2::element_blank(),
          axis.ticks.x.bottom  =  ggplot2::element_blank()
        )
      
      legend.theme <-
        ggplot2::theme(
          legend.title         =  ggplot2::element_text(size = axes.title.fontsize),
          legend.text          =  ggplot2::element_text(size = tick.fontsize))
      
      
      non.floating.theme <-  non.floating.theme + legend.theme
      floating.theme     <-  floating.theme + legend.theme
      
      
      remove.axes <-
        ggplot2::theme(
          axis.line.x          = ggplot2::element_blank(),
          axis.title.x         = ggplot2::element_blank(),
          axis.ticks.x.bottom  = ggplot2::element_blank()
        )
      
      
      
      #### Set rawdata plot ylims. ####
      if (is.null(rawplot.ylim)) {
        rawplot.ylim <- range(for.plot[[y_quoname]])
      }
      
      # Equalize the xlims across both plots.
      if (isTRUE(float.contrast)) {
        rawdata.coord_cartesian <-
          ggplot2::coord_cartesian(ylim = rawplot.ylim)
        
      } else {
        both.xlim <- c(1, length(all.groups) + 0.3)
        rawdata.coord_cartesian <-
          ggplot2::coord_cartesian(xlim = both.xlim, ylim = rawplot.ylim)
      }
      
      
      
      
      #### Plot raw data. ####
      # slopegraph.
      if (rawplot.type == "slopegraph") {
        
        rawdata.plot <-
          ggplot2::ggplot() +
          rawdata.coord_cartesian +
          ggplot2::ylab(rawplot.ylabel) +
          ggplot2::scale_x_discrete(labels = Ns$swarmticklabs,
                                    limits = all.groups)
        
        slope.line.width  <- 0.5
        
        for (subplot_groups in idx) {
          subplot  <- dplyr::filter(raw.data, !!x_enquo %in% subplot_groups)
          
          subplot[[x_quoname]] <-
            subplot[[x_quoname]] %>%
            factor(subplot_groups, ordered = TRUE)
          
          if (rlang::quo_is_null(color.col_enquo)) {
            rawdata.plot <-
              rawdata.plot +
              ggplot2::geom_line(data = subplot,
                                 size = slope.line.width,
                                 alpha = 0.8,
                                 ggplot2::aes(!!x_enquo, !!y_enquo,
                                              group = !!id.col)
              )
          } else {
            rawdata.plot <-
              rawdata.plot +
              ggplot2::geom_line(data = subplot,
                                 size = slope.line.width,
                                 alpha = 0.75,
                                 ggplot2::aes(!!x_enquo, !!y_enquo,
                                              group = !!id.col,
                                              colour = !!color.col_enquo)
              )
          }
        }
        
        
        
      } else { # swarmplot.
        rawdata.plot <-
          ggplot2::ggplot(data = for.plot,
                          ggplot2::aes(!!x_enquo, !!y_enquo)) +
          rawdata.coord_cartesian +
          #ggplot2::scale_colour_continuous(low="red",high="blue") +
          #ggplot2::scale_color_distiller(palette=palette) +
          #ggplot2::scale_color_brewer(palette = palette) +
          ggplot2::ylab(rawplot.ylabel) +
          ggplot2::scale_x_discrete(breaks = all.groups,
                                    labels = Ns$swarmticklabs)
        
        if (rawplot.type == 'swarmplot') {
          
          rawdata.plot <-
            rawdata.plot +
            do.call(ggbeeswarm::geom_quasirandom, swarmplot.params)
          
          if (isTRUE(float.contrast)) {
            rawdata.plot <- rawdata.plot + floating.theme
          } else {
            rawdata.plot <- rawdata.plot + non.floating.theme
          }
          
          # if (isTRUE(float.contrast)) {
          #   rawdata.plot <-
          #     rawdata.plot +
          #     do.call(ggbeeswarm::geom_beeswarm, swarmplot.params) +
          #     floating.theme
          # } else {
          # rawdata.plot <-
          #   rawdata.plot +
          #   do.call(ggbeeswarm::geom_quasirandom, swarmplot.params) +
          #   non.floating.theme
          # }
          
          
        } else if (rawplot.type == 'sinaplot') {
          rawdata.plot   <-
            rawdata.plot +
            do.call(ggforce::geom_sina, sinaplot.params)
          
          if (isTRUE(float.contrast)) {
            rawdata.plot <- rawdata.plot + floating.theme
          } else {
            rawdata.plot <- rawdata.plot + non.floating.theme
          }
        }
        
        #### Plot group summaries. ####
        if (isFALSE(float.contrast)) {
          line.nudge <- rawplot.groupwidth * 1.25
          if (line.nudge > 0.8) line.nudge <- 0.8
          pos.nudge = ggplot2::position_nudge(x = line.nudge)
          
          if (!is.null(group.summaries)) {
            accepted.summaries <- c('mean_sd', 'median_quartiles')
            
            not.in.g.summs <- !(group.summaries %in% accepted.summaries)
            
            if (not.in.g.summs) {
              err1 <- stringr::str_interp("${group.summaries} is not a recognized option.")
              err2 <- "Accepted `group.summaries` are 'mean_sd' or 'median_quartiles'."
              stop(paste(err1, err2))
              
            } else if (group.summaries == 'mean_sd') {
              rawdata.plot <-
                rawdata.plot +
                suppressWarnings(
                  ggplot2::geom_linerange(
                    data     = for.tufte.lines,
                    size     = 1,
                    position = pos.nudge,
                    ggplot2::aes(x = !!x_enquo, y = mean,
                                 ymin = low.sd,
                                 ymax = upper.sd)) ) +
                ggplot2::geom_point(
                  data     = for.tufte.lines,
                  size     = 0.75,
                  position = pos.nudge,
                  #colour   = "white",
                  colour="grey",
                  ggplot2::aes(x = !!x_enquo, y = mean))
              
            } else if (group.summaries == 'median_quartiles') {
              rawdata.plot <-
                rawdata.plot +
                suppressWarnings(
                  ggplot2::geom_linerange(
                    data     = for.tufte.lines,
                    size     = 1,
                    position = pos.nudge,
                    ggplot2::aes(x = !!x_enquo, y = median,
                                 ymin = low.quartile,
                                 ymax = upper.quartile)) ) +
                ggplot2::geom_point(
                  data     = for.tufte.lines,
                  size     = 0.75,
                  position = pos.nudge,
                  #colour   = "white",
                  colour="grey",
                  ggplot2::aes(x = !!x_enquo, y = median))
            }
          }
        }
      }
      
      # Plot the summary lines for each group if `float.contrast` is TRUE.
      if (isTRUE(float.contrast)) {
        func_control <- summary[[func]][1]
        func_test    <- summary[[func]][2]
        
        rawdata.plot <- rawdata.plot +
          # Plot the summary lines for the control group...
          ggplot2::geom_segment(
            #color = "white",
            color="grey",
            size  = horizontal.line.width,
            ggplot2::aes(x    = 1,
                         xend = 3,
                         y    = func_control,
                         yend = func_control)) +
          
          # ... and the test group.
          ggplot2::geom_segment(
            #color = "white",
            color="grey",
            size  = horizontal.line.width,
            ggplot2::aes(x    = 2,
                         xend = 3,
                         y    = func_test,
                         yend = func_test))
        
        # Apply appropriate theme to swarm plot.
        rawdata.plot <- rawdata.plot + floating.theme
      } else {
        rawdata.plot <- rawdata.plot + non.floating.theme
      }
      
      
      
      #### Munge bootstraps. ####
      # Munge bootstraps into tibble for easy plotting with ggplot.
      boots.for.plot <- tibble::as_tibble(data.frame(boot.result$bootstraps))
      colnames(boots.for.plot) <- boot.result$test_group
      
      if (isFALSE(float.contrast)) {
        # Add the control group as a set of NaNs.
        for (control.column in unique(boot.result$control_group)) {
          oldcols <- colnames(boots.for.plot)
          boots.for.plot <-
            boots.for.plot %>%
            tibble::add_column(placeholder = rep(NaN, nrow(boots.for.plot)))
          
          colnames(boots.for.plot) <- c(oldcols, control.column)
        }
      }
      
      boots.for.plot <-
        tidyr::gather(boots.for.plot, !!x_enquo, !!y_enquo)
      
      # Order the bootstraps so they plot in the correct order.
      boots.for.plot[[x_quoname]] <-
        factor(boots.for.plot[[x_quoname]], all.groups, ordered = TRUE)
      
      boots.for.plot <-
        dplyr::arrange(boots.for.plot, !!x_enquo)
      
      
      
      #### Set delta plot ylims. ####
      if (is.null(effsize.ylim)) {
        effsize.ylim <- range( na.omit(boots.for.plot[y_quoname]) )
      }
      
      
      
      #### Plot bootstraps. ####
      float.reflines.xstart <- 0.4
      float.reflines.xend   <- 1.6
      
      if (isTRUE(float.contrast)) {
        es0.trimming        <- 0
        flat.violin.width   <- 1
        flat.violin.adjust  <- 5
        x.start             <- float.reflines.xstart
        x.end               <- float.reflines.xend
        
      } else {
        es0.trimming        <- 0.5
        flat.violin.width   <- 0.75
        flat.violin.adjust  <- 3
        x.start             <- 0
        x.end               <- length(all.groups) + es0.trimming
      }
      
      delta.plot <-
        ggplot2::ggplot(boots.for.plot, na.rm = TRUE) +
        geom_flat_violin(
          ggplot2::aes(!!x_enquo, !!y_enquo),
          na.rm  =  TRUE,
          width  =  flat.violin.width,
          adjust =  flat.violin.adjust,
          #fill = "red",
          fill="coral",
          size   =  0 # width of the line.
        ) +
        # This is the line representing the null effect size.
        ggplot2::geom_segment(
          #color  =  "white",
          color="grey",
          size   =  horizontal.line.width,
          x      =  x.start,
          xend   =  x.end,
          y      =  0,
          yend   =  0)
      
      
      
      #### Plot effect sizes and CIs. ####
      
      delta.plot <-
        delta.plot +
        ggplot2::ylab(effsize.ylabel) +
        ggplot2::geom_point(
          data  = boot.result,
          #color = "white",
          color="grey",
          size  = effsize.markersize,
          ggplot2::aes(test_group, difference)) +
        ggplot2::geom_errorbar(
          data  = boot.result,
          #color = "white",
          color="grey",
          width = 0,
          size  = 0.75,
          ggplot2::aes(x    = test_group,
                       ymin = bca_ci_low,
                       ymax = bca_ci_high))
      
      
      
      #### Float vs nonfloat delta plots. ####
      if (isTRUE(float.contrast)) {
        
        # Shift ylims appropriately.
        if (func_control > 0) {
          new.delta.ylim <- rawplot.ylim - func_control
        } else {
          new.delta.ylim <- rawplot.ylim + func_control
        }
        
        delta.plot <-
          delta.plot +
          ggplot2::coord_cartesian(ylim = new.delta.ylim) +
          ggplot2::scale_y_continuous(position = "right") +
          # This is the delta-side effect size line,
          # that aligns with the central measure of the test group.
          ggplot2::geom_segment(#color = "white",
                                size  = horizontal.line.width,
                                x     = float.reflines.xstart,
                                xend  = float.reflines.xend,
                                y     = boot.result$difference[1],
                                yend  = boot.result$difference[1]) +
          ggplot2::scale_x_discrete(labels =
                                      c(stringr::str_interp("${all.groups[2]}\nminus ${all.groups[1]}")) ) +
          floating.theme
        
        
      } else {
        # Plot nonfloating deltas.
        # Properly concatenate the delta.plot labels.
        delta.tick.labs  <- vector("list", length(idx))
        i <- 1
        
        for (subplot_groups in idx) {
          control_group <- subplot_groups[1]
          test_groups   <- subplot_groups[2: length(subplot_groups)]
          
          labels <- c(" ",
                      paste(test_groups, stringr::str_interp("minus\n${control_group}"),
                            sep = "\n"))
          
          delta.tick.labs[[i]] = labels
          i <- i + 1
        }
        
        # Equalize the xlims across both plots, and set ylims for deltaplot.
        delta.plot <- delta.plot +
          ggplot2::coord_cartesian(xlim = both.xlim,
                                   ylim = effsize.ylim) +
          ggplot2::scale_x_discrete(breaks = all.groups,
                                    labels = delta.tick.labs) +
          non.floating.theme
      }
      
      
      
      #### Trim rawdata axes. ####
      rawdata.plot <-  rawdata.plot + remove.axes
      
      # Get the ylims.
      rawdata.plot.build       <- ggplot2::ggplot_build(rawdata.plot)
      rawdata.plot.build.panel <- rawdata.plot.build$layout$panel_params[[1]]
      rawdata.plot.ylim        <- rawdata.plot.build.panel$y.range
      
      segment.ypos             <- rawdata.plot.ylim[1]
      
      rawdata.plot.xlim        <- rawdata.plot.build.panel$x.range
      rawdata.plot.lower.xlim  <- rawdata.plot.xlim[1]
      
      
      # Set padding to add.
      start.idx          <- 1
      padding            <- 0.25
      
      # Re-draw the trimmed axes.
      for (size in plot.groups.sizes) {
        end.idx      <- start.idx + size - 1
        
        if (isTRUE(float.contrast)) {
          xstart   <- rawdata.plot.lower.xlim
        } else {
          xstart   <- start.idx - padding
        }
        
        if (isTRUE(slopegraph)) {
          
          rawdata.plot <- rawdata.plot +
            ggplot2::geom_segment(
              # size = segment.thickness,
              ggplot2::aes_(x    = xstart,
                            xend = end.idx + padding,
                            y    = segment.ypos,
                            yend = segment.ypos))
          
        } else {
          
          if (isTRUE(float.contrast)) {
            xend   <- end.idx + padding * 1.5
          } else {
            xend   <- end.idx + padding
          }
          
          rawdata.plot <- rawdata.plot +
            ggplot2::geom_segment(
              x    = xstart,
              xend = xend,
              y    = segment.ypos,
              yend = segment.ypos#,color="white"
              )
          
        }
        
        start.idx  <- start.idx + size
      }
      
      
      
      
      #### Trim deltaplot axes. ####
      delta.plot  <-  delta.plot + remove.axes
      
      # Get the ylims.
      delta.plot.build         <- ggplot2::ggplot_build(delta.plot)
      delta.plot.build.panel   <- delta.plot.build$layout$panel_params[[1]]
      delta.plot.ylim          <- delta.plot.build.panel$y.range
      segment.ypos             <- delta.plot.ylim[1]
      delta.plot.upper.ylim    <- delta.plot.ylim[2]
      
      delta.plot.xlim          <- delta.plot.build.panel$x.range
      delta.plot.lower.xlim    <- delta.plot.xlim[1]
      delta.plot.upper.xlim    <- delta.plot.xlim[2]
      
      # Set padding to add.
      start.idx      <- 1
      
      # Re-draw the trimmed axes.
      for (size in plot.groups.sizes) {
        end.idx      <- start.idx + size - 1
        
        if (isTRUE(float.contrast)) {
          xstart     <- delta.plot.lower.xlim
          xend       <- delta.plot.upper.xlim
          
        } else {
          xstart     <- start.idx - padding
          xend       <- end.idx   + padding
        }
        
        delta.plot <-
          delta.plot +
          ggplot2::geom_segment(x    = xstart,
                                xend = xend,
                                y    = segment.ypos,
                                yend = segment.ypos#,color="white"
          )
        
        start.idx  <- start.idx + size
      }
      
      
      
      #### Handle color legend. ####
      if (!rlang::quo_is_null(color.col_enquo)) {
        legend <- cowplot::get_legend(rawdata.plot)
      }
      # Remove the legend from the rawplot.
      rawdata.plot <- rawdata.plot + ggplot2::guides(color = "none")
      
      
      
      #### Equalize tick label lengths. ####
      if (isFALSE(float.contrast)) {
        rawplot.yticks.labels    <- get_tick_labels(rawdata.plot, axes="y")
        rawplot.yticks.breaks    <- as.numeric(rawplot.yticks.labels)
        max_rawplot_ticklength   <- max_nchar_ticks(rawplot.yticks.labels)
        
        deltaplot.yticks.labels  <- get_tick_labels(delta.plot, axes="y")
        deltaplot.yticks.breaks  <- as.numeric(deltaplot.yticks.labels)
        max_deltaplot_ticklength <- max_nchar_ticks(deltaplot.yticks.labels)
        
        
        if (max_rawplot_ticklength < max_deltaplot_ticklength) {
          space.diff <- max_deltaplot_ticklength - max_rawplot_ticklength
          
          suffix.spacing        <- rep(" ", space.diff)
          
          rawplot.yticks.labels <- paste(stringr::str_interp(suffix.spacing),
                                         rawplot.yticks.labels)
          rawdata.plot <-
            rawdata.plot +
            ggplot2::scale_y_continuous(breaks = rawplot.yticks.breaks,
                                        labels = rawplot.yticks.labels)
          
          
        } else if (max_rawplot_ticklength > max_deltaplot_ticklength) {
          space.diff = max_rawplot_ticklength - max_deltaplot_ticklength
          
          suffix.spacing          <- rep(" ", space.diff)
          
          deltaplot.yticks.labels <- paste(stringr::str_interp(suffix.spacing),
                                           deltaplot.yticks.labels)
          
          delta.plot <- delta.plot +
            ggplot2::scale_y_continuous(breaks = deltaplot.yticks.breaks,
                                        labels = deltaplot.yticks.labels)
        }
        
      }
      
      
      
      
      #### Determine layout. ####
      if (isTRUE(float.contrast)) {
        # Side-by-side floating plot layout.
        # plot.margin declares the top, right, bottom, left margins in order.
        rawdata.plot <-
          rawdata.plot +
          ggplot2::theme(
            plot.margin = ggplot2::unit(c(5.5, 0, 5.5, 5.5), "pt")
          )
        
        delta.plot   <-
          delta.plot +
          ggplot2::theme(
            plot.margin        = ggplot2::unit(c(5.5, 5.5, 5.5, 0), "pt"),
            axis.line.x.bottom = ggplot2::element_blank()
          )
        
        aligned_spine = 'b'
        nrows <- 1
        
        if (rlang::quo_is_null(color.col_enquo)) {
          plist <- list(rawdata.plot, delta.plot)
          ncols <- 2
          widths <- c(0.7, 0.3)
        } else {
          plist <- list(rawdata.plot, delta.plot, legend)
          ncols <- 3
          widths <- c(0.7, 0.3, 0.2)
        }
        
      } else {
        # Above-below non-floating plot layout.
        aligned_spine = 'lr'
        nrows <- 2
        
        if (rlang::quo_is_null(color.col_enquo)) {
          plist <- list(rawdata.plot, delta.plot)
          ncols <- 1
          widths <- c(1)
        } else {
          plist <- list(rawdata.plot, legend, delta.plot)
          ncols <- 2
          widths <- c(0.85, 0.15)
        }
      }
      
      
      result <- cowplot::plot_grid(
        plotlist   = plist,
        nrow       = nrows,
        ncol       = ncols,
        rel_widths = widths,
        axis       = aligned_spine)
      return(result)
      
    }
    
    geom_flat_violin <- function(
      mapping     = NULL,
      data        = NULL,
      stat        = "ydensity",
      position    = "dodge",
      trim        = TRUE,
      scale       = "area",
      show.legend = NA,
      inherit.aes = TRUE, ...) {
      
      ggplot2::layer(
        data        = data,
        mapping     = mapping,
        stat        = stat,
        geom        = GeomFlatViolin,
        position    = position,
        show.legend = show.legend,
        inherit.aes = inherit.aes,
        params      = list(trim = trim,scale = scale, ...)
      )
      
    }
    
    
    
    "%||%" <- function(a, b) {
      if (!is.null(a)) a else b
    }
    
    
    
    GeomFlatViolin <-
      ggplot2::ggproto(
        "GeomFlatViolin",
        
        ggplot2::Geom,
        
        setup_data = function(data, params) {
          data$width <- data$width %||%
            params$width %||%
            (resolution(data$x, FALSE) * 0.9)
          # ymin, ymax, xmin, and xmax define the bounding rectangle for each group.
          data %>%
            dplyr::group_by(group) %>%
            dplyr::mutate(ymin = min(y),
                          ymax = max(y),
                          xmin = x,
                          xmax = x + width / 2)
        },
        
        draw_group = function(data, panel_scales, coord) {
          # Find the points for the line to go all the way around
          data <- transform(data, xminv = x,
                            xmaxv = x + violinwidth * (xmax - x))
          
          # Make sure it's sorted properly to draw the outline
          newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
                           plyr::arrange(transform(data, x = xmaxv), -y))
          
          # Close the polygon: set first and last point the same.
          # Needed for coord_polar and such.
          newdata <- rbind(newdata, newdata[1,])
          
          ggplot2:::ggname("geom_flat_violin",
                           ggplot2::GeomPolygon$draw_panel(
                             newdata, panel_scales, coord)
          )
        },
        
        draw_key = ggplot2::draw_key_polygon,
        
        default_aes = ggplot2::aes(weight = 1, colour = "grey20",
                                   fill = "grey", size = 0.5,
                                   alpha = NA, linetype = "solid"),
        
        required_aes = c("x", "y")
      )
  }
  # cris = agg[agg$SPP=="CRISSALE",]
  # cris = cris[,c("TAIL","WHICH.SIDE.OF.CFB")]
  # cris$GROUP = substr(cris$WHICH.SIDE.OF.CFB,1,3)
  # cris = cris[complete.cases(cris),]
  # bav = dabest(cris,GROUP,TAIL, idx=(unique(cris$GROUP)))
  
  theme_transparentwhite <- function (base_size = 12, base_family = "") {
    theme_classic(base_size = base_size, base_family = base_family) +
      theme(
        axis.line = element_line(color="white"),
        axis.line.x = element_line(color="white"),
        axis.line.x.bottom = element_line(color="white"),
        axis.line.x.top = element_line(color="white"),
        axis.line.y = element_line(color="white"),
        axis.line.y.left = element_line(color="white"),
        axis.line.y.right = element_line(color="white"),
        axis.text = element_text(color = "white"),
        axis.text.x=element_text(colour="white"),
        axis.text.y=element_text(colour="white"),
        axis.ticks = element_line(color="white"),
        axis.title = element_text(color = "white"),
        axis.title.x = element_text(color="white"),
        axis.title.y = element_text(color="white"),
        axis.text.x.bottom = element_text(color="white"),
        axis.text.x.top = element_text(color="white"),
        axis.text.y.left = element_text(color="white"),
        axis.text.y.right = element_text(color="white"),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        legend.direction = 'horizontal', 
        legend.key = element_rect(fill="transparent",colour=NA),
        legend.position = 'bottom',
        legend.text = element_text(colour = 'white'),
        line = element_line(color="white"),
        panel.background = element_rect(fill = "transparent", colour = NA),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.border = element_rect(fill = "transparent", colour = NA),
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        plot.title = element_text(color="white", face="italic",hjust=0),
        rect=element_rect(color="white"),
        text=element_text(color="white")
      )
  } ;
  # g = plot.custom.dabest(x=bav,palette="Set1",group.summaries="mean_sd",
  #      rawplot.ylabel="Beak Surface Area Volume Ratio",
  #      theme=theme_transparentwhite())
  # #print(g)
  # ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/DABEST_CRISSALE_TAIL.png",
  #        g,bg="transparent",
  #        width=6.5,height=5,units="in")
  
  ## all the different variables
  {
    print("TAIL"); {
      temp = agg
      temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
      temp = temp[,c("GROUP","TAIL")]
      temp = temp[complete.cases(temp), ]
      bav = dabest(temp,GROUP,TAIL, idx=(list(c("BIL_SON","BIL_CHI"),
                                              c("BRU_SON","BRU_CHI"),
                                              c("MEL_SON","MEL_CHI"),
                                              c("NIT_SON","NIT_CHI"),
                                              c("SIN_SON","SIN_CHI"),
                                              c("BEL_SON","BEL_CHI"),
                                              c("CRI_SON","CRI_CHI"),
                                              c("CUR_SON","CUR_CHI"),
                                              c("FLA_SON","FLA_CHI"),
                                              c("FUS_SON","FUS_CHI"))))
      g=plot(x=bav,group.summaries=NULL,
             rawplot.ylabel="TAIL",
             #theme=theme_transparentwhite(),
             effsize.markersize=2)
      
      # g = plot.custom.dabest(x=bav,palette=c("red","blue"),group.summaries=NULL,
      #                        rawplot.ylabel="TAIL",
      #                        #theme=theme_transparentwhite(),
      #                        effsize.markersize=2)
      print(g)
      ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/DABEST/DABEST_TAIL.png",
             g,bg="transparent",width=23,height=7,units="in")
      
      temp = agg_res
      temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
      temp = temp[,c("GROUP","RES_TAIL")]
      temp = temp[complete.cases(temp), ]
      bav = dabest(temp,GROUP,RES_TAIL, idx=(list(c("BIL_SON","BIL_CHI"),
                                                  c("BRU_SON","BRU_CHI"),
                                                  c("MEL_SON","MEL_CHI"),
                                                  c("NIT_SON","NIT_CHI"),
                                                  c("SIN_SON","SIN_CHI"),
                                                  c("BEL_SON","BEL_CHI"),
                                                  c("CRI_SON","CRI_CHI"),
                                                  c("CUR_SON","CUR_CHI"),
                                                  c("FLA_SON","FLA_CHI"),
                                                  c("FUS_SON","FUS_CHI"))))
      g = plot.custom.dabest(x=bav,palette="Paired",group.summaries=NULL,
                             rawplot.ylabel="RES_TAIL",
                             #theme=theme_transparentwhite(),
                             effsize.markersize=2)
      #print(g)
      ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/DABEST/DABEST_RES_TAIL.png",
             g,bg="transparent",width=23,height=7,units="in")
    }
    print("BEAKAREAVOLUME_CALC"); {
      temp = agg
      temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
      temp = temp[,c("GROUP","BEAKAREAVOLUME_CALC")]
      temp = temp[complete.cases(temp), ]
      bav = dabest(temp,GROUP,BEAKAREAVOLUME_CALC, idx=(list(c("BIL_SON","BIL_CHI"),
                                                             c("BRU_SON","BRU_CHI"),
                                                             c("MEL_SON","MEL_CHI"),
                                                             c("NIT_SON","NIT_CHI"),
                                                             c("SIN_SON","SIN_CHI"),
                                                             c("BEL_SON","BEL_CHI"),
                                                             c("CRI_SON","CRI_CHI"),
                                                             c("CUR_SON","CUR_CHI"),
                                                             c("FLA_SON","FLA_CHI"),
                                                             c("FUS_SON","FUS_CHI"))))
      g = plot.custom.dabest(x=bav,palette=c("red","blue"),group.summaries=NULL,
                             rawplot.ylabel="BEAKAREAVOLUME_CALC",
                             #theme=theme_transparentwhite(),
                             effsize.markersize=2)
      print(g)
      ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/DABEST/DABEST_BEAKAREAVOLUME_CALC.png",
             g,bg="transparent",width=23,height=7,units="in")
      
      temp = agg_res
      temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
      temp = temp[,c("GROUP","RES_BEAKAREAVOLUME_CALC")]
      temp = temp[complete.cases(temp), ]
      bav = dabest(temp,GROUP,RES_BEAKAREAVOLUME_CALC, idx=(list(c("BIL_SON","BIL_CHI"),
                                                                 c("BRU_SON","BRU_CHI"),
                                                                 c("MEL_SON","MEL_CHI"),
                                                                 c("NIT_SON","NIT_CHI"),
                                                                 c("SIN_SON","SIN_CHI"),c("BEL_SON","BEL_CHI"),
                                                                 c("CRI_SON","CRI_CHI"),
                                                                 c("CUR_SON","CUR_CHI"),
                                                                 c("FLA_SON","FLA_CHI"),
                                                                 c("FUS_SON","FUS_CHI"))))
      g = plot.custom.dabest(x=bav,palette="Paired",group.summaries=NULL,
                             rawplot.ylabel="RES_BEAKAREAVOLUME_CALC",
                             #theme=theme_transparentwhite(),
                             effsize.markersize=2)
      #print(g)
      ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/DABEST/DABEST_RES_BEAKAREAVOLUME_CALC.png",
             g,bg="transparent",width=23,height=7,units="in")
      
    }
    print("LAT"); {
      
      
      temp = agg
      temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
      temp = temp[,c("GROUP","LAT")]
      temp = temp[complete.cases(temp), ]
      bav = dabest(temp,GROUP,LAT, idx=(list(c("BIL_SON","BIL_CHI"),
                                             c("BRU_SON","BRU_CHI"),
                                             c("MEL_SON","MEL_CHI"),
                                             c("NIT_SON","NIT_CHI"),
                                             c("SIN_SON","SIN_CHI"),c("BEL_SON","BEL_CHI"),
                                             c("CRI_SON","CRI_CHI"),
                                             c("CUR_SON","CUR_CHI"),
                                             c("FLA_SON","FLA_CHI"),
                                             c("FUS_SON","FUS_CHI"))))
      g = plot.custom.dabest(x=bav,palette=c("red","blue"),group.summaries=NULL,
                             rawplot.ylabel="LAT",
                             #theme=theme_transparentwhite(),
                             effsize.markersize=2)
      print(g)
      ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/DABEST/DABEST_LAT.png",
             g,bg="transparent",width=23,height=7,units="in")
      
      
    }
    print("WING.LENGTH.TO.SECONDARIES"); {
      temp = agg
      temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
      temp = temp[,c("GROUP","WING.LENGTH.TO.SECONDARIES")]
      temp = temp[complete.cases(temp), ]
      bav = dabest(temp,GROUP,WING.LENGTH.TO.SECONDARIES, idx=(list(c("BIL_SON","BIL_CHI"),
                                                                    c("BRU_SON","BRU_CHI"),
                                                                    c("MEL_SON","MEL_CHI"),
                                                                    c("NIT_SON","NIT_CHI"),
                                                                    c("SIN_SON","SIN_CHI"),c("BEL_SON","BEL_CHI"),
                                                                    c("CRI_SON","CRI_CHI"),
                                                                    c("CUR_SON","CUR_CHI"),
                                                                    c("FLA_SON","FLA_CHI"),
                                                                    c("FUS_SON","FUS_CHI"))))
      g = plot.custom.dabest(x=bav,palette=c("red","blue"),group.summaries=NULL,
                             rawplot.ylabel="WING.LENGTH.TO.SECONDARIES",
                             #theme=theme_transparentwhite(),
                             effsize.markersize=2)
      print(g)
      ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/DABEST/DABEST_WING.LENGTH.TO.SECONDARIES.png",
             g,bg="transparent",width=23,height=7,units="in")
      
      temp = agg_res
      temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
      temp = temp[,c("GROUP","RES_WING.LENGTH.TO.SECONDARIES")]
      temp = temp[complete.cases(temp), ]
      bav = dabest(temp,GROUP,RES_WING.LENGTH.TO.SECONDARIES, idx=(list(c("BIL_SON","BIL_CHI"),
                                                                        c("BRU_SON","BRU_CHI"),
                                                                        c("MEL_SON","MEL_CHI"),
                                                                        c("NIT_SON","NIT_CHI"),
                                                                        c("SIN_SON","SIN_CHI"),c("BEL_SON","BEL_CHI"),
                                                                        c("CRI_SON","CRI_CHI"),
                                                                        c("CUR_SON","CUR_CHI"),
                                                                        c("FLA_SON","FLA_CHI"),
                                                                        c("FUS_SON","FUS_CHI"))))
      g = plot.custom.dabest(x=bav,palette="Paired",group.summaries=NULL,
                             rawplot.ylabel="RES_WING.LENGTH.TO.SECONDARIES",
                             #theme=theme_transparentwhite(),
                             effsize.markersize=2)
      #print(g)
      ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/DABEST/DABEST_RES_WING.LENGTH.TO.SECONDARIES.png",
             g,bg="transparent",width=23,height=7,units="in")
      
    }
    print("BEAKVOL_CALC"); {
      temp = agg
      temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
      temp = temp[,c("GROUP","BEAKVOL_CALC")]
      temp = temp[complete.cases(temp), ]
      bav = dabest(temp,GROUP,BEAKVOL_CALC, idx=(list(c("BIL_SON","BIL_CHI"),
                                                      c("BRU_SON","BRU_CHI"),
                                                      c("MEL_SON","MEL_CHI"),
                                                      c("NIT_SON","NIT_CHI"),
                                                      c("SIN_SON","SIN_CHI"),c("BEL_SON","BEL_CHI"),
                                                      c("CRI_SON","CRI_CHI"),
                                                      c("CUR_SON","CUR_CHI"),
                                                      c("FLA_SON","FLA_CHI"),
                                                      c("FUS_SON","FUS_CHI"))))
      g = plot.custom.dabest(x=bav,palette=c("red","blue"),group.summaries=NULL,
                             rawplot.ylabel="BEAKVOL_CALC",
                             #theme=theme_transparentwhite(),
                             effsize.markersize=2)
      print(g)
      ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/DABEST/DABEST_BEAKVOL_CALC.png",
             g,bg="transparent",width=23,height=7,units="in")
      
      temp = agg_res
      temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
      temp = temp[,c("GROUP","RES_BEAKVOL_CALC")]
      temp = temp[complete.cases(temp), ]
      bav = dabest(temp,GROUP,RES_BEAKVOL_CALC, idx=(list(c("BIL_SON","BIL_CHI"),
                                                          c("BRU_SON","BRU_CHI"),
                                                          c("MEL_SON","MEL_CHI"),
                                                          c("NIT_SON","NIT_CHI"),
                                                          c("SIN_SON","SIN_CHI"),c("BEL_SON","BEL_CHI"),
                                                          c("CRI_SON","CRI_CHI"),
                                                          c("CUR_SON","CUR_CHI"),
                                                          c("FLA_SON","FLA_CHI"),
                                                          c("FUS_SON","FUS_CHI"))))
      g = plot.custom.dabest(x=bav,palette="Paired",group.summaries=NULL,
                             rawplot.ylabel="RES_BEAKVOL_CALC",
                             #theme=theme_transparentwhite(),
                             effsize.markersize=2)
      #print(g)
      ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/DABEST/DABEST_RES_BEAKVOL_CALC.png",
             g,bg="transparent",width=23,height=7,units="in")
      
    }
    print("TARSUS.LENGTH"); {
      temp = agg
      temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
      temp = temp[,c("GROUP","TARSUS.LENGTH")]
      temp = temp[complete.cases(temp), ]
      bav = dabest(temp,GROUP,TARSUS.LENGTH, idx=(list(c("BIL_SON","BIL_CHI"),
                                                       c("BRU_SON","BRU_CHI"),
                                                       c("MEL_SON","MEL_CHI"),
                                                       c("NIT_SON","NIT_CHI"),
                                                       c("SIN_SON","SIN_CHI"),c("BEL_SON","BEL_CHI"),
                                                       c("CRI_SON","CRI_CHI"),
                                                       c("CUR_SON","CUR_CHI"),
                                                       c("FLA_SON","FLA_CHI"),
                                                       c("FUS_SON","FUS_CHI"))))
      g = plot.custom.dabest(x=bav,palette=c("red","blue"),group.summaries=NULL,
                             rawplot.ylabel="TARSUS.LENGTH",
                             #theme=theme_transparentwhite(),
                             effsize.markersize=2)
      print(g)
      ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/DABEST/DABEST_TARSUS.LENGTH.png",
             g,bg="transparent",width=23,height=7,units="in")
      
      temp = agg_res
      temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
      temp = temp[,c("GROUP","RES_TARSUS.LENGTH")]
      temp = temp[complete.cases(temp), ]
      bav = dabest(temp,GROUP,RES_TARSUS.LENGTH, idx=(list(c("BIL_SON","BIL_CHI"),
                                                           c("BRU_SON","BRU_CHI"),
                                                           c("MEL_SON","MEL_CHI"),
                                                           c("NIT_SON","NIT_CHI"),
                                                           c("SIN_SON","SIN_CHI"),c("BEL_SON","BEL_CHI"),
                                                           c("CRI_SON","CRI_CHI"),
                                                           c("CUR_SON","CUR_CHI"),
                                                           c("FLA_SON","FLA_CHI"),
                                                           c("FUS_SON","FUS_CHI"))))
      g = plot.custom.dabest(x=bav,palette="Paired",group.summaries=NULL,
                             rawplot.ylabel="RES_TARSUS.LENGTH",
                             #theme=theme_transparentwhite(),
                             effsize.markersize=2)
      #print(g)
      ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/DABEST/DABEST_RES_TARSUS.LENGTH.png",
             g,bg="transparent",width=23,height=7,units="in")
      
    }
    print("BILL.LENGTH"); {
      temp = agg
      temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
      temp = temp[,c("GROUP","BILL.LENGTH")]
      temp = temp[complete.cases(temp), ]
      bav = dabest(temp,GROUP,BILL.LENGTH, idx=(list(c("BIL_SON","BIL_CHI"),
                                                     c("BRU_SON","BRU_CHI"),
                                                     c("MEL_SON","MEL_CHI"),
                                                     c("NIT_SON","NIT_CHI"),
                                                     c("SIN_SON","SIN_CHI"),c("BEL_SON","BEL_CHI"),
                                                     c("CRI_SON","CRI_CHI"),
                                                     c("CUR_SON","CUR_CHI"),
                                                     c("FLA_SON","FLA_CHI"),
                                                     c("FUS_SON","FUS_CHI"))))
      g = plot.custom.dabest(x=bav,palette=c("red","blue"),group.summaries=NULL,
                             rawplot.ylabel="BILL.LENGTH",
                             #theme=theme_transparentwhite(),
                             effsize.markersize=2)
      print(g)
      ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/DABEST/DABEST_BILL.LENGTH.png",
             g,bg="transparent",width=23,height=7,units="in")
      
      temp = agg_res
      temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
      temp = temp[,c("GROUP","RES_BILL.LENGTH")]
      temp = temp[complete.cases(temp), ]
      bav = dabest(temp,GROUP,RES_BILL.LENGTH, idx=(list(c("BIL_SON","BIL_CHI"),
                                                         c("BRU_SON","BRU_CHI"),
                                                         c("MEL_SON","MEL_CHI"),
                                                         c("NIT_SON","NIT_CHI"),
                                                         c("SIN_SON","SIN_CHI"),c("BEL_SON","BEL_CHI"),
                                                         c("CRI_SON","CRI_CHI"),
                                                         c("CUR_SON","CUR_CHI"),
                                                         c("FLA_SON","FLA_CHI"),
                                                         c("FUS_SON","FUS_CHI"))))
      g = plot.custom.dabest(x=bav,palette="Paired",group.summaries=NULL,
                             rawplot.ylabel="RES_BILL.LENGTH",
                             #theme=theme_transparentwhite(),
                             effsize.markersize=2)
      #print(g)
      ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/DABEST/DABEST_RES_BILL.LENGTH.png",
             g,bg="transparent",width=23,height=7,units="in")
      
    }
    print("BEAKTOTAREAVOLUME_CALC"); {
      
      
      temp = agg
      temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
      temp = temp[,c("GROUP","BEAKTOTAREAVOLUME_CALC")]
      temp = temp[complete.cases(temp), ]
      bav = dabest(temp,GROUP,BEAKTOTAREAVOLUME_CALC, idx=(list(c("BIL_SON","BIL_CHI"),
                                                                c("BRU_SON","BRU_CHI"),
                                                                c("MEL_SON","MEL_CHI"),
                                                                c("NIT_SON","NIT_CHI"),
                                                                c("SIN_SON","SIN_CHI"),c("BEL_SON","BEL_CHI"),
                                                                c("CRI_SON","CRI_CHI"),
                                                                c("CUR_SON","CUR_CHI"),
                                                                c("FLA_SON","FLA_CHI"),
                                                                c("FUS_SON","FUS_CHI"))))
      g = plot.custom.dabest(x=bav,palette=c("red","blue"),group.summaries=NULL,
                             rawplot.ylabel="BEAKTOTAREAVOLUME_CALC",
                             #theme=theme_transparentwhite(),
                             effsize.markersize=2)
      print(g)
      ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/DABEST/DABEST_BEAKTOTAREAVOLUME_CALC.png",
             g,bg="transparent",width=23,height=7,units="in")
      
      temp = agg_res
      temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
      temp = temp[,c("GROUP","RES_BEAKTOTAREAVOLUME_CALC")]
      temp = temp[complete.cases(temp), ]
      bav = dabest(temp,GROUP,RES_BEAKTOTAREAVOLUME_CALC, idx=(list(c("BIL_SON","BIL_CHI"),
                                                                    c("BRU_SON","BRU_CHI"),
                                                                    c("MEL_SON","MEL_CHI"),
                                                                    c("NIT_SON","NIT_CHI"),
                                                                    c("SIN_SON","SIN_CHI"),c("BEL_SON","BEL_CHI"),
                                                                    c("CRI_SON","CRI_CHI"),
                                                                    c("CUR_SON","CUR_CHI"),
                                                                    c("FLA_SON","FLA_CHI"),
                                                                    c("FUS_SON","FUS_CHI"))))
      g = plot.custom.dabest(x=bav,palette="Paired",group.summaries=NULL,
                             rawplot.ylabel="RES_BEAKTOTAREAVOLUME_CALC",
                             #theme=theme_transparentwhite(),
                             effsize.markersize=2)
      #print(g)
      ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/DABEST/DABEST_RES_BEAKTOTAREAVOLUME_CALC.png",
             g,bg="transparent",width=23,height=7,units="in")
      
    }
    print("KIPPSINDEX_CALC"); {
      
      temp = agg
      temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
      temp = temp[,c("GROUP","KIPPSINDEX_CALC")]
      temp = temp[complete.cases(temp), ]
      bav = dabest(temp,GROUP,KIPPSINDEX_CALC, idx=(list(c("BIL_SON","BIL_CHI"),
                                                         c("BRU_SON","BRU_CHI"),
                                                         c("MEL_SON","MEL_CHI"),
                                                         c("NIT_SON","NIT_CHI"),
                                                         c("SIN_SON","SIN_CHI"),c("BEL_SON","BEL_CHI"),
                                                         c("CRI_SON","CRI_CHI"),
                                                         c("CUR_SON","CUR_CHI"),
                                                         c("FLA_SON","FLA_CHI"),
                                                         c("FUS_SON","FUS_CHI"))))
      g = plot.custom.dabest(x=bav,palette=c("red","blue"),group.summaries=NULL,
                             rawplot.ylabel="KIPPSINDEX_CALC",
                             #theme=theme_transparentwhite(),
                             effsize.markersize=2)
      print(g)
      ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/DABEST/DABEST_KIPPSINDEX_CALC.png",
             g,bg="transparent",width=23,height=7,units="in")
      
      temp = agg_res
      temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
      temp = temp[,c("GROUP","RES_KIPPSINDEX_CALC")]
      temp = temp[complete.cases(temp), ]
      bav = dabest(temp,GROUP,RES_KIPPSINDEX_CALC, idx=(list(c("BIL_SON","BIL_CHI"),
                                                             c("BRU_SON","BRU_CHI"),
                                                             c("MEL_SON","MEL_CHI"),
                                                             c("NIT_SON","NIT_CHI"),
                                                             c("SIN_SON","SIN_CHI"),c("BEL_SON","BEL_CHI"),
                                                             c("CRI_SON","CRI_CHI"),
                                                             c("CUR_SON","CUR_CHI"),
                                                             c("FLA_SON","FLA_CHI"),
                                                             c("FUS_SON","FUS_CHI"))))
      g = plot.custom.dabest(x=bav,palette="Paired",group.summaries=NULL,
                             rawplot.ylabel="RES_KIPPSINDEX_CALC",
                             #theme=theme_transparentwhite(),
                             effsize.markersize=2)
      #print(g)
      ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/DABEST/DABEST_RES_KIPPSINDEX_CALC.png",
             g,bg="transparent",width=23,height=7,units="in")
    }
    print("WING.LENGTH.TO.PRIMARIES"); {
      
      temp = agg
      temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
      temp = temp[,c("GROUP","WING.LENGTH.TO.PRIMARIES")]
      temp = temp[complete.cases(temp), ]
      bav = dabest(temp,GROUP,WING.LENGTH.TO.PRIMARIES, idx=(list(c("BIL_SON","BIL_CHI"),
                                                                  c("BRU_SON","BRU_CHI"),
                                                                  c("MEL_SON","MEL_CHI"),
                                                                  c("NIT_SON","NIT_CHI"),
                                                                  c("SIN_SON","SIN_CHI"),c("BEL_SON","BEL_CHI"),
                                                                  c("CRI_SON","CRI_CHI"),
                                                                  c("CUR_SON","CUR_CHI"),
                                                                  c("FLA_SON","FLA_CHI"),
                                                                  c("FUS_SON","FUS_CHI"))))
      g = plot.custom.dabest(x=bav,palette=c("red","blue"),group.summaries=NULL,
                             rawplot.ylabel="WING.LENGTH.TO.PRIMARIES",
                             #theme=theme_transparentwhite(),
                             effsize.markersize=2)
      print(g)
      ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/DABEST/DABEST_WING.LENGTH.TO.PRIMARIES.png",
             g,bg="transparent",width=23,height=7,units="in")
      
      temp = agg_res
      temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
      temp = temp[,c("GROUP","RES_WING.LENGTH.TO.PRIMARIES")]
      temp = temp[complete.cases(temp), ]
      bav = dabest(temp,GROUP,RES_WING.LENGTH.TO.PRIMARIES, idx=(list(c("BIL_SON","BIL_CHI"),
                                                                      c("BRU_SON","BRU_CHI"),
                                                                      c("MEL_SON","MEL_CHI"),
                                                                      c("NIT_SON","NIT_CHI"),
                                                                      c("SIN_SON","SIN_CHI"),c("BEL_SON","BEL_CHI"),
                                                                      c("CRI_SON","CRI_CHI"),
                                                                      c("CUR_SON","CUR_CHI"),
                                                                      c("FLA_SON","FLA_CHI"),
                                                                      c("FUS_SON","FUS_CHI"))))
      g = plot.custom.dabest(x=bav,palette="Paired",group.summaries=NULL,
                             rawplot.ylabel="RES_WING.LENGTH.TO.PRIMARIES",
                             #theme=theme_transparentwhite(),
                             effsize.markersize=2)
      #print(g)
      ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/DABEST/DABEST_RES_WING.LENGTH.TO.PRIMARIES.png",
             g,bg="transparent",width=23,height=7,units="in")
    }
    print("LONG"); {
      temp = agg
      temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
      temp = temp[,c("GROUP","LONG")]
      temp = temp[complete.cases(temp), ]
      bav = dabest(temp,GROUP,LONG, idx=(list(c("BIL_SON","BIL_CHI"),
                                              c("BRU_SON","BRU_CHI"),
                                              c("MEL_SON","MEL_CHI"),
                                              c("NIT_SON","NIT_CHI"),
                                              c("SIN_SON","SIN_CHI"),c("BEL_SON","BEL_CHI"),
                                              c("CRI_SON","CRI_CHI"),
                                              c("CUR_SON","CUR_CHI"),
                                              c("FLA_SON","FLA_CHI"),
                                              c("FUS_SON","FUS_CHI"))))
      g = plot.custom.dabest(x=bav,palette=c("red","blue"),group.summaries=NULL,
                             rawplot.ylabel="LONG",
                             #theme=theme_transparentwhite(),
                             effsize.markersize=2)
      print(g)
      ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/DABEST/DABEST_LONG.png",
             g,bg="transparent",width=23,height=7,units="in")
      
      
      
  
      
    }
    print("BEAKLATERALSURFACE_CALC"); {
      temp = agg
      temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
      temp = temp[,c("GROUP","BEAKLATERALSURFACE_CALC")]
      temp = temp[complete.cases(temp), ]
      bav = dabest(temp,GROUP,BEAKLATERALSURFACE_CALC, idx=(list(c("BIL_SON","BIL_CHI"),
                                                                 c("BRU_SON","BRU_CHI"),
                                                                 c("MEL_SON","MEL_CHI"),
                                                                 c("NIT_SON","NIT_CHI"),
                                                                 c("SIN_SON","SIN_CHI"),c("BEL_SON","BEL_CHI"),
                                                                 c("CRI_SON","CRI_CHI"),
                                                                 c("CUR_SON","CUR_CHI"),
                                                                 c("FLA_SON","FLA_CHI"),
                                                                 c("FUS_SON","FUS_CHI"))))
      g = plot.custom.dabest(x=bav,palette=c("red","blue"),group.summaries=NULL,
                             rawplot.ylabel="BEAKLATERALSURFACE_CALC",
                             #theme=theme_transparentwhite(),
                             effsize.markersize=2)
      print(g)
      ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/DABEST/DABEST_BEAKLATERALSURFACE_CALC.png",
             g,bg="transparent",width=23,height=7,units="in")
      
      temp = agg_res
      temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
      temp = temp[,c("GROUP","RES_BEAKLATERALSURFACE_CALC")]
      temp = temp[complete.cases(temp), ]
      bav = dabest(temp,GROUP,RES_BEAKLATERALSURFACE_CALC, idx=(list(c("BIL_SON","BIL_CHI"),
                                                                     c("BRU_SON","BRU_CHI"),
                                                                     c("MEL_SON","MEL_CHI"),
                                                                     c("NIT_SON","NIT_CHI"),
                                                                     c("SIN_SON","SIN_CHI"),c("BEL_SON","BEL_CHI"),
                                                                     c("CRI_SON","CRI_CHI"),
                                                                     c("CUR_SON","CUR_CHI"),
                                                                     c("FLA_SON","FLA_CHI"),
                                                                     c("FUS_SON","FUS_CHI"))))
      g = plot.custom.dabest(x=bav,palette="Paired",group.summaries=NULL,
                             rawplot.ylabel="RES_BEAKLATERALSURFACE_CALC",
                             #theme=theme_transparentwhite(),
                             effsize.markersize=2)
      #print(g)
      ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/DABEST/DABEST_RES_BEAKLATERALSURFACE_CALC.png",
             g,bg="transparent",width=23,height=7,units="in")
    }
    print("BILL.HEIGHT"); {
      temp = agg
      temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
      temp = temp[,c("GROUP","BILL.HEIGHT")]
      temp = temp[complete.cases(temp), ]
      bav = dabest(temp,GROUP,BILL.HEIGHT, idx=(list(c("BIL_SON","BIL_CHI"),
                                                     c("BRU_SON","BRU_CHI"),
                                                     c("MEL_SON","MEL_CHI"),
                                                     c("NIT_SON","NIT_CHI"),
                                                     c("SIN_SON","SIN_CHI"),c("BEL_SON","BEL_CHI"),
                                                     c("CRI_SON","CRI_CHI"),
                                                     c("CUR_SON","CUR_CHI"),
                                                     c("FLA_SON","FLA_CHI"),
                                                     c("FUS_SON","FUS_CHI"))))
      g = plot.custom.dabest(x=bav,palette=c("red","blue"),group.summaries=NULL,
                             rawplot.ylabel="BILL.HEIGHT",
                             #theme=theme_transparentwhite(),
                             effsize.markersize=2)
      print(g)
      ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/DABEST/DABEST_BILL.HEIGHT.png",
             g,bg="transparent",width=23,height=7,units="in")
      
      temp = agg_res
      temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
      temp = temp[,c("GROUP","RES_BILL.HEIGHT")]
      temp = temp[complete.cases(temp), ]
      bav = dabest(temp,GROUP,RES_BILL.HEIGHT, idx=(list(c("BIL_SON","BIL_CHI"),
                                                         c("BRU_SON","BRU_CHI"),
                                                         c("MEL_SON","MEL_CHI"),
                                                         c("NIT_SON","NIT_CHI"),
                                                         c("SIN_SON","SIN_CHI"),c("BEL_SON","BEL_CHI"),
                                                         c("CRI_SON","CRI_CHI"),
                                                         c("CUR_SON","CUR_CHI"),
                                                         c("FLA_SON","FLA_CHI"),
                                                         c("FUS_SON","FUS_CHI"))))
      g = plot.custom.dabest(x=bav,palette="Paired",group.summaries=NULL,
                             rawplot.ylabel="RES_BILL.HEIGHT",
                             #theme=theme_transparentwhite(),
                             effsize.markersize=2)
      #print(g)
      ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/DABEST/DABEST_RES_BILL.HEIGHT.png",
             g,bg="transparent",width=23,height=7,units="in")
      
    }
    print("BILL.WIDTH"); {
      temp = agg
      temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
      temp = temp[,c("GROUP","BILL.WIDTH")]
      temp = temp[complete.cases(temp), ]
      bav = dabest(temp,GROUP,BILL.WIDTH, idx=(list(c("BIL_SON","BIL_CHI"),
                                                    c("BRU_SON","BRU_CHI"),
                                                    c("MEL_SON","MEL_CHI"),
                                                    c("NIT_SON","NIT_CHI"),
                                                    c("SIN_SON","SIN_CHI"),c("BEL_SON","BEL_CHI"),
                                                    c("CRI_SON","CRI_CHI"),
                                                    c("CUR_SON","CUR_CHI"),
                                                    c("FLA_SON","FLA_CHI"),
                                                    c("FUS_SON","FUS_CHI"))))
      g = plot.custom.dabest(x=bav,palette=c("red","blue"),group.summaries=NULL,
                             rawplot.ylabel="BILL.WIDTH",
                             #theme=theme_transparentwhite(),
                             effsize.markersize=2)
      print(g)
      ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/DABEST/DABEST_BILL.WIDTH.png",
             g,bg="transparent",width=23,height=7,units="in")
      
      temp = agg_res
      temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
      temp = temp[,c("GROUP","RES_BILL.WIDTH")]
      temp = temp[complete.cases(temp), ]
      bav = dabest(temp,GROUP,RES_BILL.WIDTH, idx=(list(c("BIL_SON","BIL_CHI"),
                                                        c("BRU_SON","BRU_CHI"),
                                                        c("MEL_SON","MEL_CHI"),
                                                        c("NIT_SON","NIT_CHI"),
                                                        c("SIN_SON","SIN_CHI"),c("BEL_SON","BEL_CHI"),
                                                        c("CRI_SON","CRI_CHI"),
                                                        c("CUR_SON","CUR_CHI"),
                                                        c("FLA_SON","FLA_CHI"),
                                                        c("FUS_SON","FUS_CHI"))))
      g = plot.custom.dabest(x=bav,palette="Paired",group.summaries=NULL,
                             rawplot.ylabel="RES_BILL.WIDTH",
                             #theme=theme_transparentwhite(),
                             effsize.markersize=2)
      #print(g)
      ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/DABEST/DABEST_RES_BILL.WIDTH.png",
             g,bg="transparent",width=23,height=7,units="in")
    }
    print("BEAKTOTALSURFACE_CALC"); {
      temp = agg
      temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
      temp = temp[,c("GROUP","BEAKTOTALSURFACE_CALC")]
      temp = temp[complete.cases(temp), ]
      bav = dabest(temp,GROUP,BEAKTOTALSURFACE_CALC, idx=(list(c("BIL_SON","BIL_CHI"),
                                                               c("BRU_SON","BRU_CHI"),
                                                               c("MEL_SON","MEL_CHI"),
                                                               c("NIT_SON","NIT_CHI"),
                                                               c("SIN_SON","SIN_CHI"),c("BEL_SON","BEL_CHI"),
                                                               c("CRI_SON","CRI_CHI"),
                                                               c("CUR_SON","CUR_CHI"),
                                                               c("FLA_SON","FLA_CHI"),
                                                               c("FUS_SON","FUS_CHI"))))
      g = plot.custom.dabest(x=bav,palette=c("red","blue"),group.summaries=NULL,
                             rawplot.ylabel="BEAKTOTALSURFACE_CALC",
                             #theme=theme_transparentwhite(),
                             effsize.markersize=2)
      print(g)
      ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/DABEST/DABEST_BEAKTOTALSURFACE_CALC.png",
             g,bg="transparent",width=23,height=7,units="in")
      
      temp = agg_res
      temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
      temp = temp[,c("GROUP","RES_BEAKTOTALSURFACE_CALC")]
      temp = temp[complete.cases(temp), ]
      bav = dabest(temp,GROUP,RES_BEAKTOTALSURFACE_CALC, idx=(list(c("BIL_SON","BIL_CHI"),
                                                                   c("BRU_SON","BRU_CHI"),
                                                                   c("MEL_SON","MEL_CHI"),
                                                                   c("NIT_SON","NIT_CHI"),
                                                                   c("SIN_SON","SIN_CHI"),c("BEL_SON","BEL_CHI"),
                                                                   c("CRI_SON","CRI_CHI"),
                                                                   c("CUR_SON","CUR_CHI"),
                                                                   c("FLA_SON","FLA_CHI"),
                                                                   c("FUS_SON","FUS_CHI"))))
      g = plot.custom.dabest(x=bav,palette="Paired",group.summaries=NULL,
                             rawplot.ylabel="RES_BEAKTOTALSURFACE_CALC",
                             #theme=theme_transparentwhite(),
                             effsize.markersize=2)
      #print(g)
      ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/DABEST/DABEST_RES_BEAKTOTALSURFACE_CALC.png",
             g,bg="transparent",width=23,height=7,units="in")
      
    }
    print("BEAKBASEAREA_CALC"); {
      temp = agg
      temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
      temp = temp[,c("GROUP","BEAKBASEAREA_CALC")]
      temp = temp[complete.cases(temp), ]
      bav = dabest(temp,GROUP,BEAKBASEAREA_CALC, idx=(list(c("BIL_SON","BIL_CHI"),
                                                           c("BRU_SON","BRU_CHI"),
                                                           c("MEL_SON","MEL_CHI"),
                                                           c("NIT_SON","NIT_CHI"),
                                                           c("SIN_SON","SIN_CHI"),c("BEL_SON","BEL_CHI"),
                                                           c("CRI_SON","CRI_CHI"),
                                                           c("CUR_SON","CUR_CHI"),
                                                           c("FLA_SON","FLA_CHI"),
                                                           c("FUS_SON","FUS_CHI"))))
      g = plot.custom.dabest(x=bav,palette=c("red","blue"),group.summaries=NULL,
                             rawplot.ylabel="BEAKBASEAREA_CALC",
                             #theme=theme_transparentwhite(),
                             effsize.markersize=2)
      print(g)
      ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/DABEST/DABEST_BEAKBASEAREA_CALC.png",
             g,bg="transparent",width=23,height=7,units="in")
      
      temp = agg_res
      temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
      temp = temp[,c("GROUP","RES_BEAKBASEAREA_CALC")]
      temp = temp[complete.cases(temp), ]
      bav = dabest(temp,GROUP,RES_BEAKBASEAREA_CALC, idx=(list(c("BIL_SON","BIL_CHI"),
                                                               c("BRU_SON","BRU_CHI"),
                                                               c("MEL_SON","MEL_CHI"),
                                                               c("NIT_SON","NIT_CHI"),
                                                               c("SIN_SON","SIN_CHI"),c("BEL_SON","BEL_CHI"),
                                                               c("CRI_SON","CRI_CHI"),
                                                               c("CUR_SON","CUR_CHI"),
                                                               c("FLA_SON","FLA_CHI"),
                                                               c("FUS_SON","FUS_CHI"))))
      g = plot.custom.dabest(x=bav,palette="Paired",group.summaries=NULL,
                             rawplot.ylabel="RES_BEAKBASEAREA_CALC",
                             #theme=theme_transparentwhite(),
                             effsize.markersize=2)
      #print(g)
      ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/DABEST/DABEST_RES_BEAKBASEAREA_CALC.png",
             g,bg="transparent",width=23,height=7,units="in")
      
    }
    
    
  }
  
  
  
  
  
  
  #pdf(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/DABEST/DABEST_","RES_TAIL",".pdf",sep=""))
}; dev.off(); names(agg)

thermo_table = rbind(c(2,1),c(3,4))
disp_table = rbind(c(4,4),c(1,1))
colnames(thermo_table) = c("YESSTR","NOSTR")
rownames(thermo_table) = c("YESTH","NOTH")
colnames(disp_table)= c("YESSTR","NOSTR")
rownames(thermo_table) = c("YESTH","NOTH")

fisher.test(thermo_table)

testtable =  rbind(c(0,5),c(5,0))
fisher.test(testtable)


## get magnitude of differences? 
print("DIFS"); {
  relative_means_full = as.data.frame(lapply(unique(agg_good$SPP),FUN=function(spp) {
    
    temp = agg_good[,c(24,37,4:10,41:47)]
    
    bellii = temp[temp$SPP==spp,]
    bellii_1 = bellii[bellii$WHICH.SIDE.OF.CFB==unique(bellii$WHICH.SIDE.OF.CFB)[1],]
    bellii_2 = bellii[bellii$WHICH.SIDE.OF.CFB==unique(bellii$WHICH.SIDE.OF.CFB)[2],]
    
    x = sapply(3:ncol(bellii),FUN=function(colnum) {
      
      #print(colnum)
      meanval = mean(bellii[,colnum],na.rm=T)
      meanval_son = mean(bellii_1[,colnum],na.rm=T)
      meanval_chi = mean(bellii_2[,colnum],na.rm=T)
      diffmean = abs(meanval_son-meanval_chi)
      diff_relative = abs(diffmean)/abs(meanval)
      names(diff_relative) = colnames(bellii)[colnum]
      return(diff_relative)
    })
    
    y = as.data.frame(x)
    colnames(y) = spp
    
    return(y)
    
    
    
  }))
  
  meanmeans = colMeans(relative_means_full,na.rm=T)
  summeans = colSums(relative_means_full,na.rm=T)
  barplot(meanmeans,las=2)
  
  barplot(rowMeans(relative_means_full,na.rm=T),las=2)
  summary(relative_means_full)
  
  
  corrplot(as.matrix(relative_means_full*100),is.corr=F,
           method="color")
  
  options(scipen=999)
  summary(unlist(relative_means_full))
  quans= quantile(unlist(relative_means_full),
                  seq(0,1,0.05),
                  na.rm=T)
  quan095 = quantile(unlist(relative_means_full),
                     0.95,
                     na.rm=T)
  boxplot(relative_means_full*100,las=2,log="y",ylab="% dif in mean (log)")
  abline(h=quans[c(11,16,20)]*100,lty=c(2:4),col=c("darkred","red","pink"))
  
  pdf("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/DABEST/percent_dif_in_mean_log.pdf")
  for (i in 1:length(colnames(relative_means_full))) {
    
    plot(relative_means_full[,i]*100,type="n",main=colnames(relative_means_full)[i],
         ylab="% dif in mean (log)")
    abline(h=quans[c(11,16,20)]*100,lty=c(2:4),col=c("darkred","red","pink"))
    text(y=relative_means_full[,i]*100,x=1:20,labels=rownames(relative_means_full),
         cex=0.5)
    
  }
  dev.off()
}; dev.off()


## work with relative means to check whether there is correlation between that and fst/dxy?
print("MEANS"); {
  x = rbind(
    c(7.89,11,7,1,0.33,0.06,43,400,0.007375754,0.058548283,0.096541807,0.077926649,0.002663856,0.020838178,0.124616436,0.105121429,0.106405876,0.174450778,0.113679375,0.11102855,0.050757185,0.055274759),
    c(2.00,1,0,0,0.3,0.03,555,120,0.005847492,0.0169948,0.003254573,0.01476654,0.035483751,0.026392264,0.018739254,0.056992487,0.007457418,0.027108887,0.02252589,0.015443793,0.013940705,0.015610247),
    c(3.53,2,1,0,0.33,0.04,119,114,0.021283565,0.025217758,0.034841221,0.070203554,0.03209628,0.011225541,0.015374856,0.148220315,0.008261637,0.033961821,0.029390245,0.023444914,0.020203154,0.020532785),
    c(10.04,12,9,2,0.31,0.15,74,93,0.08971665,0.05281363,0.08320033,0.08959619,0.03462309,0.01857874,0.05173181,0.14854596,0.16992847,0.2168172,0.13664233,0.14424397,0.08903684,0.08138562),
    c(2.36,2,1,0,0.33,0.1,19,717,0.009560529,0.026629974,0.014105447,0.008411393,0.002209541,0.017347865,0.0583799,0.137986026,0.026592289,0.0033148,0.012476909,0.001594604,0.011080025,0.001220705),
    c(3.77,7,1,0,0.35,0.06,230,181,0.005020587,0.027997082,0.025506217,0.047541158,0.037720403,0.040473093,0.087173706,0.017428578,0.039831148,0.077268087,0.047334239,0.0439726,0.013179647,0.018727452),
    c(4.27,7,1,0,0.31,0.06,224,811,0.038788728,0.071844795,0.002650374,0.063320214,0.030011037,0.050195267,0.062807151,0.112375496,0.040505802,0.031671225,0.048502929,0.004468169,0.017182628,0.024807062),
    c(2.95,3,0,0,0.32,0.04,238,81,0.015508571,0.034503047,0.022586617,0.056409532,0.006064321,0.013564112,0.074924513,0.057975234,0.003879043,0.036452455,0.032510795,0.021354798,0.016093268,0.021514494),
    c(8.83,10,7,2,0.34,0.04,134,142,0.03897788,0.04793591,0.143350701,0.002090229,0.006999664,0.005692018,0.032640636,0.087173845,0.179589332,0.224925016,0.128329107,0.157544881,0.108007725,0.072948384)
  )
  
  colnames(x)=c("Meanmorphdifference","above50q","above75q","above95q","Meandxy","Meanfst","NuminIslands","NuminSweeps","BILL.HEIGHT","BILL.LENGTH","BILL.WIDTH","TARSUS.LENGTH","WING.LENGTH.TO.PRIMARIES","WING.LENGTH.TO.SECONDARIES","TAIL","KIPPSINDEX_CALC","BEAKBASEAREA_CALC","BEAKVOL_CALC","BEAKLATERALSURFACE_CALC","BEAKTOTALSURFACE_CALC","BEAKAREAVOLUME_CALC","BEAKTOTAREAVOLUME_CALC")
  
  x=as.data.frame(x)
  plot(x)
  cor(x)
  y = cor(x)[,c(5:8)]
  y = y[c(-5:-8),]
  
  corrplot(cor(x),method="ellipse",diag=F,type="upper")
  
  for (i in 5:8) {
    for (j in c(1:4,9:ncol(x))) {
      
      mod = glm(x[,i] ~ x[,j])
      #if (coef(summary(mod))[2,4] < 0.10 ) {
      print(names(x)[i])
      print(names(x)[j])
      print(coef(summary(mod))[2,4])
      
      #print("SIGNIFICANT")
      #}
    }
  }
  
  
  
  mod = glm(x$Meandxy ~ x$Meanmorphdifference +
              x$above50q + x$above75q + x$above95q +
              x$BILL.HEIGHT + x$BILL.LENGTH + x$BILL.WIDTH +
              x$TAIL + x$TARSUS.LENGTH + x$WING.LENGTH.TO.PRIMARIES +
              x$WING.LENGTH.TO.SECONDARIES + x$KIPPSINDEX_CALC + 
              x$BEAKAREAVOLUME_CALC + x$BEAKBASEAREA_CALC +
              x$BEAKLATERALSURFACE_CALC + x$BEAKTOTALSURFACE_CALC +
              x$BEAKTOTAREAVOLUME_CALC + x$BEAKVOL_CALC)
  summary(mod)
}
}

## perform a pca 
print("PCA"); {
  summary(agg)
  
  keepcols=c("BEAKAREAVOLUME_CALC","BEAKBASEAREA_CALC","BEAKLATERALSURFACE_CALC",
             "BEAKTOTALSURFACE_CALC","BEAKTOTAREAVOLUME_CALC","BEAKVOL_CALC",
             "BILL.HEIGHT","BILL.LENGTH","BILL.WIDTH","KIPPSINDEX_CALC","TAIL",
             "TARSUS.LENGTH","WING.LENGTH.TO.PRIMARIES","WING.LENGTH.TO.SECONDARIES")
  
  metacols=c("AGE","BAD.MEASUREMENT","CATALOG.NUMBER","CONDITION","COUNTRY","COUNTY",
             "DATE","GENETIC.SIDE","GENETIC.SPLIT","GENUS","GEOREF.BY","LAT","LOCALITY",
             "LOCALITY2","LONG","MEASUREMENT","MEASURER","NOTES","NUM2","OLDCAT","SEQUENCED.",
             "SEX","SPECIES","SPP","STATE","STRUCTURE","SUBSPP","WHICH.SIDE.OF.CFB")
  
  forpca_full = (agg[,c(keepcols,metacols)])
  forpca_notail = forpca_full[,-(which(colnames(forpca_full)=="TAIL"))]
  forpca_notailsec = forpca_notail[,-(which(colnames(forpca_notail)=="WING.LENGTH.TO.SECONDARIES"))]
  
  forpca_full = forpca_full[complete.cases(forpca_full[,keepcols]),]
  forpca_notail = forpca_notail[complete.cases(forpca_notail[,keepcols[keepcols!="TAIL"]]),]
  forpca_notailsec = forpca_notailsec[complete.cases(forpca_notailsec[,keepcols!="TAIL" && keepcols!="WING.LENGTH.TO.SECONDARIES"]),]
  
  pca_full = prcomp(forpca_full[,-(which(colnames(forpca_full) %in% metacols))],scale=T,center=T)
  pca_notail = prcomp(forpca_notail[,-(which(colnames(forpca_notail) %in% metacols))],scale=T,center=T)
  pca_notailsec = prcomp(forpca_notailsec[,-(which(colnames(forpca_notailsec) %in% metacols))],scale=T,center=T)
  
  summary(pca_full)
  summary(pca_notail)
  summary(pca_notailsec)
  
  axes_full = cbind(pca_full$x,forpca_full)
  axes_notail = cbind(pca_notail$x,forpca_notail)
  axes_notailsec = cbind(pca_notailsec$x,forpca_notailsec)

  write.csv(axes_full,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/measurements/aggregated datasets/KLP_Master_Spreadsheet_Morphology_AGGREGATED_10_April_2020_PCASCORES_FULL.csv",
            row.names = F)
  write.csv(axes_notail,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/measurements/aggregated datasets/KLP_Master_Spreadsheet_Morphology_AGGREGATED_10_April_2020_PCASCORES_NOTAIL.csv",
            row.names = F)
  write.csv(axes_notailsec,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/measurements/aggregated datasets/KLP_Master_Spreadsheet_Morphology_AGGREGATED_10_April_2020_PCASCORES_NOTAILSEC.csv",
            row.names = F)
  
  
  ## DO DABESTR
  temp = axes_full
  temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
  color=as.numeric(as.factor(temp$WHICH.SIDE.OF.CFB))
  temp$COLOR=color
  temp = temp[,c("GROUP","PC1","WHICH.SIDE.OF.CFB")]
  color=color[complete.cases(temp)]
  temp = temp[complete.cases(temp), ]
  
  
    #pdf(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/DABEST_","RES_PC1",".pdf",sep=""))
  #par(mfrow=c(2,5))
  bav = dabestr::dabest(temp,GROUP,PC1, idx=(list(c("BIL_SON","BIL_CHI"),
                                                  c("BRU_SON","BRU_CHI"),
                                                  c("MEL_SON","MEL_CHI"),
                                                  c("NIT_SON","NIT_CHI"),
                                                  c("SIN_SON","SIN_CHI"),
                                                  c("BEL_SON","BEL_CHI"),
                                                  c("CRI_SON","CRI_CHI"),
                                                  c("CUR_SON","CUR_CHI"),
                                                  c("FLA_SON","FLA_CHI"),
                                                  c("FUS_SON","FUS_CHI"))))
  # g = plot.custom.dabest(x=bav,palette="Paired",group.summaries=NULL,
  #                        rawplot.ylabel="PC1",
  #                        #theme=theme_transparentwhite(),
  #                        effsize.markersize=2)
  g = plot(x=bav,palette="Paired",group.summaries=NULL,
                         rawplot.ylabel="PC1",color.column="GROUP",
                         #theme=theme_transparentwhite(),
                         effsize.markersize=2)
  #print(g)
  ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/PCA/DABEST_PCA_PC1.png",
         g,bg="transparent",width=19,height=5,units="in")
  
  temp = axes_full
  temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
  temp = temp[,c("GROUP","PC2","WHICH.SIDE.OF.CFB")]
  temp = temp[complete.cases(temp), ]
  #pdf(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/DABEST_","RES_PC2",".pdf",sep=""))
  #par(mfrow=c(2,5))
  bav = dabestr::dabest(temp,GROUP,PC2, idx=(list(c("BIL_SON","BIL_CHI"),
                                                  c("BRU_SON","BRU_CHI"),
                                                  c("MEL_SON","MEL_CHI"),
                                                  c("NIT_SON","NIT_CHI"),
                                                  c("SIN_SON","SIN_CHI"),
                                                  c("BEL_SON","BEL_CHI"),
                                                  c("CRI_SON","CRI_CHI"),
                                                  c("CUR_SON","CUR_CHI"),
                                                  c("FLA_SON","FLA_CHI"),
                                                  c("FUS_SON","FUS_CHI"))))
  g = plot(x=bav,palette="Paired",group.summaries=NULL,
                         rawplot.ylabel="PC2",color.column = "GROUP",
                         #theme=theme_transparentwhite(),
                         effsize.markersize=2)
  #print(g)
  ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/PCA/DABEST_PCA_PC2.png",
         g,bg="transparent",width=19,height=5,units="in")
  
  
  
  
  temp = axes_full
  temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
  temp = temp[,c("GROUP","PC3","WHICH.SIDE.OF.CFB")]
  temp = temp[complete.cases(temp), ]
  #pdf(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/DABEST_","RES_PC3",".pdf",sep=""))
  #par(mfrow=c(2,5))
  bav = dabestr::dabest(temp,GROUP,PC3, idx=(list(c("BIL_SON","BIL_CHI"),
                                                  c("BRU_SON","BRU_CHI"),
                                                  c("MEL_SON","MEL_CHI"),
                                                  c("NIT_SON","NIT_CHI"),
                                                  c("SIN_SON","SIN_CHI"),
                                                  c("BEL_SON","BEL_CHI"),
                                                  c("CRI_SON","CRI_CHI"),
                                                  c("CUR_SON","CUR_CHI"),
                                                  c("FLA_SON","FLA_CHI"),
                                                  c("FUS_SON","FUS_CHI"))))
  g = plot(x=bav,palette="Paired",group.summaries=NULL,
                         rawplot.ylabel="PC3",color.column = "GROUP",
                         #theme=theme_transparentwhite(),
                         effsize.markersize=2)
  #print(g)
  ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/PCA/DABEST_PCA_PC3.png",
         g,bg="transparent",width=19,height=5,units="in")
  
  pc1quan=quantile(axes_full$PC1,c(0,1))
  pc2quan=quantile(axes_full$PC2,c(0,1))
  pc3quan=quantile(axes_full$PC3,c(0,1))
  pc1diflim=c(pc1quan[1]-pc1quan[2],pc1quan[2]-pc1quan[1])
  pc2diflim=c(pc2quan[1]-pc2quan[2],pc2quan[2]-pc2quan[1])
  pc3diflim=c(pc3quan[1]-pc3quan[2],pc3quan[2]-pc3quan[1])
  
  for (species_index in 1:length(unique(axes_full$SPP))){
    spp=unique(axes_full$SPP)[species_index]
    print(as.character(spp))
    
    
    temp = axes_full[axes_full$SPP==spp,]
    temp$WHICH.SIDE.OF.CFB = substr(temp$WHICH.SIDE.OF.CFB,1,3)
    temp = temp[,c("WHICH.SIDE.OF.CFB","PC1","PC2","PC3")]
    temp = temp[complete.cases(temp), ]
    
    
    dif1=quantile(temp$PC1,1)-quantile(temp$PC1,0)
    dif1=unique(as.numeric(c(-1*abs(dif1),abs(dif1))))
    
    #pdf(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/DABEST_","RES_PC3",".pdf",sep=""))
    #par(mfrow=c(2,5))
    
    quan1=c(quantile(temp$PC1,0)-1,quantile(temp$PC1,1)+1)
    bav1 = dabestr::dabest(temp,WHICH.SIDE.OF.CFB,PC1, 
                           idx=(list(unique(temp$WHICH.SIDE.OF.CFB))))
    g1 = plot(x=bav1,palette="Paired",group.summaries=NULL,
              rawplot.ylabel="PC1",rawplot.ylim=quan1,
              # rawplot.names=c("Chihuahuan","Sonoran"),
              effsize.ylim=dif1,
              float.contrast=F,
              effsize.markersize=2)
    #print(g)
    ggsave(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/PCA/DABEST_PCA_PC1_",spp,".png",sep=""),
           g1,bg="transparent",width=4,height=4,units="in")
    
    dif2=quantile(temp$PC2,1)-quantile(temp$PC2,0)
    dif2=unique(as.numeric(c(-1*abs(dif2),abs(dif2))))
    quan2=c(quantile(temp$PC2,0)-1,quantile(temp$PC2,1)+1)
    bav2 = dabestr::dabest(temp,WHICH.SIDE.OF.CFB,PC2, idx=(list(unique(temp$WHICH.SIDE.OF.CFB))))
    g2 = plot(x=bav2,palette="Paired",group.summaries=NULL,
              rawplot.ylabel="PC2",rawplot.ylim=quan2,
              # rawplot.names=c("Chihuahuan","Sonoran"),
              effsize.ylim=dif2,
              float.contrast=F,
              effsize.markersize=2)
    #print(g)
    ggsave(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/PCA/DABEST_PCA_PC2_",spp,".png",sep=""),
           g2,bg="transparent",width=4,height=4,units="in")
    
    dif3=quantile(temp$PC3,1)-quantile(temp$PC3,0)
    dif3=unique(as.numeric(c(-1*abs(dif3),abs(dif3))))
    quan3=c(quantile(temp$PC3,0)-1,quantile(temp$PC3,1)+1)
    bav3 = dabestr::dabest(temp,WHICH.SIDE.OF.CFB,PC3, idx=(list(unique(temp$WHICH.SIDE.OF.CFB))))
    g3 = plot(x=bav3,palette="Paired",group.summaries=NULL,
             rawplot.ylabel="PC3",rawplot.ylim=quan3,
            # rawplot.names=c("Chihuahuan","Sonoran"),
             effsize.ylim=dif3,
             float.contrast=F,
             effsize.markersize=2)
    #print(g)
    ggsave(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/PCA/DABEST_PCA_PC3_",spp,".png",sep=""),
           g3,bg="transparent",width=4,height=4,units="in")
    
  }
  
  
  
  ## WHITE
  {
  {
    palette(c("white","red","orange","goldenrod","green","blue",
              "cyan","magenta","brown","pink","lightblue"))
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white"
    )
    g <- ggbiplot(pca_full, obs.scale = 1, var.scale = 1, choices = c(1,2),
                  groups = axes_full$SPP, ellipse = F, 
                  circle = F,var.axes = T, alpha=0,varname.size=1,
                  arrow.color="white")
    
    g = g + ggtitle("")
    g = g + theme(plot.title = element_text(color="white", face="italic",hjust=0))
    g <- g + theme(legend.direction = 'horizontal', 
                   legend.position = 'bottom',
                   legend.background = element_rect(fill = "transparent", colour = NA),
                   legend.text = element_text(colour = 'white'),
                   legend.key = element_rect(fill="transparent",colour=NA))
    g = g + theme(panel.background = element_rect(fill = "transparent", colour = NA))
    g = g + theme(plot.background = element_rect(fill = "transparent", colour = NA),
                  panel.grid = element_blank())
    g = g + theme(axis.text.x=element_text(colour="white"),
                  axis.text.y=element_text(colour="white"))
    g = g + theme(plot.title = element_text(color="white"),
                  axis.title.x = element_text(color="white"),
                  axis.title.y = element_text(color="white"))
    g = g + theme(line = element_line(color="white"),
                  rect= element_rect(color="white"),
                  text=element_text(color="white"))
    g = g + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
    print(g)
  }
  ggsave("Analyses and Images/PCA/pca_axis_loadings_P1P2.png", g, bg = "transparent")
  {
    palette(c("white","red","orange","goldenrod","green","blue",
              "cyan","magenta","brown","pink","lightblue"))
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white"
    )
    g <- ggbiplot(pca_full, obs.scale = 1, var.scale = 1, choices = c(1,3),
                  groups = axes_full$SPP, ellipse = F, 
                  circle = F,var.axes = T, alpha=0,varname.size=1,
                  arrow.color="white")
    
    g = g + ggtitle("")
    g = g + theme(plot.title = element_text(color="white", face="italic",hjust=0))
    g <- g + theme(legend.direction = 'horizontal', 
                   legend.position = 'bottom',
                   legend.background = element_rect(fill = "transparent", colour = NA),
                   legend.text = element_text(colour = 'white'),
                   legend.key = element_rect(fill="transparent",colour=NA))
    g = g + theme(panel.background = element_rect(fill = "transparent", colour = NA))
    g = g + theme(plot.background = element_rect(fill = "transparent", colour = NA),
                  panel.grid = element_blank())
    g = g + theme(axis.text.x=element_text(colour="white"),
                  axis.text.y=element_text(colour="white"))
    g = g + theme(plot.title = element_text(color="white"),
                  axis.title.x = element_text(color="white"),
                  axis.title.y = element_text(color="white"))
    g = g + theme(line = element_line(color="white"),
                  rect= element_rect(color="white"),
                  text=element_text(color="white"))
    g = g + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
    #print(g)
  }
  ggsave("Analyses and Images/PCA/pca_axis_loadings_P1P3.png", g, bg = "transparent")
  {
    palette(c("white","red","orange","goldenrod","green","blue",
              "cyan","magenta","brown","pink","lightblue"))
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white"
    )
    g <- ggbiplot(pca_full, obs.scale = 1, var.scale = 1, choices = c(2,3),
                  groups = axes_full$SPP, ellipse = F, 
                  circle = F,var.axes = T, alpha=0,varname.size=1,
                  arrow.color="white")
    
    g = g + ggtitle("")
    g = g + theme(plot.title = element_text(color="white", face="italic",hjust=0))
    g <- g + theme(legend.direction = 'horizontal', 
                   legend.position = 'bottom',
                   legend.background = element_rect(fill = "transparent", colour = NA),
                   legend.text = element_text(colour = 'white'),
                   legend.key = element_rect(fill="transparent",colour=NA))
    g = g + theme(panel.background = element_rect(fill = "transparent", colour = NA))
    g = g + theme(plot.background = element_rect(fill = "transparent", colour = NA),
                  panel.grid = element_blank())
    g = g + theme(axis.text.x=element_text(colour="white"),
                  axis.text.y=element_text(colour="white"))
    g = g + theme(plot.title = element_text(color="white"),
                  axis.title.x = element_text(color="white"),
                  axis.title.y = element_text(color="white"))
    g = g + theme(line = element_line(color="white"),
                  rect= element_rect(color="white"),
                  text=element_text(color="white"))
    g = g + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
    #print(g)
  }
  ggsave("Analyses and Images/PCA/pca_axis_loadings_P2P3.png", g, bg = "transparent")
  {
    palette(c("white","red","orange","goldenrod","green","blue",
              "cyan","magenta","brown","pink","lightblue"))
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white"
    )
    g <- ggbiplot(pca_full, choices = c(1,2), obs.scale = 1, var.scale = 1, 
                  groups = axes_full$SPP, ellipse = T, 
                  circle = F,var.axes = F, alpha=0)
    
    g = g +
      scale_color_manual(name="", values=c("red","orange","goldenrod","green","cyan",
                                           "blue","purple","magenta","brown","lightblue",
                                           "pink")) +
      #scale_color_manual(name="", values=c("cyan", "white","red","cyan","cyan",
      #                                     "red","cyan","red","white","white","white")) +
      #scale_color_manual(name="", values=c("cyan", "red","white","cyan","white",
      #                                     "cyan","red","cyan","cyan","red","white")) +
      scale_shape_manual(name="", values=c(1,2,3,4,5,6,7,8,9,10,11)) +
      geom_point(aes(colour=axes_full$SPP, shape=axes_full$SPP))
    g = g + ggtitle("")
    g = g + theme(plot.title = element_text(color="white", face="italic",hjust=0))
    g <- g + theme(legend.direction = 'horizontal', 
                   legend.position = 'bottom',
                   legend.background = element_rect(fill = "transparent", colour = 'white'),
                   legend.text = element_text(colour = 'white'),
                   legend.key = element_rect(fill="transparent",colour="transparent"))
    g = g + theme(panel.background = element_rect(fill = "transparent", colour = 'grey'))
    g = g + theme(plot.background = element_rect(fill = "transparent", colour = 'grey'),
                  panel.grid = element_blank())
    g = g + theme(axis.text.x=element_text(colour="white"),
                  axis.text.y=element_text(colour="white"))
    g = g + theme(plot.title = element_text(color="white"),
                  axis.title.x = element_text(color="white"),
                  axis.title.y = element_text(color="white"))
    g = g + theme(line = element_line(color="white"),
                  rect= element_rect(color="white"),
                  text=element_text(color="white"))
    g = g + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
    
    #print(g)
  }
  ggsave("Analyses and Images/PCA/morphology_pca_PC1PC2.png", g, bg = "transparent")
  {
    palette(c("white","red","orange","goldenrod","green","blue",
              "cyan","magenta","brown","pink","lightblue"))
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white"
    )
    g <- ggbiplot(pca_full, choices = c(1,3), obs.scale = 1, var.scale = 1, 
                  groups = axes_full$SPP, ellipse = T, 
                  circle = F,var.axes = F, alpha=0)
    
    g = g +
      scale_color_manual(name="", values=c("red","orange","goldenrod","green","cyan",
                                           "blue","purple","magenta","brown","lightblue",
                                           "pink")) +
      #scale_color_manual(name="", values=c("cyan", "white","red","cyan","cyan",
      #                                     "red","cyan","red","white","white","white")) +
      #scale_color_manual(name="", values=c("cyan", "red","white","cyan","white",
      #                                     "cyan","red","cyan","cyan","red","white")) +
      scale_shape_manual(name="", values=c(1,2,3,4,5,6,7,8,9,10,11)) +
      geom_point(aes(colour=axes_full$SPP, shape=axes_full$SPP))
    g = g + ggtitle("")
    g = g + theme(plot.title = element_text(color="white", face="italic",hjust=0))
    g <- g + theme(legend.direction = 'horizontal', 
                   legend.position = 'bottom',
                   legend.background = element_rect(fill = "transparent", colour = 'white'),
                   legend.text = element_text(colour = 'white'),
                   legend.key = element_rect(fill="transparent",colour="transparent"))
    g = g + theme(panel.background = element_rect(fill = "transparent", colour = 'grey'))
    g = g + theme(plot.background = element_rect(fill = "transparent", colour = 'grey'),
                  panel.grid = element_blank())
    g = g + theme(axis.text.x=element_text(colour="white"),
                  axis.text.y=element_text(colour="white"))
    g = g + theme(plot.title = element_text(color="white"),
                  axis.title.x = element_text(color="white"),
                  axis.title.y = element_text(color="white"))
    g = g + theme(line = element_line(color="white"),
                  rect= element_rect(color="white"),
                  text=element_text(color="white"))
    g = g + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
    
    #print(g)
  }
  ggsave("Analyses and Images/PCA/morphology_pca_PC1PC3.png", g, bg = "transparent")
  {
    palette(c("white","red","orange","goldenrod","green","blue",
              "cyan","magenta","brown","pink","lightblue"))
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white"
    )
    g <- ggbiplot(pca_full, choices = c(2,3), obs.scale = 1, var.scale = 1, 
                  groups = axes_full$SPP, ellipse = T, 
                  circle = F,var.axes = F, alpha=0)
    
    g = g +
      scale_color_manual(name="", values=c("red","orange","goldenrod","green","cyan",
                                           "blue","purple","magenta","brown","lightblue",
                                           "pink")) +
      #scale_color_manual(name="", values=c("cyan", "white","red","cyan","cyan",
      #                                     "red","cyan","red","white","white","white")) +
      #scale_color_manual(name="", values=c("cyan", "red","white","cyan","white",
      #                                     "cyan","red","cyan","cyan","red","white")) +
      scale_shape_manual(name="", values=c(1,2,3,4,5,6,7,8,9,10,11)) +
      geom_point(aes(colour=axes_full$SPP, shape=axes_full$SPP))
    g = g + ggtitle("")
    g = g + theme(plot.title = element_text(color="white", face="italic",hjust=0))
    g <- g + theme(legend.direction = 'horizontal', 
                   legend.position = 'bottom',
                   legend.background = element_rect(fill = "transparent", colour = 'white'),
                   legend.text = element_text(colour = 'white'),
                   legend.key = element_rect(fill="transparent",colour="transparent"))
    g = g + theme(panel.background = element_rect(fill = "transparent", colour = 'grey'))
    g = g + theme(plot.background = element_rect(fill = "transparent", colour = 'grey'),
                  panel.grid = element_blank())
    g = g + theme(axis.text.x=element_text(colour="white"),
                  axis.text.y=element_text(colour="white"))
    g = g + theme(plot.title = element_text(color="white"),
                  axis.title.x = element_text(color="white"),
                  axis.title.y = element_text(color="white"))
    g = g + theme(line = element_line(color="white"),
                  rect= element_rect(color="white"),
                  text=element_text(color="white"))
    g = g + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
    
    #print(g)
  }
  ggsave("Analyses and Images/PCA/morphology_pca_PC2PC3.png", g, bg = "transparent")
  }
  
  ## BLACK
  {
    {
      palette(c("black","red","orange","goldenrod","green","blue",
                "cyan","magenta","brown","pink","lightblue"))
      par(
        bg = NA,
        col.axis = "black",
        fg = "black",
        col.lab = "black",
        col.main = "black"
      )
      g <- ggbiplot(pca_full, obs.scale = 1, var.scale = 1, choices = c(1,2),
                    groups = axes_full$SPP, ellipse = F, 
                    circle = F,var.axes = T, alpha=0,varname.size=1,
                    arrow.color="black")
      
      g = g + ggtitle("")
      g = g + theme(plot.title = element_text(color="black", face="italic",hjust=0))
      g <- g + theme(legend.direction = 'horizontal', 
                     legend.position = 'bottom',
                     legend.background = element_rect(fill = "white", colour = NA),
                     legend.text = element_text(colour = 'white'),
                     legend.key = element_rect(fill="white",colour=NA))
      g = g + theme(panel.background = element_rect(fill = "lightgrey", colour = NA))
      g = g + theme(plot.background = element_rect(fill = "grey", colour = NA),
                    panel.grid = element_blank())
      g = g + theme(axis.text.x=element_text(colour="black"),
                    axis.text.y=element_text(colour="black"))
      g = g + theme(plot.title = element_text(color="black"),
                    axis.title.x = element_text(color="black"),
                    axis.title.y = element_text(color="black"))
      g = g + theme(line = element_line(color="black"),
                    rect= element_rect(color="black"),
                    text=element_text(color="black"))
      g = g + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
      #print(g)
    }
    ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/PCA/pca_axis_loadings_P1P2_black.png", 
           g, bg = "transparent",width=7,height=5,units = "in")
    {
      palette(c("black","red","orange","goldenrod","green","blue",
                "cyan","magenta","brown","pink","lightblue"))
      par(
        bg = NA,
        col.axis = "black",
        fg = "black",
        col.lab = "black",
        col.main = "black"
      )
      g <- ggbiplot(pca_full, obs.scale = 1, var.scale = 1, choices = c(1,3),
                    groups = axes_full$SPP, ellipse = F, 
                    circle = F,var.axes = T, alpha=0,varname.size=1,
                    arrow.color="black")
      
      g = g + ggtitle("")
      g = g + theme(plot.title = element_text(color="black", face="italic",hjust=0))
      g <- g + theme(legend.direction = 'horizontal', 
                     legend.position = 'bottom',
                     legend.background = element_rect(fill = "white", colour = NA),
                     legend.text = element_text(colour = 'white'),
                     legend.key = element_rect(fill="white",colour=NA))
      g = g + theme(panel.background = element_rect(fill = "lightgrey", colour = NA))
      g = g + theme(plot.background = element_rect(fill = "grey", colour = NA),
                    panel.grid = element_blank())
      g = g + theme(axis.text.x=element_text(colour="black"),
                    axis.text.y=element_text(colour="black"))
      g = g + theme(plot.title = element_text(color="black"),
                    axis.title.x = element_text(color="black"),
                    axis.title.y = element_text(color="black"))
      g = g + theme(line = element_line(color="black"),
                    rect= element_rect(color="black"),
                    text=element_text(color="black"))
      g = g + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
      #print(g)
    }
    ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/PCA/pca_axis_loadings_P1P3_black.png", 
           g, bg = "transparent",width=7,height=5,units = "in")
    {
      palette(c("black","red","orange","goldenrod","green","blue",
                "cyan","magenta","brown","pink","lightblue"))
      par(
        bg = NA,
        col.axis = "black",
        fg = "black",
        col.lab = "black",
        col.main = "black"
      )
      g <- ggbiplot(pca_full, obs.scale = 1, var.scale = 1, choices = c(2,3),
                    groups = axes_full$SPP, ellipse = F, 
                    circle = F,var.axes = T, alpha=0,varname.size=1,
                    arrow.color="black")
      
      g = g + ggtitle("")
      g = g + theme(plot.title = element_text(color="black", face="italic",hjust=0))
      g <- g + theme(legend.direction = 'horizontal', 
                     legend.position = 'bottom',
                     legend.background = element_rect(fill = "white", colour = NA),
                     legend.text = element_text(colour = 'white'),
                     legend.key = element_rect(fill="white",colour=NA))
      g = g + theme(panel.background = element_rect(fill = "lightgrey", colour = NA))
      g = g + theme(plot.background = element_rect(fill = "grey", colour = NA),
                    panel.grid = element_blank())
      g = g + theme(axis.text.x=element_text(colour="black"),
                    axis.text.y=element_text(colour="black"))
      g = g + theme(plot.title = element_text(color="black"),
                    axis.title.x = element_text(color="black"),
                    axis.title.y = element_text(color="black"))
      g = g + theme(line = element_line(color="black"),
                    rect= element_rect(color="black"),
                    text=element_text(color="black"))
      g = g + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
      print(g)
    }
    ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/PCA/pca_axis_loadings_P2P3_black.png", 
           g, bg = "transparent",width=7,height=5,units = "in")
    {
      palette(c("black","red","orange","goldenrod","green","blue",
                "cyan","magenta","brown","pink","lightblue"))
      par(
        bg = NA,
        col.axis = "black",
        fg = "black",
        col.lab = "black",
        col.main = "black"
      )
      g <- ggbiplot(pca_full, choices = c(1,2), obs.scale = 1, var.scale = 1, 
                    groups = axes_full$SPP, ellipse = T, 
                    circle = F,var.axes = F, alpha=0)
      
      g = g +
        scale_color_manual(name="", values=c("red","orange","goldenrod","green","cyan",
                                             "blue","purple","magenta","brown","lightblue",
                                             "pink")) +
        #scale_color_manual(name="", values=c("cyan", "black","red","cyan","cyan",
        #                                     "red","cyan","red","black","black","black")) +
        #scale_color_manual(name="", values=c("cyan", "red","black","cyan","black",
        #                                     "cyan","red","cyan","cyan","red","black")) +
        scale_shape_manual(name="", values=c(1,2,3,4,5,6,7,8,9,10,11)) +
        geom_point(aes(colour=axes_full$SPP, shape=axes_full$SPP))
      g = g + ggtitle("")
      g = g + theme(plot.title = element_text(color="black", face="italic",hjust=0))
      g <- g + theme(legend.direction = 'horizontal', 
                     legend.position = 'bottom',
                     legend.background = element_rect(fill = "white", colour = 'black'),
                     legend.text = element_text(colour = 'black'),
                     legend.key = element_rect(fill="white",colour="white"))
      g = g + theme(panel.background = element_rect(fill = "white", colour = 'grey'))
      g = g + theme(plot.background = element_rect(fill = "white", colour = 'grey'),
                    panel.grid = element_blank())
      g = g + theme(axis.text.x=element_text(colour="black"),
                    axis.text.y=element_text(colour="black"))
      g = g + theme(plot.title = element_text(color="black"),
                    axis.title.x = element_text(color="black"),
                    axis.title.y = element_text(color="black"))
      g = g + theme(line = element_line(color="black"),
                    rect= element_rect(color="black"),
                    text=element_text(color="black"))
      g = g + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
      
      print(g)
    }
    ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/PCA/morphology_pca_PC1PC2_black.png", 
           g, bg = "transparent",width=7,height=5,units = "in")
    {
      palette(c("black","red","orange","goldenrod","green","blue",
                "cyan","magenta","brown","pink","lightblue"))
      par(
        bg = NA,
        col.axis = "black",
        fg = "black",
        col.lab = "black",
        col.main = "black"
      )
      g <- ggbiplot(pca_full, choices = c(1,3), obs.scale = 1, var.scale = 1, 
                    groups = axes_full$SPP, ellipse = T, 
                    circle = F,var.axes = F, alpha=0)
      
      g = g +
        scale_color_manual(name="", values=c("red","orange","goldenrod","green","cyan",
                                             "blue","purple","magenta","brown","lightblue",
                                             "pink")) +
        #scale_color_manual(name="", values=c("cyan", "black","red","cyan","cyan",
        #                                     "red","cyan","red","black","black","black")) +
        #scale_color_manual(name="", values=c("cyan", "red","black","cyan","black",
        #                                     "cyan","red","cyan","cyan","red","black")) +
        scale_shape_manual(name="", values=c(1,2,3,4,5,6,7,8,9,10,11)) +
        geom_point(aes(colour=axes_full$SPP, shape=axes_full$SPP))
      g = g + ggtitle("")
      g = g + theme(plot.title = element_text(color="black", face="italic",hjust=0))
      g <- g + theme(legend.direction = 'horizontal', 
                     legend.position = 'bottom',
                     legend.background = element_rect(fill = "white", colour = 'black'),
                     legend.text = element_text(colour = 'black'),
                     legend.key = element_rect(fill="white",colour="white"))
      g = g + theme(panel.background = element_rect(fill = "white", colour = 'grey'))
      g = g + theme(plot.background = element_rect(fill = "white", colour = 'grey'),
                    panel.grid = element_blank())
      g = g + theme(axis.text.x=element_text(colour="black"),
                    axis.text.y=element_text(colour="black"))
      g = g + theme(plot.title = element_text(color="black"),
                    axis.title.x = element_text(color="black"),
                    axis.title.y = element_text(color="black"))
      g = g + theme(line = element_line(color="black"),
                    rect= element_rect(color="black"),
                    text=element_text(color="black"))
      g = g + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
      
      print(g)
      
      #print(g)
    }
    ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/PCA/morphology_pca_PC1PC3_black.png", 
           g, bg = "transparent",width=7,height=5,units = "in")
    {
      palette(c("black","red","orange","goldenrod","green","blue",
                "cyan","magenta","brown","pink","lightblue"))
      par(
        bg = NA,
        col.axis = "black",
        fg = "black",
        col.lab = "black",
        col.main = "black"
      )
      g <- ggbiplot(pca_full, choices = c(2,3), obs.scale = 1, var.scale = 1, 
                    groups = axes_full$SPP, ellipse = T, 
                    circle = F,var.axes = F, alpha=0)
      
      g = g +
        scale_color_manual(name="", values=c("red","orange","goldenrod","green","cyan",
                                             "blue","purple","magenta","brown","lightblue",
                                             "pink")) +
        #scale_color_manual(name="", values=c("cyan", "black","red","cyan","cyan",
        #                                     "red","cyan","red","black","black","black")) +
        #scale_color_manual(name="", values=c("cyan", "red","black","cyan","black",
        #                                     "cyan","red","cyan","cyan","red","black")) +
        scale_shape_manual(name="", values=c(1,2,3,4,5,6,7,8,9,10,11)) +
        geom_point(aes(colour=axes_full$SPP, shape=axes_full$SPP))
      g = g + ggtitle("")
      g = g + theme(plot.title = element_text(color="black", face="italic",hjust=0))
      g <- g + theme(legend.direction = 'horizontal', 
                     legend.position = 'bottom',
                     legend.background = element_rect(fill = "white", colour = 'black'),
                     legend.text = element_text(colour = 'black'),
                     legend.key = element_rect(fill="white",colour="white"))
      g = g + theme(panel.background = element_rect(fill = "white", colour = 'grey'))
      g = g + theme(plot.background = element_rect(fill = "white", colour = 'grey'),
                    panel.grid = element_blank())
      g = g + theme(axis.text.x=element_text(colour="black"),
                    axis.text.y=element_text(colour="black"))
      g = g + theme(plot.title = element_text(color="black"),
                    axis.title.x = element_text(color="black"),
                    axis.title.y = element_text(color="black"))
      g = g + theme(line = element_line(color="black"),
                    rect= element_rect(color="black"),
                    text=element_text(color="black"))
      g = g + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
      
      print(g)
      
      #print(g)
    }
    ggsave("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/Analyses and Images/PCA/morphology_pca_PC2PC3_black.png", 
           g, bg = "transparent",width=7,height=5,units = "in")
  }
  
  png("all_pcs_plot.png",width=900,height=400); {
  par(mfrow=c(1,4),mar=c(4,4,0.2,0.2))
    m <- rbind(c(1,1,1), c(2,3,4), c(2,3,4), c(2,3,4), c(2,3,4), c(2,3,4))
    #print(m)
    layout(m)
    #layout.show(4)
  palette(c("black","orange","goldenrod","green","pink",
            "magenta","purple","blue","grey","red"))
  plot.new()
  legend("center",legend=unique(axes_full$SPP),ncol=5,bty="n",cex=1.7,pt.cex=2,
         col=as.factor(axes_full$SPP),pch=as.numeric(unique(axes_full$SPP)))
  
  plot(axes_full$PC2,axes_full$PC1,col=as.factor(axes_full$SPP),
       pch=as.numeric(as.factor(axes_full$SPP)),
       xlim=c(-3,5),ylab="PC1",xlab="PC2")
  for(i in 1:length(unique(axes_full$SPP))) {
    spp = sort(unique(axes_full$SPP))[i]
    toplot = axes_full[axes_full$SPP==spp,c("PC2","PC1")]
    toplot = toplot[complete.cases(toplot),]
    car::dataEllipse(toplot$PC2, toplot$PC1, levels=c(0.75),lty=1,
                     add=T,col=i,plot.points=F,center.pch = F)
  }
  plot(axes_full$PC3,axes_full$PC1,col=as.factor(axes_full$SPP),
       pch=as.numeric(as.factor(axes_full$SPP)),
       xlim=c(-3,3),ylab="PC1",xlab="PC3")
  for(i in 1:length(unique(axes_full$SPP))) {
    spp = sort(unique(axes_full$SPP))[i]
    toplot = axes_full[axes_full$SPP==spp,c("PC3","PC1")]
    toplot = toplot[complete.cases(toplot),]
    car::dataEllipse(toplot$PC3, toplot$PC1, levels=c(0.75),lty=1,
                     add=T,col=i,plot.points=F,center.pch = F)
  }
  plot(axes_full$PC3,axes_full$PC2,col=as.factor(axes_full$SPP),
       pch=as.numeric(as.factor(axes_full$SPP)),
       xlim=c(-3,3),ylab="PC2",xlab="PC3")
  for(i in 1:length(unique(axes_full$SPP))) {
    spp = sort(unique(axes_full$SPP))[i]
    toplot = axes_full[axes_full$SPP==spp,c("PC3","PC2")]
    toplot = toplot[complete.cases(toplot),]
    car::dataEllipse(toplot$PC3, toplot$PC2, levels=c(0.75),lty=1,
                     add=T,col=i,plot.points=F,center.pch = F)
  }
    }; dev.off()
  
  
  
  
  
  
  plot(axes_full$PC1,axes_full$PC2,col=as.numeric(axes_full$SPP),
       pch=as.numeric(axes_full$SPP))
  unique(axes_full$SPP)
  legend("topright",legend=unique(axes_full$SPP),ncol=2,bty="n",
         col=as.numeric(axes_full$SPP),pch=as.numeric(unique(axes_full$SPP)))
  
  plot(axes_full$PC3,axes_full$PC2,col=as.numeric(axes_full$SPP),
       pch=as.character(axes_full$SPP))
  unique(axes_full$SPP)
  
  plot(axes_notail$PC1,axes_notail$PC2,col=as.numeric(axes_notail$SPP),
       pch=as.character(axes_notail$SPP))
  unique(axes_notail$SPP)
  plot(axes_notail$PC3,axes_notail$PC2,col=as.numeric(axes_notail$SPP),
       pch=as.character(axes_notail$SPP))
  unique(axes_notail$SPP)
  
  palette(c("black","red","cyan"))
  plot(axes_notailsec$PC1,axes_notailsec$PC2,col=as.numeric(axes_notailsec$WHICH.SIDE.OF.CFB),
       pch=as.numeric(axes_notailsec$SPP))
  unique(axes_notailsec$SPP)
  
  palette(c("black","red","orange","goldenrod","green","blue",
            "cyan","magenta","brown","pink","lightblue"))
  plot(axes_notailsec$PC3,axes_notailsec$PC2,col=as.numeric(axes_notailsec$SPP),
       pch=as.character(axes_notailsec$SPP))
  unique(axes_notailsec$SPP)
  
  
  p <- plot_ly(axes_full, x = ~PC1, y = ~PC2, z = ~PC3, color = ~SPP,symbol = ~WHICH.SIDE.OF.CFB, colors = c("black","red","orange","goldenrod","green","blue",
                                                                                                             "cyan","magenta","brown","pink","lightblue"),
               symbols=c(1,16)) %>%
    add_markers() %>%
    layout(scene = list(xaxis = list(title = 'PC1'),
                        yaxis = list(title = 'PC2'),
                        zaxis = list(title = 'PC3')))
  p
  
  corrd = cor(axes_full[,c(1:3,15:28)])
  corrd = corrd[1:3,4:17]
  colnames(corrd) = c("BILL HEIGHT","BILL LENGTH","BILL WIDTH","TARSUS","PRIMARIES",
                      "SECONDARIES","TAIL","KIPPS","BASE AREA","VOLUME","LATERAL AREA","TOTAL AREA","LATERAL VOL RATIO","TOTAL VOL RATIO")
  
  png("corrplot_pcs_morphology.png")
  corrplot::corrplot(corrd,is.corr = F,method="ellipse",diag=T,
                     cl.lim=c(-1,1), col=colorRampPalette(c("cornflowerblue","lightgrey","coral"))(200))
  corrplot::corrplot(corrd,is.corr = F,method="number",diag=T,add=T,bg=rgb(0,0,0,0),cl.pos="n",
                       cl.lim=c(-1,1), col=colorRampPalette(c("black","black","black"))(200))
    dev.off()
  
  corrplot::corrplot(cor(axes_full[,c(1:3,15:28)]),
           method="ellipse",diag=F,order="hclust")
  
  corrplot(cor(axes_notail[,c(1:3,14:26)]),
           method="ellipse",diag=F,type="upper")
  
  corrplot(cor(axes_notailsec[,c(1:3,12:22)]),
           method="ellipse",diag=F,type="upper")
  
  ## more t-tests with pca
  
  for (species in unique(axes_full$SPP)) {
    test_axes = axes_full[axes_full$SPP==species,]
    son_axes = test_axes[test_axes$WHICH.SIDE.OF.CFB=="SONORAN",]
    chi_axes = test_axes[test_axes$WHICH.SIDE.OF.CFB=="CHIHUAHUAN",]
    for (colnum in c(1:3,4:14,15:24)) {
      colname = names(test_axes)[colnum]
      son_data = son_axes[,colnum]
      chi_data = chi_axes[,colnum]
      tested = t.test(son_data,chi_data)
      if (tested$p.value <= 0.05) {
        print(paste("SIGNIF:",species,colname,tested$p.value,sep=" "))
        #} else {
        #print(paste("          NOT:",species,colname,tested$p.value,sep=" "))
      }
    }
  }
  
  for (species in sample(unique(agg_res$SPP))) {
    print(species)
    test_axes = agg_res[agg_res$SPP==species,]
    son_axes = test_axes[test_axes$WHICH.SIDE.OF.CFB=="SONORAN",]
    chi_axes = test_axes[test_axes$WHICH.SIDE.OF.CFB=="CHIHUAHUAN",]
    #for (colnum in c(28:33,41:47)) { ## skip 34 for nitens and cris
    for (colnum in c(34)) { ## skip 34 for nitens and cris
      #print(colnum)
      colname = names(test_axes)[colnum]
      son_data = son_axes[,colnum]
      chi_data = chi_axes[,colnum]
      tested = t.test(son_data,chi_data)
      if (tested$p.value <= 0.05) {
        print(paste("SIGNIF:",species,colname,tested$p.value,sep=" "))
        #} else {
        #print(paste("          NOT:",species,colname,tested$p.value,sep=" "))
      }
    }
  }
}; dev.off()

## do a distance matrix of morphology? 
print("DIST"); {
  distances = dist(axes_full[axes_full$SPP=="BELLII",c("PC1","PC2","PC3")])
  plot(distances)
}; dev.off()

## discriminant function analysis
print("DFA"); {
  
  
  together = merge(agg[,c(2,4:12,24,37,41:47)],agg_res[,c(2,28:47)],by="CATALOG.NUMBER")
  together=together[together$SPP!="CARDINALIS",]
  
  
  
  for_dfa = together[,c(11:12,2:8,13:39)]
  for_dfa = for_dfa[complete.cases(for_dfa),]
  
  set.seed(123)
  training.samples <- for_dfa$SPP %>%
    createDataPartition(p = 0.8, list = FALSE)
  train.data <- for_dfa[training.samples, ]
  test.data <- for_dfa[-training.samples, ]
  
  # Estimate preprocessing parameters
  preproc.param <- train.data %>% 
    preProcess(method = c("center", "scale"))
  # Transform the data using the estimated parameters
  train.transformed <- preproc.param %>% predict(train.data)
  test.transformed <- preproc.param %>% predict(test.data)
  
  # Fit the model
  model <- lda(SPP~., data = train.transformed)
  model
  #plot(model)
  # Make predictions
  predictions <- model %>% predict(test.transformed)
  # Model accuracy
  mean(predictions$class==test.transformed$SPP)
  lda.data <- cbind(train.transformed, predict(model)$x)
  png("Analyses and Images/DFA/dfa_all_spp.png")
  g = ggplot(lda.data, aes(LD1, LD2)) +
    geom_point(aes(color = SPP))
  #print(g)
  dev.off()
  
  ## now do dfa per single species 
  par(ask=F)
  
  for (spp in sort(unique(together$SPP))) {
    png(paste("Analyses and Images/DFA/dfa_single",spp,".png",sep=""))
    print(spp)
    
    for_dfa = together[together$SPP==spp,c(11:12,2:8,13:39)]
    for_dfa = for_dfa[complete.cases(for_dfa),]
    
    set.seed(123)
    training.samples <- for_dfa$WHICH.SIDE.OF.CFB %>%
      createDataPartition(p = 0.8, list = FALSE)
    train.data <- for_dfa[training.samples, ]
    test.data <- for_dfa[-training.samples, ]
    
    # Estimate preprocessing parameters
    preproc.param <- train.data %>% 
      preProcess(method = c("center", "scale"))
    # Transform the data using the estimated parameters
    train.transformed <- preproc.param %>% predict(train.data)
    test.transformed <- preproc.param %>% predict(test.data)
    
    # Fit the model
    model <- lda(WHICH.SIDE.OF.CFB~., data = train.transformed)
    model
    #plot(model)
    # Make predictions
    predictions <- model %>% predict(test.transformed)
    # Model accuracy
    mean(predictions$class==test.transformed$WHICH.SIDE.OF.CFB)
    lda.data <- cbind(train.transformed, predict(model)$x)
    g = ggplot(lda.data, aes(LD1, LD1)) +
      geom_point(aes(color = WHICH.SIDE.OF.CFB))
    #print(g)
    dev.off()
    
  }
}

#### dispersal and thermoregulation
print("DISPERSAL") {
pal = colorRampPalette(c("red","blue"))
colors = findInterval(agg$KIPPSINDEX_CALC, sort(agg$KIPPSINDEX_CALC))
cols=pal(nrow(agg))[colors]
plot(agg$WING.LENGTH.TO.PRIMARIES,agg$WING.LENGTH.TO.SECONDARIES,
     col=cols,pch=as.numeric(agg$SPP))
abline(a=0,b=1)

png("Analyses and Images/DISPERSE-THERMO/primaries_vs_kipps.png"); {
  palette(c("red","orange","goldenrod","green","cyan",
            "blue","purple","magenta","brown","lightblue",
            "pink"))
  par(
    bg = NA,
    col.axis = "white",
    fg = "white",
    col.lab = "white",
    col.main = "white"
  )
  plot(agg$WING.LENGTH.TO.PRIMARIES,agg$KIPPSINDEX_CALC,
       col=as.numeric(agg$SPP),pch=as.numeric(agg$WHICH.SIDE.OF.CFB),
       type="n",ylab="Kipps Index",xlab="Primary Length (mm)")
  for(i in 1:length(unique(agg$SPP))) {
    spp = sort(unique(agg$SPP))[i]
    toplot = agg[agg$SPP==spp,c("WING.LENGTH.TO.PRIMARIES","KIPPSINDEX_CALC")]
    toplot = toplot[complete.cases(toplot),]
    car::dataEllipse(toplot$WING.LENGTH.TO.PRIMARIES, toplot$KIPPSINDEX_CALC, levels=c(0.5),
                     add=T,col=i,plot.points=F)
  }
  legend("top",legend=sort(unique(agg$SPP)),ncol=2,
         fill=as.numeric(sort(unique(agg$SPP))),bty="n",
         border="white")
}; dev.off()

png("Analyses and Images/DISPERSE-THERMO/primaries_vs_secondaries.png"); {
  palette(c("red","orange","goldenrod","green","cyan",
            "blue","purple","magenta","brown","lightblue",
            "pink"))
  par(
    bg = NA,
    col.axis = "white",
    fg = "white",
    col.lab = "white",
    col.main = "white"
  )
  plot(agg$WING.LENGTH.TO.PRIMARIES,agg$WING.LENGTH.TO.SECONDARIES,
       col=as.numeric(agg$SPP),pch=as.numeric(agg$WHICH.SIDE.OF.CFB),
       type="n",ylab="Secondary Length (mm)",xlab="Primary Length (mm)")
  abline(a=0,b=1,lty=2)
  for(i in 1:length(unique(agg$SPP))) {
    spp = sort(unique(agg$SPP))[i]
    toplot = agg[agg$SPP==spp,c("WING.LENGTH.TO.PRIMARIES","WING.LENGTH.TO.SECONDARIES")]
    toplot = toplot[complete.cases(toplot),]
    car::dataEllipse(toplot$WING.LENGTH.TO.PRIMARIES, toplot$WING.LENGTH.TO.SECONDARIES, levels=c(0.5),
                     add=T,col=i,plot.points=F)
  }
  legend("bottomright",legend=sort(unique(agg$SPP)),ncol=2,
         fill=as.numeric(sort(unique(agg$SPP))),bty="n",
         border="white")
}; dev.off()

png("Analyses and Images/DISPERSE-THERMO/tail_vs_kipps.png",width=750,height=480); {
  palette(c("red","orange","goldenrod","green","cyan",
            "blue","purple","magenta","brown","lightblue",
            "pink"))
  par(
    bg = NA,
    col.axis = "white",
    fg = "white",
    col.lab = "white",
    col.main = "white"
  )
  plot(agg$TAIL,agg$KIPPSINDEX_CALC,
       col=as.numeric(agg$SPP),pch=as.numeric(agg$WHICH.SIDE.OF.CFB),
       type="n",ylab="Kipps Index",xlab="Tail Length (mm)")
  for(i in 1:length(unique(agg$SPP))) {
    spp = sort(unique(agg$SPP))[i]
    if (spp != "CARDINALIS") {
      toplot = agg[agg$SPP==spp,c("TAIL","KIPPSINDEX_CALC")]
      toplot = toplot[complete.cases(toplot),]
      car::dataEllipse(toplot$TAIL, toplot$KIPPSINDEX_CALC, levels=c(0.5),
                       add=T,col=i,plot.points=F)
    }
  }
  legend("top",legend=sort(unique(agg$SPP)),ncol=2,
         fill=as.numeric(sort(unique(agg$SPP))),bty="n",
         border="white")
}; dev.off()

png("Analyses and Images/DISPERSE-THERMO/tail_vs_primaries.png"); {
  palette(c("red","orange","goldenrod","green","cyan",
            "blue","purple","magenta","brown","lightblue",
            "pink"))
  par(
    bg = NA,
    col.axis = "white",
    fg = "white",
    col.lab = "white",
    col.main = "white"
  )
  plot(agg$TAIL,agg$WING.LENGTH.TO.PRIMARIES,
       col=as.numeric(agg$SPP),pch=as.numeric(agg$WHICH.SIDE.OF.CFB),
       type="n",ylab="Primary Length (mm)",xlab="Tail Length (mm)")
  abline(a=0,b=1,lty=2)
  
  for(i in 1:length(unique(agg$SPP))) {
    spp = sort(unique(agg$SPP))[i]
    if (spp != "CARDINALIS") {
      toplot = agg[agg$SPP==spp,c("TAIL","WING.LENGTH.TO.PRIMARIES")]
      toplot = toplot[complete.cases(toplot),]
      car::dataEllipse(toplot$TAIL, toplot$WING.LENGTH.TO.PRIMARIES, levels=c(0.5),
                       add=T,col=i,plot.points=F)
    }
  }
  legend("bottomright",legend=sort(unique(agg$SPP)),ncol=2,
         fill=as.numeric(sort(unique(agg$SPP))),bty="n",
         border="white")
}; dev.off()

print("THERMOREG")
png("Analyses and Images/DISPERSE-THERMO/beaksavr_vs_tarsus.png",width=750,height=480); {
  palette(c("red","orange","goldenrod","green","cyan",
            "blue","purple","magenta","brown","lightblue",
            "pink"))
  par(
    bg = NA,
    col.axis = "white",
    fg = "white",
    col.lab = "white",
    col.main = "white"
  )
  plot(agg$TARSUS.LENGTH,agg$BEAKAREAVOLUME_CALC,
       col=as.numeric(agg$SPP),pch=as.numeric(agg$WHICH.SIDE.OF.CFB),
       type="n",ylab="Beak Surface Area-Volume Ratio",xlab="Tarsus Length (mm)")
  for(i in 1:length(unique(agg$SPP))) {
    spp = sort(unique(agg$SPP))[i]
    toplot = agg[agg$SPP==spp,c("TARSUS.LENGTH","BEAKAREAVOLUME_CALC")]
    toplot = toplot[complete.cases(toplot),]
    car::dataEllipse(toplot$TARSUS.LENGTH, toplot$BEAKAREAVOLUME_CALC, levels=c(0.5),
                     add=T,col=i,plot.points=F)
  }
  legend("topright",legend=sort(unique(agg$SPP)),ncol=1,
         fill=as.numeric(sort(unique(agg$SPP))),bty="n",
         border="white")
}; dev.off()

png("Analyses and Images/DISPERSE-THERMO/beakvol_vs_tarsus.png"); {
  palette(c("red","orange","goldenrod","green","cyan",
            "blue","purple","magenta","brown","lightblue",
            "pink"))
  par(
    bg = NA,
    col.axis = "white",
    fg = "white",
    col.lab = "white",
    col.main = "white"
  )
  plot(agg$TARSUS.LENGTH,agg$BEAKVOL_CALC,
       col=as.numeric(agg$SPP),pch=as.numeric(agg$WHICH.SIDE.OF.CFB),
       type="n",ylab="Beak Volume",xlab="Tarsus Length (mm)")
  for(i in 1:length(unique(agg$SPP))) {
    spp = sort(unique(agg$SPP))[i]
    toplot = agg[agg$SPP==spp,c("TARSUS.LENGTH","BEAKVOL_CALC")]
    toplot = toplot[complete.cases(toplot),]
    car::dataEllipse(toplot$TARSUS.LENGTH, toplot$BEAKVOL_CALC, levels=c(0.5),
                     add=T,col=i,plot.points=F)
  }
  legend("topleft",legend=sort(unique(agg$SPP)),ncol=1,
         fill=as.numeric(sort(unique(agg$SPP))),bty="n",
         border="white")
}; dev.off()

png("Analyses and Images/DISPERSE-THERMO/beakvlatsurf_vs_tarsus.png"); {
  palette(c("red","orange","goldenrod","green","cyan",
            "blue","purple","magenta","brown","lightblue",
            "pink"))
  par(
    bg = NA,
    col.axis = "white",
    fg = "white",
    col.lab = "white",
    col.main = "white"
  )
  plot(agg$TARSUS.LENGTH,agg$BEAKLATERALSURFACE_CALC,
       col=as.numeric(agg$SPP),pch=as.numeric(agg$WHICH.SIDE.OF.CFB),
       type="n",ylab="Beak Surface Area",xlab="Tarsus Length (mm)")
  for(i in 1:length(unique(agg$SPP))) {
    spp = sort(unique(agg$SPP))[i]
    toplot = agg[agg$SPP==spp,c("TARSUS.LENGTH","BEAKLATERALSURFACE_CALC")]
    toplot = toplot[complete.cases(toplot),]
    car::dataEllipse(toplot$TARSUS.LENGTH, toplot$BEAKLATERALSURFACE_CALC, levels=c(0.5),
                     add=T,col=i,plot.points=F)
  }
  legend("bottomright",legend=sort(unique(agg$SPP)),ncol=1,
         fill=as.numeric(sort(unique(agg$SPP))),bty="n",
         border="white")
}; dev.off()

print("BOTH")

png("Analyses and Images/DISPERSE-THERMO/beaksavr_vs_kipps.png"); {
  palette(c("red","orange","goldenrod","green","cyan",
            "blue","purple","magenta","brown","lightblue",
            "pink"))
  par(
    bg = NA,
    col.axis = "white",
    fg = "white",
    col.lab = "white",
    col.main = "white"
  )
  plot(agg$BEAKAREAVOLUME_CALC,agg$KIPPSINDEX_CALC,
       col=as.numeric(agg$SPP),pch=as.numeric(agg$WHICH.SIDE.OF.CFB),
       type="n",xlab="Beak Surface-Area Volume Ratio",ylab="Kipp's Index")
  for(i in 1:length(unique(agg$SPP))) {
    spp = sort(unique(agg$SPP))[i]
    toplot = agg[agg$SPP==spp,c("KIPPSINDEX_CALC","BEAKAREAVOLUME_CALC")]
    toplot = toplot[complete.cases(toplot),]
    car::dataEllipse(toplot$BEAKAREAVOLUME_CALC, toplot$KIPPSINDEX_CALC, levels=c(0.5),
                     add=T,col=i,plot.points=F)
  }
  legend("topright",legend=sort(unique(agg$SPP)),ncol=2,
         fill=as.numeric(sort(unique(agg$SPP))),bty="n",
         border="white")
}; dev.off()

png("Analyses and Images/DISPERSE-THERMO/beaksavr_vs_tail.png"); {
  palette(c("red","orange","goldenrod","green","cyan",
            "blue","purple","magenta","brown","lightblue",
            "pink"))
  par(
    bg = NA,
    col.axis = "white",
    fg = "white",
    col.lab = "white",
    col.main = "white"
  )
  plot(agg$TAIL,agg$BEAKAREAVOLUME_CALC,
       col=as.numeric(agg$SPP),pch=as.numeric(agg$WHICH.SIDE.OF.CFB),
       type="n",ylab="Beak Surface-Area Volume Ratio",xlab="Tail Length (mm)")
  for(i in 1:length(unique(agg$SPP))) {
    spp = sort(unique(agg$SPP))[i]
    print(i)
    toplot = agg[agg$SPP==spp,c("TAIL","BEAKAREAVOLUME_CALC")]
    toplot = toplot[complete.cases(toplot),]
    car::dataEllipse(toplot$TAIL,toplot$BEAKAREAVOLUME_CALC, levels=c(0.5),
                     add=T,col=i,plot.points=F)
  }
  legend("topright",legend=sort(unique(agg$SPP)),ncol=1,
         fill=as.numeric(sort(unique(agg$SPP))),bty="n",
         border="white")
}; dev.off()


#hist(abs(agg$BEAKTOTAREAVOLUME-agg$BEAKTOTAREAVOLUME_CALC))
#quantile(abs(agg$BEAKTOTAREAVOLUME-agg$BEAKTOTAREAVOLUME_CALC),seq(0.05,0.95,0.05),na.rm=T)
#unique(agg[which(abs(agg$BEAKTOTAREAVOLUME-agg$BEAKTOTAREAVOLUME_CALC)>0.0006),"CATALOG.NUMBER"])


## brif islands detour
names = c("bel","cri","cur","fla","fus",
          "bil","bru","mel","nit","sin")
names2 = c("V","c","T","f","M","b","A","m","n","s")
islan = c(25,49,14,92,212,
          257,72,123,72,78)
sweep = c(127,127,200,123,72,
          55,66,44,85,69)
str=c(1,1,1,1,1,0,0,0,0,0)

x = data.frame(names,names2,islan,sweep,str)
rownames(x) = names

mean_vals = aggregate(cbind(BILL.HEIGHT,BILL.LENGTH,BILL.WIDTH,TARSUS.LENGTH,WING.LENGTH.TO.PRIMARIES,
                            WING.LENGTH.TO.SECONDARIES,TAIL,KIPPSINDEX_CALC,BEAKBASEAREA_CALC,
                            BEAKVOL_CALC,BEAKLATERALSURFACE_CALC,BEAKTOTALSURFACE_CALC,BEAKAREAVOLUME_CALC)~SPP,data=agg,FUN=mean)

mean_vals2 = aggregate(cbind(BEAKAREAVOLUME_CALC,BEAKBASEAREA_CALC,BEAKLATERALSURFACE_CALC,BEAKTOTALSURFACE_CALC,
                             BEAKTOTAREAVOLUME_CALC,BEAKVOL_CALC,BILL.HEIGHT,BILL.LENGTH,BILL.WIDTH,KIPPSINDEX_CALC,
                             PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,
                             TAIL,TARSUS.LENGTH,WING.LENGTH.TO.PRIMARIES,WING.LENGTH.TO.SECONDARIES)~SPP,data=axes_full,FUN=mean)

mean_vals3 = aggregate(cbind(RES_BILL.HEIGHT,RES_BILL.LENGTH,RES_BILL.WIDTH,RES_TARSUS.LENGTH,RES_WING.LENGTH.TO.PRIMARIES,
                             RES_WING.LENGTH.TO.SECONDARIES,RES_TAIL,RES_KIPPSINDEX_CALC,RES_BEAKBASEAREA_CALC,
                             RES_BEAKVOL_CALC,RES_BEAKLATERALSURFACE_CALC,RES_BEAKTOTALSURFACE_CALC,RES_BEAKAREAVOLUME_CALC)~SPP,data=agg_res,FUN=mean)


rownames(mean_vals2) = sort(names)

together =merge(x,mean_vals2,by=0)
cors = cor(together[,c(4:6,8:20,32:35)])

corrplot(cors[c(1:3),c(4:20)],
         method="ellipse",diag=T)

png("Analyses and Images/GENOMES/morphology_sweeps_correlated.png",
    width=600,height=450)
palette(c("red","cyan"))
par(mfrow=c(2,3))
par(
  bg = NA,
  col.axis = "white",
  fg = "white",
  col.lab = "white",
  col.main = "white"
)
par(mar=c(4.5,4,0,0))
plot(together$sweep,together$PC3,pch=together$str+16,col=together$str+1,cex=2,xlab="Number Sweeps",ylab="Mean PC3")
plot(together$sweep,together$BILL.LENGTH,pch=together$str+16,col=together$str+1,cex=2,xlab="Number Sweeps",ylab="Mean Bill Length")
plot(together$islan,together$BILL.LENGTH,pch=together$str+16,col=together$str+1,cex=2,xlab="Number Islands",ylab="Mean Bill Length")
plot(together$sweep,together$BEAKLATERALSURFACE_CALC,pch=together$str+16,col=together$str+1,cex=2,xlab="Number Sweeps",ylab="Mean Beak Surface Area")
plot(together$sweep,together$TARSUS.LENGTH,pch=together$str+16,col=together$str+1,cex=2,xlab="Number Sweeps",ylab="Mean Tarsus Length")
boxplot(together$BILL.LENGTH~together$str,col=c("red","cyan"),cex=2,names=c("No","Yes"),xlab="Structure",ylab="Mean Bill Length")
dev.off()


plot(sweep,islan,pch=str+16,col=str+1,cex=2,
     xlab="Number Sweeps",
     ylab="Number Islands") # plot with body size on x-axis and survival (0 or 1) on y-axis
#points(islan, predict(modelC, newdata=x, type="response"),pch=21,bg=as.numeric(str)+2)
legend("topright",legend=c("Structure",
                           "No Structure"),
       col=c("cyan","red"),
       pch=c(16,17),cex=2,
       border="white",bg="grey")


## calculate zscores for pc1-pc3
by1 = axes_full$SPP
by2 = axes_full$WHICH.SIDE.OF.CFB
bys = list(by1,by2)
attach(agg_res)
attach(axes_full)
mean_vals_bydes = aggregate(x=cbind(PC1,PC2,PC3),by=bys,FUN=function(x) {
  y = mean(x,na.rm=T)
  return(y)
})

scaled_full = data.frame(matrix(ncol=5,nrow=0))
colnames(scaled_full) = c("SPP","WHICH.SIDE.OF.CFB",
                          "PC1","PC2","PC3")
for (spp in unique(axes_full$SPP)) {
  sppdf = axes_full[axes_full$SPP==spp,c("SPP","WHICH.SIDE.OF.CFB","PC1","PC2","PC3")]
  
  toscale = sppdf[,c("PC1","PC2","PC3")]
  scaledf = data.frame(scale(toscale))
  scaledf$SPP = substr(as.character(sppdf$SPP),1,3)
  scaledf$WHICH.SIDE.OF.CFB = substr(as.character(sppdf$WHICH.SIDE.OF.CFB),1,3)
  scaled_full = rbind(scaled_full,scaledf)
  
}

scaled_full$POP = paste(scaled_full$SPP,scaled_full$WHICH.SIDE.OF.CFB,sep="-")
# by1 = scaled_full$SPP
# by2 = scaled_full$WHICH.SIDE.OF.CFB
# bys = list(by1,by2)
head(scaled_full)
mean_vals_bydes_scale = aggregate(x=cbind(scaled_full$PC1,scaled_full$PC2,scaled_full$PC3),by=list(scaled_full$POP),FUN=function(x) {
  y = mean(x,na.rm=T)
  return(y)
})
names(mean_vals_bydes) = c("POP","PC1","PC2","PC3") 
names(mean_vals_bydes_scale) = c("POP","PC1","PC2","PC3") 

mean_vals_bydes_scale = mean_vals_bydes_scale[order(mean_vals_bydes_scale$POP),]
mean_vals_bydes = mean_vals_bydes[order(mean_vals_bydes_scale$POP),]

mean_vals_bydes_scale_chi = as.matrix(mean_vals_bydes_scale[seq(1,nrow(mean_vals_bydes_scale),2),])
mean_vals_bydes_scale_son = as.matrix(mean_vals_bydes_scale[seq(2,nrow(mean_vals_bydes_scale),2),])

pc1dif = as.numeric(mean_vals_bydes_scale_chi[,"PC1"])-as.numeric(mean_vals_bydes_scale_son[,"PC1"])
pc2dif = as.numeric(mean_vals_bydes_scale_chi[,"PC2"])-as.numeric(mean_vals_bydes_scale_son[,"PC2"])
pc3dif = as.numeric(mean_vals_bydes_scale_chi[,"PC3"])-as.numeric(mean_vals_bydes_scale_son[,"PC3"])
difmatrix = as.data.frame(cbind(pc1dif,pc2dif,pc3dif))
absmatrix = abs(difmatrix)

difmatrix$species=c("bel","bil","bru","cri","cur",
                    "fla","fus","mel","nit","sin")
absmatrix$species=c("bel","bil","bru","cri","cur",
                    "fla","fus","mel","nit","sin")


barplot(absmatrix[,1])

}; dev.off();

