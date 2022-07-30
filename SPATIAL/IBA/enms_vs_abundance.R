library(raster)

abd_path = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/BBS_data/abundance_asciis/clipped/"

abd_list = c("BELLII_Idw_abundance_clipped.asc",
             "CURVIROSTRE_Idw_abundance_clipped.asc",
             "CRISSALE_Idw_abundance_clipped.asc",
             "MELANURA_Idw_abundance_clipped.asc",
             "NITENS_Idw_abundance_clipped.asc",
             "FUSCA_Idw_abundance_clipped.asc",
             "SINUATUS_Idw_abundance_clipped.asc",
             "BRUNNEICAPILLUS_Idw_abundance_clipped.asc",
             "FLAVICEPS_Idw_abundance_clipped.asc",
             "BILINEATA_Idw_abundance_clipped.asc")

enm_list = c("Vireo/AddedPredThin_BestModel_Vireo bellii_addedLayers_USA.asc_-4.56468915939331_-18.812442779541.rescaled.asc",
             "Toxostoma/cur_full/AddedPredThin_BestModel_Toxostoma curvirostre_addedLayers_USA.asc_-5.62186193466187_-12.8309497833252.rescaled.asc",
             "Toxostoma/crissale/AddedPredThin_BestModel_Toxostoma crissale_addedLayers_USA.asc_-4.63216781616211_-17.7025375366211.rescaled.asc",
             "Polioptila/AddedPredThin_BestModel_Polioptila melanura_addedLayers_USA.asc_-3.93294405937195_-21.4272289276123.rescaled.asc",
             "Phainopepla/AddedPredThin_BestModel_Phainopepla nitens_addedLayers_USA.asc_-3.93656158447266_-22.6519756317139.rescaled.asc",
             "Melozone/AddedPredThin_BestModel_Melozone fusca_addedLayers_USA.asc_-5.73480558395386_-13.4299039840698.rescaled.asc",
             "Cardinalis/AddedPredThin_BestModel_Cardinalis sinuatus_addedLayers_USA.asc_-5.24065685272217_-19.4379730224609.rescaled.asc",
             "Campylorhynchus/AddedPredThin_BestModel_Campylorhynchus brunneicapillus_addedLayers_USA.asc_-5.1320424079895_-18.1571197509766.rescaled.asc",
             "Auriparus/AddedPredThin_BestModel_Auriparus flaviceps_addedLayers_USA.asc_-4.34592342376709_-18.2621974945068.rescaled.asc",
             "Amphispiza/AddedPredThin_BestModel_Amphispiza bilineata_addedLayers_USA.asc_-5.90643692016602_-13.0102529525757.rescaled.asc")

for (i in 1:10) {
print(i)
abd = raster(paste(abd_path,abd_list[i],sep=""))
enm = raster(paste(enm_path,enm_list[i],sep=""))

print(abd)
print(enm)

ext_a = extent(abd)
ext_e = extent(enm)

abd_c = crop(abd,ext_e)

enm_c = aggregate(enm,20)
#abd_cr = disaggregate(abd_c,20)

abd_val = values(abd_c)
enm_val = values(enm_c)

val = data.frame(abd_val,enm_val)
val = unique(val)
val = val[complete.cases(val),]
names(val) = c("abd","enm")

nrow(val)
cor(val)

subone = val[val$abd < 1,]
moreone = val[val$abd >= 1,]

#subonequan = quantile(subone$enm,probs=c(0.25,0.50,0.75))
moreonequan = quantile(moreone$enm,probs=c(0.25,0.50,0.75))

always = max(subone$enm)

never = min(moreone$enm)

#png(paste(abd_path,abd_list[i],".HIST.png",sep=""))
# hist(val$enm,main=abd_list[i],xlab="Environmental Value")
# abline(v=always,col="red")
# abline(v=never,col="blue")
#dev.off()

#png(paste(abd_path,abd_list[i],".RASTER.png",sep=""))
# breakpoints <- c(0,never,moreonequan,always,1)
# colors <- c("black","red","pink","yellow","lightblue","blue")
# plot(enm,breaks=breakpoints,col=colors,main=abd_list[i])
# legend("bottomleft",legend=c("Never","<25%","25-50%","50-75%","75+%","Always"),
#        fill=colors)
#dev.off()

#png(paste(abd_path,abd_list[i],".png",sep=""))
# plot(val$abd,val$enm,col=rgb(0,0,0,0.1),
#      main=abd_list[i],xlab="Abundance",
#      ylab="Habitat Suitability")
# mod = lm(val$enm~sqrt(sqrt(val$abd)))
# summary(mod)
# rsq = summary(mod)$adj.r.squared
# newx = (val$abd)
# newy = predict(mod,data.frame(abd=newx))
# points(newx,newy, col=rgb(1,0,0,0.1))
# legend("bottomright",legend=c(paste("Rsq: ",round(rsq,2),sep=""),
#                               paste("Always: ",round(always,2),sep=""),
#                               paste("Never: ",round(never,2),sep="")),
#        bty="n")
# abline(h=always,col="blue",lty=2,lwd=2)
# abline(h=never,col="blue",lty=2,lwd=2)
#dev.off()

png(paste(abd_path,abd_list[i],"_LOG.png",sep=""))
plot(log10(val$abd),log10(val$enm),col=rgb(0,0,0,0.2),
     main=abd_list[i],xlab="Log Abundance",
     ylab="Log Habitat Suitability")
mod = lm(log10(val$enm)~log10(val$abd))
summary(mod)
newx = log10(val$abd)
newy = predict(mod,data.frame(abd=newx))
points(10^newx,10^newy, col=rgb(1,0,0,0.1))
abline(mod,col="red",lwd=3)
rsq = summary(mod)$adj.r.squared
legend("bottomright",legend=c(paste("Rsq: ",round(rsq,2),sep=""),
                              paste("Always: ",round(always,2),sep=""),
                              paste("Never: ",round(never,2),sep="")),
       bty="n")
abline(h=log10(always),col="blue",lty=2,lwd=2)
abline(h=log10(never),col="blue",lty=2,lwd=2)
dev.off()
png(paste(abd_path,abd_list[i],"_LOGMODEL_notlogplot.png",sep=""))
plot(val$abd,val$enm,col=rgb(0,0,0,0.2),
     main=abd_list[i],xlab="Abundance",
     ylab="Habitat Suitability")
newx = log10(val$abd)
newy = predict(mod,data.frame(abd=newx))
points(10^newx,10^newy, col=rgb(1,0,0,0.1))
#abline(mod,col="red",lwd=3)
#rsq = summary(mod)$adj.r.squared
legend("bottomright",legend=c(paste("Rsq: ",round(rsq,2),sep=""),
                              paste("Always: ",round(always,2),sep=""),
                              paste("Never: ",round(never,2),sep="")),
       bty="n")
abline(h=(always),col="blue",lty=2,lwd=2)
abline(h=(never),col="blue",lty=2,lwd=2)
dev.off()



}

## abundance through space

enm = raster(paste(enm_path,enm_list[i],sep=""))
enm_c = aggregate(enm,20)
ext_e = extent(enm)

for (i in 1:10) {
  print(i)
  abd2 = raster(paste(abd_path,abd_list[i],sep=""))
  
  print(abd2)
  #plot(abd)
  
  print(enm)

  abd_c = crop(abd2,ext_e)
  
  values(abd_c)[is.na(values(enm_c))] = NA
  

  ## going to do one with all data, one with 
  
    abd_c2 = crop(abd_c,extent(-120,-97,24,35))
    #abd_c2s = abd_c2
    #values(abd_c2s) = scales::rescale(values(abd_c2s),to=c(0,1))

    df = as.matrix(abd_c2)
    latmeans=colMeans(df,na.rm=T)
  latsds = matrixStats::colSds(df,na.rm=T)

  spts <- as.data.frame(rasterToPoints(abd_c2, spatial = TRUE))
  spts = spts[,c("x","y")]
  spts = unique(spts)
  lats = unique(spts[,"x"])
  longs = unique(spts[,"y"])

  min1=min((latmeans-latsds),na.rm=T)
  max1=max((latmeans+latsds),na.rm=T)

  png(paste("~/",basename(abd_list[i]),"_ABUN.png",sep=""))
  plot(lats,latmeans,type="l",ylim=c(0,max1),
      ylab="Abundance",xlab="Longitude")
  lines(lats,latmeans+latsds,type="l",col="red")
  lines(lats,latmeans-latsds,type="l",col="red")
  dev.off()
  
} 


enm_path2 = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/OUTPUTS/3.enms/WITH_MEXICO/"

enm_list2 = c("Vireo/AddedPredThin_BestModel_Vireo bellii_addedLAyers.asc_-5.37957525253296_-21.9553413391113.rescaled.asc",
             "Toxostoma/AddedPredThin_BestModel_Toxostoma curvirostre_addedLAyers.asc_-5.85683298110962_-20.8998832702637.rescaled.asc",
             "Toxostoma/ascii/AddedPredThin_BestModel_Toxostoma crissale_addedLAyers.asc_-2.73102569580078_-44.6700744628906.rescaled.asc",
             "Polioptila/AddedPredThin_BestModel_Polioptila melanura_addedLAyers.asc_-3.06035566329956_-54.4188537597656.rescaled.asc",
             "Phainopepla/AddedPredThin_BestModel_Phainopepla nitens_addedLAyers.asc_-2.63155722618103_-64.5665817260742.rescaled.asc",
             "Melozone/AddedPredThin_BestModel_Melozone fusca_addedLAyers.asc_-5.53850936889648_-20.7584438323975.rescaled.asc",
             "Cardinalis/AddedPredThin_BestModel_Cardinalis sinuatus_addedLAyers.asc_-5.41119766235352_-32.4833374023438.rescaled.asc",
             "Campylorhynchus/AddedPredThin_BestModel_Campylorhynchus brunneicapillus_addedLAyers.asc_-4.57353591918945_-82.043083190918.rescaled.asc",
             "Auriparus/AddedPredThin_BestModel_Auriparus flaviceps_addedLAyers.asc_-4.41374969482422_-73.6424255371094.rescaled.asc",
             "Amphispiza/AddedPredThin_BestModel_Amphispiza bilineata_addedLAyers.asc_-3.78692746162415_-58.8649635314941.rescaled.asc")

#pdf("~/niche_suitability_2dplot_allspecies.pdf")
pdf("~/niche_suitability_2dplot_zoomedin_allspecies.pdf")
palette(c("black","green","goldenrod","magenta","darkred","cyan","grey","orange","blue","red"))
#plot(0,ylim=c(0,1),xlim=c(-120,-97),type="n",ylab="Scaled Relative Suitability",xlab="Longitude")
plot(0,ylim=c(0,1),xlim=c(-114,-100),type="n",ylab="Scaled Relative Suitability",xlab="Longitude")

enm_means_sds = NULL

library(raster)
for (i in 1:10) {
  print(i)
  enm2 = raster(paste(enm_path2,enm_list2[i],sep=""))
  
  print(enm2)
  
  ext_e2 = extent(enm2)
  
  ## going to do one with all data, one with 
  
#   enm_c2 = crop(enm2,extent(-120,-97,24,35))
#   enm_c2s = enm_c2
#   values(enm_c2s) = scales::rescale(values(enm_c2s),to=c(0,1))
#   
#   df = as.matrix(enm_c2s)
#   latmeans=colMeans(df,na.rm=T)
# latsds = matrixStats::colSds(df,na.rm=T)
# 
# spts <- as.data.frame(rasterToPoints(enm_c2s, spatial = TRUE)) 
# spts = spts[,c("x","y")]
# spts = unique(spts)  
# lats = unique(spts[,"x"])
# longs = unique(spts[,"y"])
# 
# min1=min((latmeans-latsds),na.rm=T)
# max1=max((latmeans+latsds),na.rm=T)
# 
# lines(lats,latmeans,type="l",ylim=c(0,1),col=i)
  #png(paste("~/",basename(enm_list[i]),"_SUITABILITY.png",sep=""))
  #plot(lats,latmeans,type="l",ylim=c(0,1),
  #     ylab="Suitability",xlab="Longitude")
  #lines(lats,latmeans+latsds,type="l",col="red")
  #lines(lats,latmeans-latsds,type="l",col="red")
  #dev.off()
  
## now just zoomed into the contact zone 
 enm_c3 = crop(enm2,extent(-114,-100,30,35))
  enm_c3s = enm_c3
  values(enm_c3s) = scales::rescale(values(enm_c3s),to=c(0,1))

  df = as.matrix(enm_c3s)
  latmeans=colMeans(df,na.rm=T)
  latsds = matrixStats::colSds(df,na.rm=T)

  allmean = mean(df,na.rm=T)
  allsd = sd(df,na.rm=T)
  if(is.null(enm_means_sds)){
    enm_means_sds = cbind(name=enm_list2[i],allmean,allsd)
  } else {
    enm_means_sds = rbind(enm_means_sds,cbind(name=enm_list2[i],allmean,allsd))
  }
  
  spts <- as.data.frame(rasterToPoints(enm_c3s, spatial = TRUE))
  spts = spts[,c("x","y")]
  spts = unique(spts)
  lats = unique(spts[,"x"])
  longs = unique(spts[,"y"])
  
  df2 = cbind(lats,latmeans,latsds)
  write.table(df2,paste(enm_path2,enm_list2[i],"_SUITABILITY_CONTACT_TABLE.txt",sep=""),
              sep="\t",row.names = F,quote = F)
  #lines(lats,latmeans,type="l",col=i)
 #  
 #  min1=min((latmeans-latsds),na.rm=T)
 #  max1=max((latmeans+latsds),na.rm=T)
 #  
 #  png(paste("~/",basename(enm_list[i]),"_SUITABILITY_CONTACT.png",sep=""))
 #  plot(lats,latmeans,type="l",ylim=c(0,1),
 #       ylab="Suitability",xlab="Longitude")
 #  lines(lats,latmeans+latsds,type="l",col="red")
 #  lines(lats,latmeans-latsds,type="l",col="red")
 #  dev.off()
  
} 
dev.off()
write.table(enm_means_sds,file=paste(enm_path2,"species_contact_zone_suitability_mean_sd.txt"),sep="\t",
            quote = F,row.names = F)



