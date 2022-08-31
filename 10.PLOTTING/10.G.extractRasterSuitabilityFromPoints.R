library(dismo)
library(rgdal)
library(raster)
points=read.csv("cell_locations_raster_enms.csv")
points=round(points,5)
points=unique(points)
write.csv(points,"cell_locations_raster_enms.csv")
coor = coordinates(points[c(2,1)])
coor = coordinates(points[c(1,2)])
bilineata = raster("Amphispiza bilineataAllSpecies_NotThinned_NoUncertains_WITHINPOLYGON_out_2017_10_28_00_16_23/amphispiza_maxent.asc")
fusca = raster("Melozone fuscaAllSpecies_NotThinned_NoUncertains_WITHINPOLYGON_out_2017_10_28_04_38_29/fusca_maxent.asc")
nitens = raster("Phainopepla nitensAllSpecies_NotThinned_NoUncertains_WITHINPOLYGON_out_2017_10_27_20_43_40/phainopepla_maxent.asc")
melanura = raster('Polioptila melanuraAllSpecies_NotThinned_NoUncertains_WITHINPOLYGON_out_2017_10_27_20_18_09/polioptila_maxent.asc')
bellii = raster("Vireo belliiAllSpecies_NotThinned_NoUncertains_WITHINPOLYGON_out_2017_10_28_06_42_26/vireo_maxent.asc")
sinuatus = raster("Cardinalis sinuatusAllSpecies_NotThinned_NoUncertains_WITHINPOLYGON_out_2017_10_27_20_46_31/sinuatus_maxent.asc")
flaviceps = raster("Auriparus flavicepsAllSpecies_NotThinned_NoUncertains_WITHINPOLYGON_out_2017_10_28_06_42_26/flaviceps_maxent.asc")
brunneicapillus = raster("Campylorhynchus brunneicapillusAllSpecies_NotThinned_NoUncertains_WITHINPOLYGON_out_2017_10_28_06_42_26/brunneicapillus_maxent_threshval_0.12123203_.asc")
curvirostre = raster("Toxostoma curvirostreAllSpecies_NotThinned_NoUncertains_WITHINPOLYGON_out_2017_10_28_03_11_12/curvirostre_maxent_threshval_0.09289573_.asc")
plot(fusca,xlim=c(-115,-100),ylim=c(29,34))
## extract requires LONG LAT, not LAT LONG
bil_suit = data.frame(coor,extract(bilineata,coor,small=T)); names(bil_suit) = c("LONG","LAT","SUIT_BIL")
fus_suit = data.frame(coor,extract(fusca,coor,small=T)); names(fus_suit) = c("LONG","LAT","SUIT_fus")
nit_suit = data.frame(coor,extract(nitens,coor,small=T)); names(nit_suit) = c("LONG","LAT","SUIT_NIT")
mel_suit = data.frame(coor,extract(melanura,coor,small=T)); names(mel_suit) = c("LONG","LAT","SUIT_MEL")
bel_suit = data.frame(coor,extract(bellii,coor,small=T)); names(bel_suit) = c("LONG","LAT","SUIT_BEL")
sin_suit = data.frame(coor,extract(sinuatus,coor,small=T)); names(sin_suit) = c("LONG","LAT","SUIT_SIN")
fla_suit = data.frame(coor,extract(flaviceps,coor,small=T)); names(fla_suit) = c("LONG","LAT","SUIT_fla")
bru_suit = data.frame(coor,extract(brunneicapillus,coor,small=T)); names(bru_suit) = c("LONG","LAT","SUIT_bru")
cur_suit = data.frame(coor,extract(curvirostre,coor,small=T)); names(cur_suit) = c("LONG","LAT","SUIT_cur")
bil_suit = bil_suit[order(bil_suit$LONG),]
fus_suit = fus_suit[order(fus_suit$LONG),]
nit_suit = nit_suit[order(nit_suit$LONG),]
mel_suit = mel_suit[order(mel_suit$LONG),]
bel_suit = bel_suit[order(bel_suit$LONG),]
sin_suit = sin_suit[order(sin_suit$LONG),]
fla_suit = fla_suit[order(fla_suit$LONG),]
bru_suit = bru_suit[order(bru_suit$LONG),]
cur_suit = cur_suit[order(cur_suit$LONG),]
write.table(bil_suit,"bilineata_suitability.csv")
write.table(bel_suit,"bellii_suitability.csv")
write.table(mel_suit,"melanura_suitability.csv")
write.table(sin_suit,"sinuatus_suitability.csv")
write.table(fus_suit,"fusca_suitability.csv")
write.table(fla_suit,"flaviceps_suitability.csv")
write.table(nit_suit,"nitens_suitability.csv")
write.table(bru_suit,"brunneicapillus_suitability.csv")
write.table(cur_suit,"curvirostre_suitability.csv")
bil_suit = read.table("bilineata_suitability.csv")
bel_suit = read.table("bellii_suitability.csv")
mel_suit = read.table("melanura_suitability.csv")
sin_suit = read.table("sinuatus_suitability.csv")
fus_suit = read.table("fusca_suitability.csv")
fla_suit = read.table("flaviceps_suitability.csv")
nit_suit = read.table("nitens_suitability.csv")
bru_suit = read.table("brunneicapillus_suitability.csv")
cur_suit = read.table("curvirostre_suitability.csv")
## GET SUITS 
spp = c("flaviceps","bellii","curvirostre","melanura","brunneicapillus","fusca","bilineata","sinuatus","nitens")
thresh = c("0.01578953","0.14259912","0.09289573","0.02024653",
           "0.12123203","0.09119378","0.19570312","0.05652053","0.18654452")
bil_barr = bil_suit[bil_suit$LONG<=-108,]; bil_barr = bil_barr[bil_barr$LONG>=-112,]
bel_barr = bel_suit[bel_suit$LONG<=-108,]; bel_barr = bel_barr[bel_barr$LONG>=-112,]
mel_barr = mel_suit[mel_suit$LONG<=-108,]; mel_barr = mel_barr[mel_barr$LONG>=-112,]
sin_barr = sin_suit[sin_suit$LONG<=-108,]; sin_barr = sin_barr[sin_barr$LONG>=-112,]
fla_barr = fla_suit[fla_suit$LONG<=-108,]; fla_barr = fla_barr[fla_barr$LONG>=-112,]
fus_barr = fus_suit[fus_suit$LONG<=-108,]; fus_barr = fus_barr[fus_barr$LONG>=-112,]
nit_barr = nit_suit[nit_suit$LONG<=-108,]; nit_barr = nit_barr[nit_barr$LONG>=-112,]
bru_barr = bru_suit[bru_suit$LONG<=-108,]; bru_barr = bru_barr[bru_barr$LONG>=-112,]
cur_barr = cur_suit[cur_suit$LONG<=-108,]; cur_barr = cur_barr[cur_barr$LONG>=-112,]
sum(bil_barr$SUIT_BIL>=0.2)
sum(bil_barr$SUIT_BIL<0.2) 
bil_prop20 = sum(bil_barr$SUIT_BIL<0.2) / sum(bil_barr$SUIT_BIL>0)
bel_prop20 = sum(bel_barr$SUIT_BEL<0.2) / sum(bel_barr$SUIT_BEL>0)
mel_prop20 = sum(mel_barr$SUIT_MEL<0.2) / sum(mel_barr$SUIT_MEL>0)
sin_prop20 = sum(sin_barr$SUIT_SIN<0.2) / sum(sin_barr$SUIT_SIN>0)
fus_prop20 = sum(fus_barr$SUIT_fus<0.2) / sum(fus_barr$SUIT_fus>0)
fla_prop20 = sum(fla_barr$SUIT_fla<0.2) / sum(fla_barr$SUIT_fla>0)
nit_prop20 = sum(nit_barr$SUIT_NIT<0.2) / sum(nit_barr$SUIT_NIT>0)
bru_prop20 = sum(bru_barr$SUIT_bru<0.2) / sum(bru_barr$SUIT_bru>0)
cur_prop20 = sum(cur_barr$SUIT_cur<0.2) / sum(cur_barr$SUIT_cur>0)
bil_prop20 #0.03841932
bel_prop20 #0.2755214
mel_prop20 #0.2184413
sin_prop20 #0.07025247
fus_prop20 #0.004390779
fla_prop20 #0.06586169
nit_prop20 #0.07903403
bru_prop20 #0.04390779
cur_prop20 #0.03512623
## thresh
"bellii          0.14259912
bilineata       0.19570312
brunneicapillus 0.12123203
curvirostre     0.09289573       
flaviceps       0.01578953
fusca           0.09119378
melanura        0.02024653
nitens          0.18654452
sinuatus        0.05652053"
bil_propThresh = sum(bil_barr$SUIT_BIL<0.19570312) / sum(bil_barr$SUIT_BIL>0)
bel_propThresh = sum(bel_barr$SUIT_BEL<0.14259912) / sum(bel_barr$SUIT_BEL>0)
mel_propThresh = sum(mel_barr$SUIT_MEL<0.02024653) / sum(mel_barr$SUIT_MEL>0)
sin_propThresh = sum(sin_barr$SUIT_SIN<0.05652053) / sum(sin_barr$SUIT_SIN>0)
fus_propThresh = sum(fus_barr$SUIT_fus<0.09119378) / sum(fus_barr$SUIT_fus>0)
fla_propThresh = sum(fla_barr$SUIT_fla<0.01578953) / sum(fla_barr$SUIT_fla>0)
nit_propThresh = sum(nit_barr$SUIT_NIT<0.18654452) / sum(nit_barr$SUIT_NIT>0)
bru_propThresh = sum(bru_barr$SUIT_bru<0.12123203) / sum(bru_barr$SUIT_bru>0)
cur_propThresh = sum(cur_barr$SUIT_cur<0.09289573) / sum(cur_barr$SUIT_cur>0)
bil_propThresh #0.03841932
bel_propThresh #0.219539
mel_propThresh #0.02744237
sin_propThresh #0.04500549
fus_propThresh #0
fla_propThresh #0.008781559
nit_propThresh #0.07903403
bru_propThresh #0.03732162
cur_propThresh #0.03293085
maxsuit = max(c(bil_suit$SUIT_BIL,
                nit_suit$SUIT_fus,
                nit_suit$SUIT_NIT,
                mel_suit$SUIT_MEL,
                bel_suit$SUIT_BEL,
                sin_suit$SUIT_SIN))
minsuit = min(c(bil_suit$SUIT_BIL,
                nit_suit$SUIT_fus,
                nit_suit$SUIT_NIT,
                mel_suit$SUIT_MEL,
                bel_suit$SUIT_BEL,
                sin_suit$SUIT_SIN))
plot(bil_suit$LONG,bil_suit$SUIT_BIL,xlim=c(-115,-100),ylim=c(0,1),
     col="black",pch=1,type="n")
points(bel_suit$LONG,bel_suit$SUIT_BEL,xlim=c(-115,-100),ylim=c(0,1),col="magenta",pch=5)
points(bil_suit$LONG,bil_suit$SUIT_BIL,xlim=c(-115,-100),ylim=c(0,1),col="black",pch=1)
points(fus_suit$LONG,fus_suit$SUIT_fus,xlim=c(-115,-100),ylim=c(0,1),col="blue",pch=3)
points(mel_suit$LONG,mel_suit$SUIT_MEL,xlim=c(-115,-100),ylim=c(0,1),col="green",pch=4)
points(nit_suit$LONG,nit_suit$SUIT_NIT,xlim=c(-115,-100),ylim=c(0,1),col="red",pch=2)
points(sin_suit$LONG,sin_suit$SUIT_SIN,xlim=c(-115,-100),ylim=c(0,1),col="orange",pch=6)
smooth_bel = smooth.spline(bel_suit$LONG, bel_suit$SUIT_BEL, spar=0.5); lines(smooth_bel,xlim=c(-115,-100), ylim=c(0,1),col="magenta")
smooth_bil = smooth.spline(bil_suit$LONG, bil_suit$SUIT_BIL, spar=0.5); lines(smooth_bil,xlim=c(-115,-100), ylim=c(0,1),col="black")
smooth_fus = smooth.spline(fus_suit$LONG, fus_suit$SUIT_FUS, spar=0.5); lines(smooth_fus,xlim=c(-115,-100), ylim=c(0,1),col="blue")
smooth_mel = smooth.spline(mel_suit$LONG, mel_suit$SUIT_MEL, spar=0.5); lines(smooth_mel,xlim=c(-115,-100), ylim=c(0,1),col="green")
smooth_nit = smooth.spline(nit_suit$LONG, nit_suit$SUIT_NIT, spar=0.5); lines(smooth_nit,xlim=c(-115,-100), ylim=c(0,1),col="red")
smooth_sin = smooth.spline(sin_suit$LONG, sin_suit$SUIT_SIN, spar=0.5); lines(smooth_sin,xlim=c(-115,-100), ylim=c(0,1),col="blue")
legend("topleft",col=c("black","blue","red","green","orange","magenta"),
       legend=c("BIL","FUS","NIT","MEL","SIN","BEL"),
       pch=c(1,3,2,4,6,5))
png("nicheModelSuitability_fullRange_lines.png",
    width=12,height=8,pointsize=8,units="cm",
    res=500)
par(mfrow=c(2,4),bg="black",col.axis="white",col.lab="white",
    col.main="white",col.sub="white",
    lwd=1,font.main=4,mar=c(2,2,2,1))
{
  plot(bel_suit$LONG,bel_suit$SUIT_BEL,xlim=c(-115,-100),ylim=c(0,1),col=adjustcolor("white", alpha=0.2),pch=1,main="Vireo bellii",ylab="",xlab="")
  smooth_bel = smooth.spline(bel_suit$LONG, bel_suit$SUIT_BEL, spar=1); lines(smooth_bel,xlim=c(-115,-100), ylim=c(0,1),col="red")
  axis(1,col="white",col.ticks="white",lwd=2); axis(2,col="white",col.ticks="white",lwd=2)
  abline(v=-108,lty=2,col="cyan")
  abline(v=-112,lty=2,col="cyan")
  plot(fla_suit$LONG,fla_suit$SUIT_fla,xlim=c(-115,-100),ylim=c(0,1),col=adjustcolor("white", alpha=0.2),pch=1,main="Auriparus flaviceps",ylab="",xlab="")
  smooth_bel = smooth.spline(fla_suit$LONG, fla_suit$SUIT_fla, spar=1); lines(smooth_bel,xlim=c(-115,-100), ylim=c(0,1),col="red")
  axis(1,col="white",col.ticks="white",lwd=2); axis(2,col="white",col.ticks="white",lwd=2)
  abline(v=-108,lty=2,col="cyan")
  abline(v=-112,lty=2,col="cyan")
  plot(bil_suit$LONG,bil_suit$SUIT_BIL,xlim=c(-115,-100),ylim=c(0,1),col=adjustcolor("white", alpha=0.2),pch=1,main="Amphispiza bilineata",ylab="",xlab="")
  smooth_bil = smooth.spline(bil_suit$LONG, bil_suit$SUIT_BIL, spar=1); lines(smooth_bil,xlim=c(-115,-100), ylim=c(0,1),col="red")
  axis(1,col="white",col.ticks="white",lwd=2); axis(2,col="white",col.ticks="white",lwd=2)
  abline(v=-108,lty=2,col="cyan")
  abline(v=-112,lty=2,col="cyan")
  plot(fus_suit$LONG,fus_suit$SUIT_fus,xlim=c(-115,-100),ylim=c(0,1),col=adjustcolor("white", alpha=0.2),pch=1,main="Melozone fusca",ylab="",xlab="")
  smooth_fus = smooth.spline(fus_suit$LONG, fus_suit$SUIT_fus, spar=1); lines(smooth_fus,xlim=c(-115,-100), ylim=c(0,1),col="red")
  axis(1,col="white",col.ticks="white",lwd=2); axis(2,col="white",col.ticks="white",lwd=2)
  abline(v=-108,lty=2,col="cyan")
  abline(v=-112,lty=2,col="cyan")
  plot(mel_suit$LONG,mel_suit$SUIT_MEL,xlim=c(-115,-100),ylim=c(0,1),col=adjustcolor("white", alpha=0.2),pch=1,main="Polioptila melanura",ylab="",xlab="")
  smooth_mel = smooth.spline(mel_suit$LONG, mel_suit$SUIT_MEL, spar=1); lines(smooth_mel,xlim=c(-115,-100), ylim=c(0,1),col="red")
  axis(1,col="white",col.ticks="white",lwd=2); axis(2,col="white",col.ticks="white",lwd=2)
  abline(v=-108,lty=2,col="cyan")
  abline(v=-112,lty=2,col="cyan")
  plot(nit_suit$LONG,nit_suit$SUIT_NIT,xlim=c(-115,-100),ylim=c(0,1),col=adjustcolor("white", alpha=0.2),pch=1,main="Phainopepla nitens",ylab="",xlab="")
  smooth_nit = smooth.spline(nit_suit$LONG, nit_suit$SUIT_NIT, spar=1); lines(smooth_nit,xlim=c(-115,-100), ylim=c(0,1),col="red")
  axis(1,col="white",col.ticks="white",lwd=2); axis(2,col="white",col.ticks="white",lwd=2)
  abline(v=-108,lty=2,col="cyan")
  abline(v=-112,lty=2,col="cyan")
  plot(sin_suit$LONG,sin_suit$SUIT_SIN,xlim=c(-115,-100),ylim=c(0,1),col=adjustcolor("white", alpha=0.2),pch=1,main="Cardinalis sinuatus",ylab="",xlab="")
  smooth_sin = smooth.spline(sin_suit$LONG, sin_suit$SUIT_SIN, spar=1); lines(smooth_sin,xlim=c(-115,-100), ylim=c(0,1),col="red")
  axis(1,col="white",col.ticks="white",lwd=2); axis(2,col="white",col.ticks="white",lwd=2)
  abline(v=-108,lty=2,col="cyan")
  abline(v=-112,lty=2,col="cyan")
  plot(bru_suit$LONG,bru_suit$SUIT_bru,xlim=c(-115,-100),ylim=c(0,1),col=adjustcolor("white", alpha=0.2),pch=1,main="Cardinalis bruuatus",ylab="",xlab="")
  smooth_bru = smooth.spline(bru_suit$LONG, bru_suit$SUIT_bru, spar=1); lines(smooth_bru,xlim=c(-115,-100), ylim=c(0,1),col="red")
  axis(1,col="white",col.ticks="white",lwd=2); axis(2,col="white",col.ticks="white",lwd=2)
  abline(v=-108,lty=2,col="cyan")
  abline(v=-112,lty=2,col="cyan")
}
dev.off()
plot(bil_suit$LAT,bil_suit$SUIT_BIL,xlim=c(29,35),ylim=c(0,1),
     col="black",pch=1,type="n")
points(bel_suit$LAT,bel_suit$SUIT_BEL,xlim=c(29,35),ylim=c(0,1),col="magenta",pch=5)
points(bil_suit$LAT,bil_suit$SUIT_BIL,xlim=c(29,35),ylim=c(0,1),col="black",pch=1)
points(fus_suit$LAT,fus_suit$SUIT_FUS,xlim=c(29,35),ylim=c(0,1),col="blue",pch=3)
points(mel_suit$LAT,mel_suit$SUIT_MEL,xlim=c(29,35),ylim=c(0,1),col="green",pch=4)
points(nit_suit$LAT,nit_suit$SUIT_NIT,xlim=c(29,35),ylim=c(0,1),col="red",pch=2)
points(sin_suit$LAT,sin_suit$SUIT_SIN,xlim=c(29,35),ylim=c(0,1),col="orange",pch=6)
smooth_bel = smooth.spline(bel_suit$LAT, bel_suit$SUIT_BEL, spar=1); lines(smooth_bel,xlim=c(29,35), ylim=c(0,1),col="magenta")
smooth_bil = smooth.spline(bil_suit$LAT, bil_suit$SUIT_BIL, spar=1); lines(smooth_bil,xlim=c(29,35), ylim=c(0,1),col="black")
smooth_fus = smooth.spline(fus_suit$LAT, fus_suit$SUIT_FUS, spar=1); lines(smooth_fus,xlim=c(29,35), ylim=c(0,1),col="blue")
smooth_mel = smooth.spline(mel_suit$LAT, mel_suit$SUIT_MEL, spar=1); lines(smooth_mel,xlim=c(29,35), ylim=c(0,1),col="green")
smooth_nit = smooth.spline(nit_suit$LAT, nit_suit$SUIT_NIT, spar=1); lines(smooth_nit,xlim=c(29,35), ylim=c(0,1),col="red")
smooth_sin = smooth.spline(sin_suit$LAT, sin_suit$SUIT_SIN, spar=1); lines(smooth_sin,xlim=c(29,35), ylim=c(0,1),col="blue")
plot(bil_suit$LONG,bil_suit$SUIT_BIL)
plot(nit_suit$LONG,nit_suit$SUIT_nit)
plot(nit_suit$LONG,nit_suit$SUIT_NIT)
plot(mel_suit$LONG,mel_suit$SUIT_MEL)
plot(bel_suit$LONG,bel_suit$SUIT_BEL)
plot(sin_suit$LONG,sin_suit$SUIT_SIN)
plot(bil_suit$LONG, bil_suit$SUIT_BIL, xlim=range(bil_suit$LONG), ylim=range(bil_suit$SUIT_BIL), xlab="x", ylab="y", 
     main = "noise-less data",pch=16)
lines(bil_suit$LONG[order(bil_suit$LONG)], bil_suit$SUIT_BIL[order(bil_suit$LONG)], xlim=range(bil_suit$LONG), ylim=range(bil_suit$SUIT_BIL), pch=16)
)
##########
## aggregated data for each species
bil_suit$LONG_CHAR = as.character(bil_suit$LONG)
bil_med = aggregate(bil_suit$SUIT_BIL,list(bil_suit$LONG_CHAR),quantile,probs=c(0.05,0.25,0.50,0.75,0.95))
colnames(bil_med$x) = c("q5","q25","q50","q75","q95")
plot(bil_med$Group.1,bil_med$x[,3],ylim=c(0,1))
lines(bil_med$Group.1,bil_med$x[,4],col="red")
lines(bil_med$Group.1,bil_med$x[,2],col="blue")
lines(bil_med$Group.1,bil_med$x[,5],col="pink")
lines(bil_med$Group.1,bil_med$x[,1],col="lightblue")
fus_suit$LONG_CHAR = as.character(fus_suit$LONG)
fus_med = aggregate(fus_suit$SUIT_fus,list(fus_suit$LONG_CHAR),quantile,probs=c(0.05,0.25,0.50,0.75,0.95))
colnames(fus_med$x) = c("q5","q25","q50","q75","q95")
plot(fus_med$Group.1,fus_med$x[,3],ylim=c(0,1))
lines(fus_med$Group.1,fus_med$x[,4],col="red")
lines(fus_med$Group.1,fus_med$x[,2],col="blue")
lines(fus_med$Group.1,fus_med$x[,5],col="pink")
lines(fus_med$Group.1,fus_med$x[,1],col="lightblue")
###
nit_suit$LONG_CHAR = as.character(nit_suit$LONG)
nit_med = aggregate(nit_suit$SUIT_NIT,list(nit_suit$LONG_CHAR),quantile,probs=c(0.05,0.25,0.50,0.75,0.95))
colnames(nit_med$x) = c("q5","q25","q50","q75","q95")
plot(nit_med$Group.1,nit_med$x[,3],ylim=c(0,1))
lines(nit_med$Group.1,nit_med$x[,4],col="red")
lines(nit_med$Group.1,nit_med$x[,2],col="blue")
lines(nit_med$Group.1,nit_med$x[,5],col="pink")
lines(nit_med$Group.1,nit_med$x[,1],col="lightblue")
####
mel_suit$LONG_CHAR = as.character(mel_suit$LONG)
mel_med = aggregate(mel_suit$SUIT_MEL,list(mel_suit$LONG_CHAR),quantile,probs=c(0.05,0.25,0.50,0.75,0.95))
colnames(mel_med$x) = c("q5","q25","q50","q75","q95")
plot(mel_med$Group.1,mel_med$x[,3],ylim=c(0,1))
lines(mel_med$Group.1,mel_med$x[,4],col="red")
lines(mel_med$Group.1,mel_med$x[,2],col="blue")
lines(mel_med$Group.1,mel_med$x[,5],col="pink")
lines(mel_med$Group.1,mel_med$x[,1],col="lightblue")
#####
bel_suit$LONG_CHAR = as.character(bel_suit$LONG)
bel_med = aggregate(bel_suit$SUIT_BEL,list(bel_suit$LONG_CHAR),quantile,probs=c(0.05,0.25,0.50,0.75,0.95))
colnames(bel_med$x) = c("q5","q25","q50","q75","q95")
plot(bel_med$Group.1,bel_med$x[,3],ylim=c(0,1))
lines(bel_med$Group.1,bel_med$x[,4],col="red")
lines(bel_med$Group.1,bel_med$x[,2],col="blue")
lines(bel_med$Group.1,bel_med$x[,5],col="pink")
lines(bel_med$Group.1,bel_med$x[,1],col="lightblue")
#####
sin_suit$LONG_CHAR = as.character(sin_suit$LONG)
sin_med = aggregate(sin_suit$SUIT_SIN,list(sin_suit$LONG_CHAR),quantile,probs=c(0.05,0.25,0.50,0.75,0.95))
colnames(sin_med$x) = c("q5","q25","q50","q75","q95")
plot(sin_med$Group.1,sin_med$x[,3],ylim=c(0,1))
lines(sin_med$Group.1,sin_med$x[,4],col="red")
lines(sin_med$Group.1,sin_med$x[,2],col="blue")
lines(sin_med$Group.1,sin_med$x[,5],col="pink")
lines(sin_med$Group.1,sin_med$x[,1],col="lightblue")
###
fla_suit$LONG_CHAR = as.character(fla_suit$LONG)
fla_med = aggregate(fla_suit$SUIT_fla,list(fla_suit$LONG_CHAR),quantile,probs=c(0.05,0.25,0.50,0.75,0.95))
colnames(fla_med$x) = c("q5","q25","q50","q75","q95")
plot(fla_med$Group.1,fla_med$x[,3],ylim=c(0,1))
lines(fla_med$Group.1,fla_med$x[,4],col="red")
lines(fla_med$Group.1,fla_med$x[,2],col="blue")
lines(fla_med$Group.1,fla_med$x[,5],col="pink")
lines(fla_med$Group.1,fla_med$x[,1],col="lightblue")
bru_suit$LONG_CHAR = as.character(bru_suit$LONG)
bru_med = aggregate(bru_suit$SUIT_bru,list(bru_suit$LONG_CHAR),quantile,probs=c(0.05,0.25,0.50,0.75,0.95))
colnames(bru_med$x) = c("q5","q25","q50","q75","q95")
plot(bru_med$Group.1,bru_med$x[,3],ylim=c(0,1))
lines(bru_med$Group.1,bru_med$x[,4],col="red")
lines(bru_med$Group.1,bru_med$x[,2],col="blue")
lines(bru_med$Group.1,bru_med$x[,5],col="pink")
lines(bru_med$Group.1,bru_med$x[,1],col="lightblue")
##
plot(bel_med$Group.1,bel_med$x[,3],ylim=c(0,1),type="n")
points(bel_med$Group.1,bel_med$x[,3],ylim=c(0,1),col="red",pch=1)
points(bil_med$Group.1,bil_med$x[,3],ylim=c(0,1),col="orange",pch=2)
points(fla_med$Group.1,fla_med$x[,3],ylim=c(0,1),col="green",pch=3)
points(fus_med$Group.1,fus_med$x[,3],ylim=c(0,1),col="blue",pch=4)
points(mel_med$Group.1,mel_med$x[,3],ylim=c(0,1),col="purple",pch=5)
points(nit_med$Group.1,nit_med$x[,3],ylim=c(0,1),col="magenta",pch=6)
points(sin_med$Group.1,sin_med$x[,3],ylim=c(0,1),col="black",pch=7)
plot(bel_med$Group.1,bel_med$x[,3],ylim=c(0,1),type="n")
lines(bel_med$Group.1,bel_med$x[,3],ylim=c(0,1),col="red",pch=1)
lines(bil_med$Group.1,bil_med$x[,3],ylim=c(0,1),col="orange",pch=2)
lines(fla_med$Group.1,fla_med$x[,3],ylim=c(0,1),col="green",pch=3)
lines(fus_med$Group.1,fus_med$x[,3],ylim=c(0,1),col="blue",pch=4)
lines(mel_med$Group.1,mel_med$x[,3],ylim=c(0,1),col="purple",pch=5)
lines(nit_med$Group.1,nit_med$x[,3],ylim=c(0,1),col="magenta",pch=6)
lines(sin_med$Group.1,sin_med$x[,3],ylim=c(0,1),col="black",pch=7)
####
library(ENMTools)
raster.breadth(bilineata) ## b1, 0.8533968. b2, 0.05962488
raster.breadth(nitens) # 0.839922, 0.04222851
raster.breadth(fusca) # 0.8515305, 0.05971484
raster.breadth(melanura) # 0.8389893, 0.04339477
raster.breadth(sinuatus) # 0.8458246, 0.05258115
raster.breadth(flaviceps) # 0.848704, 0.05331433
raster.breadth(bellii) # 0.8465935, 0.04825829
breadth = read.table(text="bilineata	0.8533968	0.05962488
nitens	0.839922	0.04222851
fusca	0.8515305	0.05971484
melanura	0.8389893	0.04339477
sinuatus	0.8458246	0.05258115
flaviceps	0.848704	0.05331433
bellii	0.8465935	0.04825829")
colnames(breadth) = c("spp","b1","b2")
plot(breadth$b1,breadth$b2)
barplot(breadth$b1)
barplot(breadth$b2)
####
png("nicheModelSuitability_fullRange_median_lines.png",
    width=16,height=8,pointsize=8,units="cm",
    res=500)
par(mfrow=c(2,4),bg="black",col.axis="white",col.lab="white",
    col.main="white",col.sub="white",
    lwd=1,font.main=4,mar=c(2,2,2,1))
{
  plot(bel_suit$LONG,bel_suit$SUIT_BEL,xlim=c(-115,-100),ylim=c(0,1),col=adjustcolor("white", alpha=0.2),pch=1,main="Vireo bellii",ylab="",xlab=""); axis(1,col="white",col.ticks="white",lwd=2); axis(2,col="white",col.ticks="white",lwd=2)
  lines(bel_med$Group.1,bel_med$x[,3],ylim=c(0,1),col="red",pch=1)
  abline(v=-108,lty=2,col="cyan")
  abline(v=-112,lty=2,col="cyan")
  plot(bil_suit$LONG,bil_suit$SUIT_BIL,xlim=c(-115,-100),ylim=c(0,1),col=adjustcolor("white", alpha=0.2),pch=1,main="Amphispiza bilineata",ylab="",xlab=""); axis(1,col="white",col.ticks="white",lwd=2); axis(2,col="white",col.ticks="white",lwd=2)
  lines(bil_med$Group.1,bil_med$x[,3],ylim=c(0,1),col="red",pch=2)
  abline(v=-108,lty=2,col="cyan")
  abline(v=-112,lty=2,col="cyan")
  plot(fla_suit$LONG,fla_suit$SUIT_fla,xlim=c(-115,-100),ylim=c(0,1),col=adjustcolor("white", alpha=0.2),pch=1,main="Auriparus flaviceps",ylab="",xlab=""); axis(1,col="white",col.ticks="white",lwd=2); axis(2,col="white",col.ticks="white",lwd=2)
  lines(fla_med$Group.1,fla_med$x[,3],ylim=c(0,1),col="red",pch=3)
  abline(v=-108,lty=2,col="cyan")
  abline(v=-112,lty=2,col="cyan")
  plot(fus_suit$LONG,fus_suit$SUIT_fus,xlim=c(-115,-100),ylim=c(0,1),col=adjustcolor("white", alpha=0.2),pch=1,main="Melozone fusca",ylab="",xlab=""); axis(1,col="white",col.ticks="white",lwd=2); axis(2,col="white",col.ticks="white",lwd=2)
  lines(fus_med$Group.1,fus_med$x[,3],ylim=c(0,1),col="red",pch=4)
  abline(v=-108,lty=2,col="cyan")
  abline(v=-112,lty=2,col="cyan")
  plot(mel_suit$LONG,mel_suit$SUIT_MEL,xlim=c(-115,-100),ylim=c(0,1),col=adjustcolor("white", alpha=0.2),pch=1,main="Polioptila melanura",ylab="",xlab=""); axis(1,col="white",col.ticks="white",lwd=2); axis(2,col="white",col.ticks="white",lwd=2)
  lines(mel_med$Group.1,mel_med$x[,3],ylim=c(0,1),col="red",pch=5)
  abline(v=-108,lty=2,col="cyan")
  abline(v=-112,lty=2,col="cyan")
  plot(nit_suit$LONG,nit_suit$SUIT_NIT,xlim=c(-115,-100),ylim=c(0,1),col=adjustcolor("white", alpha=0.2),pch=1,main="Phainopepla nitens",ylab="",xlab=""); axis(1,col="white",col.ticks="white",lwd=2); axis(2,col="white",col.ticks="white",lwd=2)
  lines(nit_med$Group.1,nit_med$x[,3],ylim=c(0,1),col="red",pch=6)
  abline(v=-108,lty=2,col="cyan")
  abline(v=-112,lty=2,col="cyan")
  plot(sin_suit$LONG,sin_suit$SUIT_SIN,xlim=c(-115,-100),ylim=c(0,1),col=adjustcolor("white", alpha=0.2),pch=1,main="Cardinalis sinuatus",ylab="",xlab=""); axis(1,col="white",col.ticks="white",lwd=2); axis(2,col="white",col.ticks="white",lwd=2)
  lines(sin_med$Group.1,sin_med$x[,3],ylim=c(0,1),col="red",pch=7)
  abline(v=-108,lty=2,col="cyan")
  abline(v=-112,lty=2,col="cyan")
  plot(bru_suit$LONG,bru_suit$SUIT_bru,xlim=c(-115,-100),ylim=c(0,1),col=adjustcolor("white", alpha=0.2),pch=1,main="Campylorhynchus brunneicapillus",ylab="",xlab=""); axis(1,col="white",col.ticks="white",lwd=2); axis(2,col="white",col.ticks="white",lwd=2)
  lines(bru_med$Group.1,bru_med$x[,3],ylim=c(0,1),col="red",pch=7)
  abline(v=-108,lty=2,col="cyan")
  abline(v=-112,lty=2,col="cyan")
}
dev.off()
png("nicheModelSuitability_fullRange_median_vir_fus_lines.png",
    width=4,height=8,pointsize=8,units="cm",
    res=500)
par(mfrow=c(2,1),bg="black",col.axis="white",col.lab="white",
    col.main="white",col.sub="white",
    lwd=1,font.main=4,mar=c(2,2,2,1))
{
  plot(bel_suit$LONG,bel_suit$SUIT_BEL,xlim=c(-115,-100),ylim=c(0,1),col=adjustcolor("white", alpha=0.2),pch=1,main="Vireo bellii",ylab="",xlab=""); axis(1,col="white",col.ticks="white",lwd=2); axis(2,col="white",col.ticks="white",lwd=2)
  lines(bel_med$Group.1,bel_med$x[,3],ylim=c(0,1),col="red",pch=1)
  abline(v=-108,lty=2,col="cyan")
  abline(v=-112,lty=2,col="cyan")
  plot(fus_suit$LONG,fus_suit$SUIT_fus,xlim=c(-115,-100),ylim=c(0,1),col=adjustcolor("white", alpha=0.2),pch=1,main="Melozone fusca",ylab="",xlab=""); axis(1,col="white",col.ticks="white",lwd=2); axis(2,col="white",col.ticks="white",lwd=2)
  lines(fus_med$Group.1,fus_med$x[,3],ylim=c(0,1),col="red",pch=4)
  abline(v=-108,lty=2,col="cyan")
  abline(v=-112,lty=2,col="cyan")
}
dev.off()
png("nicheModelSuitability_fullRange_hists_vir_fus.png",
    width=4,height=8,pointsize=8,units="cm",
    res=500)
par(mfrow=c(2,1),bg="black",col.axis="white",col.lab="white",
    col.main="white",col.sub="white",
    lwd=1,font.main=4,mar=c(2,2,2,1))
{
  hist(bel_suit$SUIT_BEL,xlim=c(0,1),border="white",pch=1,breaks=seq(0,1,1/30),main="Vireo bellii",ylab="",xlab="")
  axis(1,col="white",col.ticks="white",lwd=2); axis(2,col="white",col.ticks="white",lwd=2)
  hist(fus_suit$SUIT_fus,xlim=c(0,1),border="white",pch=1,breaks=seq(0,1,1/30),main="Melozone fusca",ylab="",xlab="")
  axis(1,col="white",col.ticks="white",lwd=2); axis(2,col="white",col.ticks="white",lwd=2)
}
dev.off()
png("nicheModelSuitability_fullRange_hists.png",
    width=16,height=8,pointsize=8,units="cm",
    res=500)
par(mfrow=c(2,4),bg="black",col.axis="white",col.lab="white",
    col.main="white",col.sub="white",
    lwd=1,font.main=4,mar=c(2,2,2,1))
{
  hist(bel_suit$SUIT_BEL,xlim=c(0,1),border="white",pch=1,breaks=seq(0,1,1/30),main="Vireo bellii",ylab="",xlab="")
  axis(1,col="white",col.ticks="white",lwd=2); axis(2,col="white",col.ticks="white",lwd=2,ylim=c(0,400))
  hist(bil_suit$SUIT_BIL,xlim=c(0,1),border="white",pch=1,breaks=seq(0,1,1/30),main="Amphispiza bilineata",ylab="",xlab="")
  axis(1,col="white",col.ticks="white",lwd=2); axis(2,col="white",col.ticks="white",lwd=2,ylim=c(0,400))
  hist(fla_suit$SUIT_fla,xlim=c(0,1),border="white",pch=1,breaks=seq(0,1,1/30),main="Auriparus flaviceps",ylab="",xlab="")
  axis(1,col="white",col.ticks="white",lwd=2); axis(2,col="white",col.ticks="white",lwd=2,ylim=c(0,400))
  hist(fus_suit$SUIT_fus,xlim=c(0,1),border="white",pch=1,breaks=seq(0,1,1/30),main="Melozone fusca",ylab="",xlab="")
  axis(1,col="white",col.ticks="white",lwd=2); axis(2,col="white",col.ticks="white",lwd=2,ylim=c(0,400))
  hist(mel_suit$SUIT_MEL,xlim=c(0,1),border="white",pch=1,breaks=seq(0,1,1/30),main="Polioptila melanura",ylab="",xlab="")
  axis(1,col="white",col.ticks="white",lwd=2); axis(2,col="white",col.ticks="white",lwd=2,ylim=c(0,400))
  hist(nit_suit$SUIT_NIT,xlim=c(0,1),border="white",pch=1,breaks=seq(0,1,1/30),main="Phainopepla nitens",ylab="",xlab="")
  axis(1,col="white",col.ticks="white",lwd=2); axis(2,col="white",col.ticks="white",lwd=2,ylim=c(0,400))
  hist(sin_suit$SUIT_SIN,xlim=c(0,1),border="white",pch=1,breaks=seq(0,1,1/30),main="Cardinalis sinuatus",ylab="",xlab="")
  axis(1,col="white",col.ticks="white",lwd=2); axis(2,col="white",col.ticks="white",lwd=2,ylim=c(0,400))
  hist(bru_suit$SUIT_bru,xlim=c(0,1),border="white",pch=1,breaks=seq(0,1,1/30),main="Campylorhynchus brunneicapillus",ylab="",xlab="")
  axis(1,col="white",col.ticks="white",lwd=2); axis(2,col="white",col.ticks="white",lwd=2,ylim=c(0,400))
}
dev.off()
