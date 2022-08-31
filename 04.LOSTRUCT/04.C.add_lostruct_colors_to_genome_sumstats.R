sumstats_file = read.table("rec_taj_dxy_fst_islswp_miss_lostruct.temp",
                           header=T,fill=T)
aggregate(sumstats_file$Fst~sumstats_file$color+sumstats_file$species,FUN=function(x){mean(x,na.rm=T)})
bel=read.table("Vireo.bellii.called.geno.PseudoNC.all.sorted.nospace.coords_SMALL.csv",
               header=T,fill=T,sep=",",comment.char="")
bil=read.table("Amphispiza.bilineata.called.geno.PseudoNC.all.sorted.nospace.coords_SMALL.csv",
               header=T,fill=T,sep=",",comment.char="")
fla=read.table("Auriparus.flaviceps.called.geno.PseudoNC.all.sorted.nospace.coords_SMALL.csv",
               header=T,fill=T,sep=",",comment.char="")
bru=read.table("Campylorhynchus.brunneicapillus.called.geno.PseudoNC.all.fixedchroms.converted.sorted.nospace.coords_SMALL.csv",
               header=T,fill=T,sep=",",comment.char="")
cri=read.table("Toxostoma.crissale.called.geno.PseudoNC.all.fixedchroms.converted.sorted.nospace.coords_SMALL.csv",
               header=T,fill=T,sep=",",comment.char="")
cur=read.table("Toxostoma.curvirostre.called.geno.PseudoNC.all.fixedchroms.converted.sorted.nospace.coords_SMALL.csv",
               header=T,fill=T,sep=",",comment.char="")
fus=read.table("Melozone.fusca_fixed.sorted.nospace.coords_SMALL.csv",
               header=T,fill=T,sep=",",comment.char="")
mel=read.table("Polioptila.melanura.called.geno.PseudoNC.all.fixedchroms.converted.sorted.nospace.coords_SMALL.csv",
               header=T,fill=T,sep=",",comment.char="")
nit=read.table("Phainopepla.nitens.called.geno.PseudoNC.all.fixedchroms.converted.sorted.nospace.coords_SMALL.csv",
               header=T,fill=T,sep=",",comment.char="")
sin=read.table("Cardinalis.sinuatus.called.geno.PseudoNC.all.fixedchroms.converted.sorted.nospace.coords_SMALL.csv",
               header=T,fill=T,sep=",",comment.char="")
bil$species = "bil"
bel$species = "bel"
bru$species = "bru"
fla$species = "fla"
fus$species = "fus"
cri$species = "cri"
cur$species = "cur"
mel$species = "mel"
nit$species = "nit"
sin$species = "sin"
colors = rbind(bel,bil,bru,cri,cur,fla,fus,mel,nit,sin)
colors=unique(colors)
for(spp in c("bel","bil","bru","cri","cur","fla","fus","mel","nit","sin")){
  color_sub = colors[colors$species==spp,]
  rows=nrow(color_sub)
  for(i in 1:rows) {
    if(i %% 50 == 0) {print(paste(i,spp,rows))}
    this_start = color_sub$start[i]
    this_stop = color_sub$end[i]
    this_chrom = color_sub$chrom[i]
    this_species = color_sub$species[i]
    this_col = color_sub$ccols[i]
    sumstats_file$color[sumstats_file$species==this_species & sumstats_file$chr==this_chrom & sumstats_file$midPos>= this_start & sumstats_file$midPos<= this_stop & 
                          sumstats_file$windowstarts<=this_stop & sumstats_file$windowstops>=this_start & !(is.na(sumstats_file$midPos)) & sumstats_file$color!=this_col] = this_col
  }
  sumstats_file=unique(sumstats_file)
  write.table(sumstats_file,"rec_taj_dxy_fst_islswp_miss_lostruct.txt",
              col.names = T,row.names = F)
}
for(i in 1:nrow(colors)) {
  print(paste(i,nrow(colors)))
  this_start = colors$start[i]
  this_stop = colors$end[i]
  this_chrom = colors$chrom[i]
  this_species = colors$species[i]
  this_col = colors$ccols[i]
  sumstats_file$color[sumstats_file$species==this_species & sumstats_file$chr==this_chrom & sumstats_file$midPos>= this_start & sumstats_file$midPos<= this_stop & 
                        sumstats_file$windowstarts<=this_stop & sumstats_file$windowstops>=this_start] = this_col
}
sumstats_file$color[is.na(sumstats_file$color)] = "empty"
sumstats_file=unique(sumstats_file)
write.table(sumstats_file,"rec_taj_dxy_fst_islswp_miss_lostruct.temp",
            col.names = T,row.names = F)
sumstats_file = read.table("rec_taj_dxy_fst_islswp_miss_lostruct.temp",header=T)
pdf("temp.boxplot.colors.lostruct.pdf",width=10)
for(spp in unique(sumstats_file$species)){
  par(mfrow=c(3,5),mar=c(1,4,0,0))
  print(spp)
  temp = sumstats_file[sumstats_file$species==spp,]
  colors=sort(unique(temp$color))
  colors[colors=="black"]="grey"
  colors[colors=="empty"]="white"
  for(col_number in c(4,5,8,9,16,17,20:27)) {
    print(colnames(temp)[col_number])
    if (sum(complete.cases(temp[,c(col_number,28)]))>0) {
      boxplot(temp[,col_number]~temp$color,main=spp,ylab=colnames(temp)[col_number],
              las=1,col=colors,xlab="")
    }
  }
}
dev.off()
chromstarts=aggregate(sumstats_file$plotorder~sumstats_file$chr,FUN=function(x){min(x,na.rm=T)})
setwd("~")
for(spp in (unique(sumstats_file$species))){
  sppdf = sumstats_file[sumstats_file$species==spp,]
  sppdf=sppdf[sppdf$color!="empty",]
  print(spp)
  x = table(sppdf[,c("color","chr")])
  size=colSums(x,na.rm=T)
  size = size/100
  size[size<5] = 5
  x_norm = t(x)/colSums(x,na.rm = T)
  png(paste(spp,".lostructplot.png",sep=""),height=3,width=12,units="in",res=600)
  par(mar=c(3,2,1,0))
  barplot(t(x_norm),col=colnames(x_norm),las=2,width = (size),border=NA,beside=F,ylab="",cex.names = 0.75)
  dev.off()
}
png(paste("all.lostructplot.png",sep=""),height=3,width=10,units="in",res=600)
x = table(sumstats_file[,c("color","chr")])
size=colSums(x,na.rm=T)
size = size/1000
size[size<10] = 10
x_norm = t(x)/colSums(x,na.rm = T)
colnames(x_norm)[colnames(x_norm)=="empty"] = "grey"
chroms = c(1,"1A","1B",2:4,"4A",5:28,"LG2","LG5","LGE22","mtDNA","Z")
chroms= chroms[chroms %in% rownames(x_norm)]
par(mar=c(3,4,1,0))
barplot(t(x_norm),col=colnames(x_norm),las=2,width = (size),border=NA,beside=F,ylab="All Species",cex.names = 0.5)
dev.off()
chromstarts=aggregate(sumstats_file$plotorder~sumstats_file$chr,FUN=function(x){min(x,na.rm=T)})
for(spp in (unique(sumstats_file$species))) {
  sppdf = sumstats_file[sumstats_file$species==spp,]
  sppdf$color[sppdf$color=="empty"] = "grey"
  print(spp)
  pdf(paste(spp,".lostruct_chroms.pdf",sep=""))
  for(chrom in unique(sppdf$chr)) {
    chr = sppdf[sppdf$chr==chrom,]
    plot(x=chr$plotorder,y=as.numeric(as.factor(chr$color)),col=chr$color,cex=1,pch=16,main=chrom)
  }
  dev.off()
}
dev.off()
for(spp in (unique(sumstats_file$species))) {
  sppdf = sumstats_file[sumstats_file$species==spp,]
  sppdf$color[sppdf$color=="empty"] = "grey"
  print(spp)
  png(paste(spp,".lostruct_fst.png",sep=""),height=3,width=10,units="in",res=600)
  par(mar=c(3,4,1,0))
  if(spp=="sppdf") {
    plotrix::gap.plot(x=sppdf$plotorder[!is.na(sppdf$plotorder) & !(is.na(sppdf$Fst))], y=sppdf$Fst[!is.na(sppdf$plotorder)& !(is.na(sppdf$Fst))],gap.axis="y",gap=c(0.25,0.55),col=sppdf$color[!is.na(sppdf$plotorder)& !(is.na(sppdf$Fst))],cex=0.5,main="",pch=16,ylim=c(0,0.6),ytics=c(0,0.1,0.2,0.25,0.55,0.6),ylab="Fst")
    plotrix::gap.plot(x=sppdf$plotorder[!is.na(sppdf$plotorder) & !(is.na(sppdf$Fst)) & sppdf$color!="grey"], y=sppdf$Fst[!is.na(sppdf$plotorder) & !(is.na(sppdf$Fst)) & sppdf$color!="grey"],gap.axis="y",gap=c(0.25,0.55),add=T,col=sppdf$color[!is.na(sppdf$plotorder) & !(is.na(sppdf$Fst)) & sppdf$color!="grey"],cex=0.5,main="",pch=16,ylim=c(0,0.6),ylab="Fst")
    axis(2,at=0.25)
  } else {
    plot(x=sppdf$plotorder,y=sppdf$Fst,col=sppdf$color,cex=0.5,main=spp,pch=16,ylab="Fst")
    points(x=sppdf$plotorder[sppdf$color!="grey"],y=sppdf$Fst[sppdf$color!="grey"],col=sppdf$color[sppdf$color!="grey"],cex=0.5,main="",pch=16,ylab="Fst")
  }
  abline(v=as.numeric(chromstarts[,2]),col="red",lwd=0.5,lty=3)
  dev.off()
}
x=table(sumstats_file[,c("color","chr")])
(proportions(x[1:4,],margin=2)[4,])
proportions(table(sumstats_file$color)[1:4])
library(lme4)
pdf("~/color-fst.pdf")
for(spp in unique(sumstats_file$species)) {
  subset = sumstats_file[sumstats_file$species==spp,c("Fst","color")]
  subset=subset[subset$color!="empty",]
  mod = aov(subset$Fst~subset$color)
  print(summary(mod))
  TukeyHSD(mod)
  boxplot(subset$Fst~subset$color,main=spp)
}
dev.off()
