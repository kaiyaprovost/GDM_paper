rm(list=ls())
setwd("~")
library(RColorBrewer)
date=format(Sys.time(), "%d%b%Y")
df = read.table("bivariate_gdm_results_NGSDIST_20July2022.csv",
                sep=",",
                header=T,fill=T,
                stringsAsFactors = F,
                skip = 0)
colnames(df) = toupper(colnames(df))
plot(log10(df$CHROM_LENGTH)[df$DATASET!="GENOME"],df$HZAR_width[df$DATASET!="GENOME"],
     col=as.numeric(as.factor(df$SPECIES[df$DATASET!="GENOME"])))
mod=lm(df$HZAR_width[df$DATASET!="GENOME"]~log10(df$CHROM_LENGTH[df$DATASET!="GENOME"]))
mod=glm(df$HZAR_width[df$DATASET!="GENOME"]
        ~log10(df$CHROM_LENGTH[df$DATASET!="GENOME"])+df$SPECIES[df$DATASET!="GENOME"])
abline(mod,col="red")
summary(mod) 
par(mfrow=c(2,5))
for(spp in sort(unique(df$SPECIES))){
  print(spp)
  temp = df[df$SPECIES==spp,]
  temp = temp[,c("CHROM_LENGTH","DATASET","HZAR_WIDTH")]
  temp = temp[temp$DATASET!="GENOME",]
  temp=temp[complete.cases(temp),]
  plot(log10(temp$CHROM_LENGTH),temp$HZAR_WIDTH,main=spp)
  mod=lm(temp$HZAR_WIDTH~log10(temp$CHROM_LENGTH))
  abline(mod,col="red")
  print(summary(mod))
}
sumstats=c("CHROM_LENGTH","GDMCOLOR","HZAR_CENTER","HZAR_WIDTH","MEAN_DXY","MEAN_FST","MEAN_RECOMB","NUMBER.DXY.LOWS",
           "NUMBER.DXY.PEAKS","NUMBER.FST.PEAKS","NUMBER.ISLANDS","NUMBER.SWEEPS",
           "PERCENT_MISSING","PROP.DXY.LOW","PROP.DXY.PEAK","PROP.FST.PEAK","PROP.ISL","PROP.SWP","TAJIMAS")
sumstats = intersect(sumstats,colnames(df))
pdf(paste("summary_stats_by_species_gdm_",date,".pdf",sep="")); {
  cols=rainbow(10)
  for (stat in sumstats){
    print(stat)
    boxplot(df[,stat] ~ substr(df$SPECIES,1,3),col=cols,las=2,ylab=stat,xlab="Species",main="normal")
  }
}
dev.off()
pdf(paste("summary_stats_by_mostamodel2_gdm_",date,".pdf",sep="")); {
  cols=c(brewer.pal(8,"Dark2")[c(1,3,7,4)],"grey")
  for (stat in sumstats){
    print(stat)
    data_to_plot = df[,stat]
    boxplot(data_to_plot[df$MOSTA.1MODEL2!=""] ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],col=cols,las=2,ylab=stat,xlab="Species",main="normal")
  }
}
dev.off()
pdf(paste("summary_stats_by_mostamodel1_gdm_",date,".pdf",sep="")); {
  cols=c(brewer.pal(8,"Dark2")[c(1,2,3,7,4)],"grey")
  for (stat in sumstats){
    data_to_plot = df[,stat]
    boxplot(data_to_plot[df$MOSTA.1MODEL1!=""] ~ df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""],col=cols,las=2,ylab=stat,xlab="Species",main="normal")
  }
}
dev.off()
pdf(paste("summary_stats_by_mostamodel3_gdm_",date,".pdf")); {
  cols=c(brewer.pal(8,"Dark2")[c(1,7,4)],"grey")
  for (stat in sumstats){
    data_to_plot = df[,stat]
    boxplot(data_to_plot[df$MOSTA.1MODEL3!=""] ~ df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""],col=cols,las=2,ylab=stat,xlab="Species",main="normal")
  }
}
dev.off()
## check hzar 
boxplot(df$HZAR_CENTER[df$MOSTA.1MODEL1!=""]~df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""])
boxplot(df$HZAR_WIDTH[df$MOSTA.1MODEL1!=""]~df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""])
model1=aov(df$HZAR_WIDTH[df$MOSTA.1MODEL1!=""] ~ df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""])
summary(model1) 
TukeyHSD(model1) 
model1=aov(df$HZAR_CENTER[df$MOSTA.1MODEL1!=""] ~ df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""])
summary(model1) 
TukeyHSD(model1) 
png(paste("summary_stats_by_species_gdm_",date,".png",sep=""),height=4.5,width=6,units = "in",
    res=300); {
      par(mfrow=c(2,2),
          cex.axis=0.75,cex.lab=0.75,
          mar=c(0.5,4,0.5,0))
      cols=rainbow(10)
      plotrix::gap.boxplot(df$MEAN_FST_100~substr(df$SPECIES,1,3),
                           gap=list(top=c(0.22,0.56),bottom=c(NA,NA)),
                           axes=F,ylim=c(0.01,0.58),ylab="Mean Fst",
                           xaxt="n",col=cols)
      title(ylab="Mean Fst")
      axis(2,labels=c(seq(0,0.2,0.05),0.57),
           at=c(seq(0,0.2,0.05),0.235),tick=T)
      boxplot(df$MEAN_DXY_100~substr(df$SPECIES,1,3),ylab="Mean Dxy",
              xaxt="n",col=cols)
      boxplot(df$PERCENT_MISSING~substr(df$SPECIES,1,3),
              ylab="Percent Missing",xaxt="n",col=cols)
      boxplot((df$MEAN_RECOMB_100/(1e-9))~substr(df$SPECIES,1,3),
              ylab="Mean Recombination Rate (x 1^-9)",
              xaxt="n",col=cols)
    }
dev.off()
png(paste("summary_stats_by_species_gdmlostruct_",date,".png"),height=4.5,width=6,units = "in",
    res=300); {
      par(mfrow=c(2,2),
          cex.axis=0.75,cex.lab=0.75,
          mar=c(0.5,4,0.5,0))
      cols=rainbow(10)
      boxplot(lostruct$MEAN_FST~substr(lostruct$SPECIES,1,3),
              ylab="Mean Fst",
              xaxt="n",col=cols)
      title(ylab="Mean Fst")
      boxplot(lostruct$MEAN_DXY~substr(lostruct$SPECIES,1,3),ylab="Mean Dxy",
              xaxt="n",col=cols)
      boxplot(lostruct$PERCENT_MISSING~substr(lostruct$SPECIES,1,3),
              ylab="Percent Missing",xaxt="n",col=cols)
      boxplot((lostruct$MEAN_RECOMB/(1e-9))~substr(lostruct$SPECIES,1,3),
              ylab="Mean Recombination Rate (x 1^-9)",
              xaxt="n",col=cols)
    }
dev.off()
numonly = df[,c("MEAN_FST_100","MEAN_DXY_100","MEAN_RECOMB_100","PERCENT_MISSING")]
numonly_nofstout = numonly[numonly$MEAN_FST_100<=0.5,]
plot(numonly_nofstout)
cor(numonly,use="pairwise.complete.obs")
cor(numonly_nofstout,use="pairwise.complete.obs")
par(mfrow=c(1,2))
corrplot::corrplot(cor(numonly,use="pairwise.complete.obs"),
                   method = "number",diag=T,type="upper")
corrplot::corrplot(cor(numonly_nofstout,use="pairwise.complete.obs"),
                   method="number",diag=T,type="upper")
par(mfrow=c(2,3),mar=c(4,4,0,0))
plot(numonly$PERCENT_MISSING,numonly$MEAN_FST_100)
abline(lm(numonly$MEAN_FST_100~numonly$PERCENT_MISSING),col="red")
summary(lm(numonly$MEAN_FST_100~numonly$PERCENT_MISSING))$adj.r.squared
plot(numonly$PERCENT_MISSING,numonly$MEAN_DXY_100)
abline(lm(numonly$MEAN_DXY_100~numonly$PERCENT_MISSING),col="red")
summary(lm(numonly$MEAN_DXY_100~numonly$PERCENT_MISSING))$adj.r.squared
plot(numonly$PERCENT_MISSING,numonly$MEAN_RECOMB_100)
abline(lm(numonly$MEAN_RECOMB_100~numonly$PERCENT_MISSING),col="red")
summary(lm(numonly$MEAN_RECOMB_100~numonly$PERCENT_MISSING))$adj.r.squared
plot(numonly$MEAN_DXY_100,numonly$MEAN_FST_100T)
abline(lm(numonly$MEAN_FST_100~numonly$MEAN_DXY_100),col="red")
summary(lm(numonly$MEAN_FST_100~numonly$MEAN_DXY_100))$adj.r.squared
plot(numonly$MEAN_RECOMB_100,numonly$MEAN_FST_100)
abline(lm(numonly$MEAN_FST_100~numonly$MEAN_RECOMB_100),col="red")
summary(lm(numonly$MEAN_FST_100~numonly$MEAN_RECOMB_100))$adj.r.squared
plot(numonly$MEAN_RECOMB_100,numonly$MEAN_DXY_100)
abline(lm(numonly$MEAN_DXY_100~numonly$MEAN_RECOMB_100),col="red")
summary(lm(numonly$MEAN_DXY_100~numonly$MEAN_RECOMB_100))$adj.r.squared
par(mfrow=c(2,3),mar=c(4,4,0,0))
plot(numonly_nofstout$PERCENT_MISSING,numonly_nofstout$MEAN_FST_100)
abline(lm(numonly_nofstout$MEAN_FST_100~numonly_nofstout$PERCENT_MISSING),col="red")
summary(lm(numonly_nofstout$MEAN_FST_100~numonly_nofstout$PERCENT_MISSING))$adj.r.squared
plot(numonly_nofstout$PERCENT_MISSING,numonly_nofstout$MEAN_DXY_100)
abline(lm(numonly_nofstout$MEAN_DXY_100~numonly_nofstout$PERCENT_MISSING),col="red")
summary(lm(numonly_nofstout$MEAN_DXY_100~numonly_nofstout$PERCENT_MISSING))$adj.r.squared
plot(numonly_nofstout$PERCENT_MISSING,numonly_nofstout$MEAN_RECOMB_100)
abline(lm(numonly_nofstout$MEAN_RECOMB_100~numonly_nofstout$PERCENT_MISSING),col="red")
summary(lm(numonly_nofstout$MEAN_RECOMB_100~numonly_nofstout$PERCENT_MISSING))$adj.r.squared
plot(numonly_nofstout$MEAN_DXY_100,numonly_nofstout$MEAN_FST_100)
abline(lm(numonly_nofstout$MEAN_FST_100~numonly_nofstout$MEAN_DXY_100),col="red")
summary(lm(numonly_nofstout$MEAN_FST_100~numonly_nofstout$MEAN_DXY_100))$adj.r.squared
plot(numonly_nofstout$MEAN_RECOMB_100,numonly_nofstout$MEAN_FST_100T)
abline(lm(numonly_nofstout$MEAN_FST_100~numonly_nofstout$MEAN_RECOMB_100),col="red")
summary(lm(numonly_nofstout$MEAN_FST_100~numonly_nofstout$MEAN_RECOMB_100))$adj.r.squared
plot(numonly_nofstout$MEAN_RECOMB_100,numonly_nofstout$MEAN_DXY_100)
abline(lm(numonly_nofstout$MEAN_DXY_100~numonly_nofstout$MEAN_RECOMB_100),col="red")
summary(lm(numonly_nofstout$MEAN_DXY_100~numonly_nofstout$MEAN_RECOMB_100))$adj.r.squared
plot(df$MEAN_FST_100,df$PERCENT_MISSING)
plot(df$MEAN_DXY_100,df$PERCENT_MISSING)
plot(df$MEAN_RECOMB_100,df$PERCENT_MISSING)
cor(df$MEAN_FST_100,df$PERCENT_MISSING,use="pairwise.complete.obs") 
cor(df$MEAN_DXY_100,df$PERCENT_MISSING,use="pairwise.complete.obs")
cor(df$MEAN_RECOMB_100,df$PERCENT_MISSING,use="pairwise.complete.obs") 
library(RColorBrewer)
cols=c(brewer.pal(8,"Dark2")[c(1,3,7,4)],"grey")
par(mar=c(4,4,0.1,0.1),
    mfrow=c(1,3))
boxplot(df$MEAN_RECOMB_100[df$MOSTA.1MODEL2!=""]/(1e-9)~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
        col=cols,
        ylab="Mean Recombination Rate (x 1^-9)",
        xlab="Best Model",las=2)
boxplot(df$MEAN_RECOMB_100[df$MOSTA.1MODEL1!=""]/(1e-9)~df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""],
        col=cols,
        ylab="Mean Recombination Rate (x 1^-9)",
        xlab="Best Model",las=2)
boxplot(df$MEAN_RECOMB_100[df$MOSTA.1MODEL3!=""]/(1e-9)~df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""],
        col=cols,
        ylab="Mean Recombination Rate (x 1^-9)",
        xlab="Best Model",las=2)
model=aov(df$MEAN_RECOMB_100[df$MOSTA.1MODEL1!=""] ~ df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""])
summary(model) 
TukeyHSD(model)
agg=aggregate(df$PERCENT_MISSING~df$SPECIES,FUN=function(x){sd(x,na.rm=T)})
png(paste("four_panel_figure_chapter_3_ibd2mix_",date,".png"),height=4,width=6,units = "in",
    res=300); {
      library(RColorBrewer)
      cols=c(brewer.pal(8,"Dark2")[c(1,3,7,4)],"grey")
      par(mfrow=c(2,2),cex.axis=1)
      par(mar=c(0.1,4,0.3,0.1))
      plotrix::gap.boxplot(df$MEAN_FST_100[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
                           gap=list(top=c(0.22,0.56),bottom=c(NA,NA)),
                           axes=F,ylim=c(0.01,0.58),ylab="Mean Fst",xlab="Best Model",
                           xaxt="n",col=cols)
      legend("topright",legend=c("IBA","IBD","IBE","IBH","MIX"),
             col=cols,fill=cols,cex=0.8,ncol=1)
      title(ylab="Mean Fst")
      axis(2,labels=c(seq(0,0.2,0.05),0.57),
           at=c(seq(0,0.2,0.05),0.235),tick=T,las=1)
      boxplot(df$MEAN_DXY_100[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
              col=cols,xaxt="n",
              ylab="Mean Dxy",
              xlab="Best Model",las=1)
      boxplot(df$PERCENT_MISSING[df$MOSTA.1MODEL2!=""]*100~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
              col=cols,xaxt="n",
              ylab="Percent Missing",
              xlab="Best Model",las=1)
      boxplot(df$MEAN_RECOMB_100[df$MOSTA.1MODEL2!=""]/(1e-9)~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
              col=cols,xaxt="n",
              ylab="Mean Recombination Rate",
              xlab="Best Model",las=2)
    }
dev.off()
png(paste("four_panel_figure_chapter_3_ibd2mix_lostruct_",date,".png"),height=4,width=6,units = "in",
    res=300); {
      library(RColorBrewer)
      cols=c(brewer.pal(8,"Dark2")[c(1,3,7,4)],"grey")
      par(mfrow=c(2,2),cex.axis=1)
      par(mar=c(0.1,4,0.3,0.1))
      boxplot(lostruct$MEAN_FST[lostruct$MOSTA.1MODEL2!=""]~lostruct$MOSTA.1MODEL2[lostruct$MOSTA.1MODEL2!=""],
              ylab="Mean Fst",xlab="Best Model",
              xaxt="n",col=cols,border=c(rep("black",4),cols[5]))
      legend("topright",legend=c("IBA","IBD","IBE","IBH","MIX"),
             col=cols,fill=cols,cex=0.8,ncol=1)
      b=boxplot(lostruct$MEAN_DXY[lostruct$MOSTA.1MODEL2!=""]~lostruct$MOSTA.1MODEL2[lostruct$MOSTA.1MODEL2!=""],
                col=cols,border=c(rep("black",4),cols[5]),xaxt="n",
                ylab="Mean Dxy",
                xlab="Best Model",las=1)
      boxplot(lostruct$PERCENT_MISSING[lostruct$MOSTA.1MODEL2!=""]*100~lostruct$MOSTA.1MODEL2[lostruct$MOSTA.1MODEL2!=""],
              col=cols[1:5],border=c("black","black","black","black",cols[5]),
              na.action="na.pass",
              ylab="Percent Missing",xaxt="n",
              xlab="Best Model",las=1)
      boxplot(lostruct$MEAN_RECOMB[lostruct$MOSTA.1MODEL2!=""]/(1e-9)~lostruct$MOSTA.1MODEL2[lostruct$MOSTA.1MODEL2!=""],
              col=cols,border=c(rep("black",4),cols[5]),xaxt="n",
              ylab="Mean Recombination Rate",
              xlab="Best Model",las=2)
    }
dev.off()
png(paste("four_panel_figure_chapter_3_univariate_ibd2mix_lostruct_",date,".png"),height=4,width=6,units = "in",
    res=300); {
      library(RColorBrewer)
      cols=c(brewer.pal(8,"Dark2")[c(1,2,3,7,4)],"grey")
      par(mfrow=c(2,2),cex.axis=1)
      par(mar=c(0.1,4,0.3,0.1))
      boxplot(lostruct$MEAN_FST[lostruct$MOSTA.1MODEL1!=""]~lostruct$MOSTA.1MODEL1[lostruct$MOSTA.1MODEL1!=""],
              ylab="Mean Fst",xlab="Best Model",
              xaxt="n",col=cols,border=c(rep("black",4),cols[5:6]))
      legend("topright",legend=c("IBA","IBB","IBD","IBE","IBH","MIX"),
             col=cols,fill=cols,cex=0.8,ncol=1)
      b=boxplot(lostruct$MEAN_DXY[lostruct$MOSTA.1MODEL1!=""]~lostruct$MOSTA.1MODEL1[lostruct$MOSTA.1MODEL1!=""],
                col=cols,border=c(rep("black",5),cols[6]),xaxt="n",
                ylab="Mean Dxy",
                xlab="Best Model",las=1)
      boxplot(lostruct$PERCENT_MISSING[lostruct$MOSTA.1MODEL1!=""]*100~lostruct$MOSTA.1MODEL1[lostruct$MOSTA.1MODEL1!=""],
              col=cols[1:5],border=c("black","black","black","black","black",cols[6]),
              na.action="na.pass",
              ylab="Percent Missing",xaxt="n",
              xlab="Best Model",las=1)
      boxplot(lostruct$MEAN_RECOMB[lostruct$MOSTA.1MODEL1!=""]/(1e-9)~lostruct$MOSTA.1MODEL1[lostruct$MOSTA.1MODEL1!=""],
              col=cols,border=c(rep("black",5),cols[6]),xaxt="n",
              ylab="Mean Recombination Rate",
              xlab="Best Model",las=2)
    }
dev.off()
png(paste("four_panel_figure_chapter_3_univariate_ibd2mix_",date,".png"),height=4,width=6,units = "in",
    res=300); {
      library(RColorBrewer)
      cols=c(brewer.pal(8,"Dark2")[c(1,2,3,7,4)],"grey")
      par(mfrow=c(2,2))
      par(mar=c(0.1,4,0.3,0.1))
      plotrix::gap.boxplot(df$MEAN_FST[df$MOSTA.1MODEL1!=""]~df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""],
                           gap=list(top=c(0.22,0.56),bottom=c(NA,NA)),
                           axes=F,ylim=c(0.01,0.58),ylab="Mean Fst",xlab="Best Model",
                           xaxt="n",col=cols)
      legend("topright",legend=c("IBA","IBB","IBD","IBE","IBH","MIX"),
             col=cols,fill=cols,cex=0.8,ncol=1)
      title(ylab="Mean Fst")
      axis(2,labels=c(seq(0,0.2,0.05),0.57),
           at=c(seq(0,0.2,0.05),0.235),tick=T,las=1)
      boxplot(df$MEAN_DXY[df$MOSTA.1MODEL1!=""]~df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""],
              col=cols,xaxt="n",
              ylab="Mean Dxy",
              xlab="Best Model",las=2)
      boxplot(df$PERCENT_MISSING[df$MOSTA.1MODEL1!=""]*100~df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""],
              col=cols,xaxt="n",
              ylab="Percent Missing",
              xlab="Best Model",las=2)
      boxplot(df$MEAN_RECOMB[df$MOSTA.1MODEL1!=""]/(1e-9)~df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""],
              col=cols,xaxt="n",
              ylab="Mean Recombination Rate (x 1^-9)",
              xlab="Best Model",las=2)
    }
dev.off()
png(paste("four_panel_figure_chapter_3_multivariate_ibd2mix_lostruct_",date,".png"),height=4,width=6,units = "in",
    res=300); {
      library(RColorBrewer)
      cols=c(brewer.pal(8,"Dark2")[c(1,7,4)],"grey")
      par(mfrow=c(2,2),cex.axis=1)
      par(mar=c(0.1,4,0.3,0.1))
      boxplot(lostruct$MEAN_FST[lostruct$MOSTA.1MODEL3!=""]~lostruct$MOSTA.1MODEL3[lostruct$MOSTA.1MODEL3!=""],
              ylab="Mean Fst",xlab="Best Model",
              xaxt="n",col=cols,border=c(rep("black",3),cols[4]))
      legend("topright",legend=c("IBA","IBE","IBH","MIX"),
             col=cols,fill=cols,cex=0.8,ncol=1)
      b=boxplot(lostruct$MEAN_DXY[lostruct$MOSTA.1MODEL3!=""]~lostruct$MOSTA.1MODEL3[lostruct$MOSTA.1MODEL3!=""],
                col=cols,border=c(rep("black",5),cols[6]),xaxt="n",
                ylab="Mean Dxy",
                xlab="Best Model",las=1)
      boxplot(lostruct$PERCENT_MISSING[lostruct$MOSTA.1MODEL3!=""]*100~lostruct$MOSTA.1MODEL3[lostruct$MOSTA.1MODEL3!=""],
              col=cols[1:5],border=c("black","black","black",cols[4]),
              na.action="na.pass",
              ylab="Percent Missing",xaxt="n",
              xlab="Best Model",las=1)
      boxplot(lostruct$MEAN_RECOMB[lostruct$MOSTA.1MODEL3!=""]/(1e-9)~lostruct$MOSTA.1MODEL3[lostruct$MOSTA.1MODEL3!=""],
              col=cols,border=c(rep("black",5),cols[6]),xaxt="n",
              ylab="Mean Recombination Rate",
              xlab="Best Model",las=2)
    }
dev.off()
png(paste("four_panel_figure_chapter_3_multivariate_",date,".png"),height=4,width=6,units = "in",
    res=300); {
      library(RColorBrewer)
      cols=c(brewer.pal(8,"Dark2")[c(1,7,4)],"grey")
      par(mfrow=c(2,2))
      par(mar=c(0.1,4,0.3,0.1))
      plotrix::gap.boxplot(df$MEAN_FST[df$MOSTA.1MODEL3!=""]~df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""],
                           gap=list(top=c(0.22,0.56),bottom=c(NA,NA)),
                           axes=F,ylim=c(0.01,0.58),ylab="Mean Fst",xlab="Best Model",
                           xaxt="n",col=cols)
      legend("topright",legend=c("IBA","IBE","IBH","MIX"),
             col=cols,fill=cols,cex=0.8,ncol=1)
      title(ylab="Mean Fst")
      axis(2,labels=c(seq(0,0.2,0.05),0.57),
           at=c(seq(0,0.2,0.05),0.235),tick=T,las=1)
      boxplot(df$MEAN_DXY[df$MOSTA.1MODEL3!=""]~df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""],
              col=cols,xaxt="n",
              ylab="Mean Dxy",
              xlab="Best Model",las=2)
      boxplot(df$PERCENT_MISSING[df$MOSTA.1MODEL3!=""]*100~df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""],
              col=cols,xaxt="n",
              ylab="Percent Missing",
              xlab="Best Model",las=2)
      boxplot(df$MEAN_RECOMB[df$MOSTA.1MODEL3!=""]/(1e-9)~df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""],
              col=cols,xaxt="n",
              ylab="Mean Recombination Rate (x 1^-9)",
              xlab="Best Model",las=2)
    }
dev.off()
## aov model -- BIVARIATE
model1=aov(df$MEAN_RECOMB[df$MOSTA.1MODEL2!=""] ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
model2=aov(df$MEAN_FST[df$MOSTA.1MODEL2!=""] ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
model3=aov(df$MEAN_DXY[df$MOSTA.1MODEL2!=""] ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
model4=aov(df$PERCENT_MISSING[df$MOSTA.1MODEL3!=""] ~ df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""])
model4n=aov(df$PERCENT_MISSING[df$MOSTA.1MODEL3!="" & df$SPECIES!="NITENS"] ~ df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!="" & df$SPECIES!="NITENS"])
model5=aov(log(df$CHROM_LENGTH[df$MOSTA.1MODEL2!=""]) ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
summary(model1) 
summary(model2) 
summary(model3) 
summary(model4)
summary(model4n) 
summary(model5) 
TukeyHSD(model1) 
TukeyHSD(model2) 
TukeyHSD(model3) 
TukeyHSD(model4) 
TukeyHSD(model4n)
TukeyHSD(model5)
## aov model -- UNIVARIATE
model1=aov(df$MEAN_RECOMB_100[df$MOSTA.1MODEL1!=""] ~ df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""])
model2=aov(df$MEAN_FST_100[df$MOSTA.1MODEL1!=""] ~ df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""])
model3=aov(df$MEAN_DXY_100[df$MOSTA.1MODEL1!=""] ~ df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""])
model4=aov(df$PERCENT_MISSING[df$MOSTA.1MODEL1!=""] ~ df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""])
model5=aov(log(df$CHROM_LENGTH[df$MOSTA.1MODEL1!=""]) ~ df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""])
summary(model1) 
summary(model2) 
summary(model3) 
summary(model4) 
summary(model5) 
TukeyHSD(model1) 
TukeyHSD(model2) 
TukeyHSD(model3) 
TukeyHSD(model4) 
TukeyHSD(model5) 
boxplot(df$CONTACT_SUITABILITY_SD[df$MOSTA.1MODEL1!=""]~ df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""])
aggregate(df$CONTACT_SUITABILITY_SD~df$MOSTA.1MODEL1,FUN=function(x){mean(x,na.rm=T)})
tab=table(df$CONTACT_SUITABILITY_SD,df$MOSTA.1MODEL1)
par(mfrow=c(2,3))
plot(as.numeric(rownames(tab)),tab[,"IBA"],type="b")
plot(as.numeric(rownames(tab)),tab[,"IBB"],type="b")
plot(as.numeric(rownames(tab)),tab[,"IBD"],type="b")
plot(as.numeric(rownames(tab)),tab[,"IBE"],type="b")
plot(as.numeric(rownames(tab)),tab[,"IBH"],type="b")
plot(as.numeric(rownames(tab)),tab[,"MIXED"],type="b")
summary(lm(tab[,"IBA"]~as.numeric(rownames(tab))))
summary(lm(tab[,"IBB"]~as.numeric(rownames(tab))))
summary(lm(tab[,"IBD"]~as.numeric(rownames(tab))))
summary(lm(tab[,"IBE"]~as.numeric(rownames(tab))))
summary(lm(tab[,"IBH"]~as.numeric(rownames(tab))))
summary(lm(tab[,"MIXED"]~as.numeric(rownames(tab))))
model6=aov((df$CONTACT_SUITABILITY_MEAN[df$MOSTA.1MODEL1!=""]) ~ df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""])
summary(model6) 
TukeyHSD(model6)
model6=aov((df$CONTACT_SUITABILITY_SD[df$MOSTA.1MODEL1!=""]) ~ df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""])
summary(model6) 
TukeyHSD(model6)
model6=aov((df$TAJIMAS[df$MOSTA.1MODEL1!=""]) ~ df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""])
summary(model6) 
TukeyHSD(model6)
dfpca=df[,c("CONTACT_SUITABILITY_MEAN","CONTACT_SUITABILITY_SD",
            "HZAR_CENTER","HZAR_WIDTH","MEAN_DXY_100","MEAN_FST_100",
            "MEAN_RECOMB_100","PERCENT_MISSING","TAJIMAS")]
dfpca=dfpca[complete.cases(dfpca),]
prcomp()
## aov model -- TRIVARIATE
model1=aov(df$MEAN_RECOMB[df$MOSTA.1MODEL3!=""] ~ df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""])
model2=aov(df$MEAN_FST[df$MOSTA.1MODEL3!=""] ~ df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""])
model3=aov(df$MEAN_DXY[df$MOSTA.1MODEL3!=""] ~ df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""])
model4=aov(df$PERCENT_MISSING[df$MOSTA.1MODEL3!=""] ~ df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""])
model5=aov(log(df$CHROM_LENGTH[df$MOSTA.1MODEL3!=""]) ~ df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""])
summary(model1) 
summary(model2) 
summary(model3) 
summary(model4)
summary(model5) 
TukeyHSD(model1) 
TukeyHSD(model2) 
TukeyHSD(model3) 
TukeyHSD(model4) 
TukeyHSD(model5) 
png(paste("four_panel_figure_chapter_3_nothypotheses_",date,".png"),height=4,width=6,units = "in",
    res=300); {
      library(RColorBrewer)
      cols=c(brewer.pal(8,"Dark2")[c(1,7,3,4)],"grey",brewer.pal(8,"Dark2")[c(7)])
      par(mfrow=c(2,2))
      par(mar=c(4,4,0,0))
      boxplot(df$MEAN_FST[df$MOST.1MODEL2!=""]~df$MOST.1MODEL2[df$MOST.1MODEL2!=""],
              col=cols,
              ylab="Mean Fst (cut off)",ylim=c(0,0.21),
              xlab="Best Model",las=2)
      boxplot(df$MEAN_DXY[df$MOST.1MODEL2!=""]~df$MOST.1MODEL2[df$MOST.1MODEL2!=""],
              col=cols,
              ylab="Mean Dxy",
              xlab="Best Model",las=2)
      boxplot(df$PERCENT_MISSING[df$MOST.1MODEL2!=""]*100~df$MOST.1MODEL2[df$MOST.1MODEL2!=""],
              col=cols,
              ylab="Percent Missing",
              xlab="Best Model",las=2)
      boxplot(df$MEAN_RECOMB[df$MOST.1MODEL2!=""]/(1e-9)~df$MOST.1MODEL2[df$MOST.1MODEL2!=""],
              col=cols,
              ylab="Mean Recombination Rate (x 1^-9)",
              xlab="Best Model",las=2)
    }
dev.off()
## aov model -- bivariate no hypotheses
model1=aov(df$MEAN_RECOMB[df$MOST.1MODEL2!=""] ~ df$MOST.1MODEL2[df$MOST.1MODEL2!=""])
model2=aov(df$MEAN_FST[df$MOST.1MODEL2!=""] ~ df$MOST.1MODEL2[df$MOST.1MODEL2!=""])
MODEL3=aov(df$MEAN_DXY[df$MOST.1MODEL2!=""] ~ df$MOST.1MODEL2[df$MOST.1MODEL2!=""])
model4=aov(df$PERCENT_MISSING[df$MOST.1MODEL2!=""] ~ df$MOST.1MODEL2[df$MOST.1MODEL2!=""])
summary(model1) 
summary(model2) 
summary(MODEL3) 
summary(model4)
TukeyHSD(model1) 
TukeyHSD(model2) 
TukeyHSD(MODEL3) 
TukeyHSD(model4) 
boxplot(df$PERCENT_MISSING[df$MOSTA.1MODEL2!="" & df$SPECIES!="NITENS"]*100
        ~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!="" & df$SPECIES!="NITENS"],
        col=cols,
        ylab="Percent Missing",
        xlab="Best Model",las=2)
model4=aov(df$PERCENT_MISSING[df$MOSTA.1MODEL2!="" & df$SPECIES!="NITENS"] 
           ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!="" & df$SPECIES!="NITENS"])
summary(model4) 
TukeyHSD(model4) 
df_numb = df[,c("PERCENT_MISSING","MEAN_FST","MEAN_DXY","MEAN_RECOMB")]
df_meta = df[,c("DATASET","SPECIES","MOST.1MODEL2","MOSTA.1MODEL2")]
df_meta = df_meta[complete.cases(df_numb),]
df_numb = df_numb[complete.cases(df_numb),]
pca = prcomp(df_numb,center = T,scale. = T)
summary(pca)
df_pca = cbind(df_meta,pca$x)
palette(c(brewer.pal(8,"Dark2")[c(1,7,3,4)],"grey"))
plot(df_pca$PC2,df_pca$PC1,col=as.numeric(as.factor(df_pca$MOSTA.1MODEL2)),
     pch=as.numeric(as.factor(df_pca$MOSTA.1MODEL2)),ylim=c(-3,2))
## look at tajima's 
par(mfrow=c(1,3),mar=c(4,4,0,0))
cols=c(brewer.pal(8,"Dark2")[c(1,7,3,4)],"grey",brewer.pal(8,"Dark2")[c(7)])
boxplot(df$TAJIMAS[df$MOST.1MODEL1!=""]~df$MOST.1MODEL1[df$MOST.1MODEL1!=""],
        col=cols,
        ylab="Tajima's D",
        xlab="",las=2)
boxplot(df$TAJIMAS[df$MOST.1MODEL2!=""]~df$MOST.1MODEL2[df$MOST.1MODEL2!=""],
        col=cols,
        ylab="",
        xlab="Best Model",las=2)
boxplot(df$TAJIMAS[df$MOST.1MODEL3!=""]~df$MOST.1MODEL3[df$MOST.1MODEL3!=""],
        col=cols,
        ylab="",
        xlab="",las=2)
png("tajimatest.png")
par(mfrow=c(1,3),mar=c(4,4,0,0))
cols=c(brewer.pal(8,"Dark2")[c(1,3,7,4)],"grey")
boxplot(df$TAJIMAS[df$MOSTA.1MODEL1!=""]~df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""],
        col=cols,
        ylab="Tajima's D",
        xlab="",las=2)
abline(h=0)
boxplot(df$TAJIMAS[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
        col=cols,
        ylab="",
        xlab="Best Model",las=2)
abline(h=0)
boxplot(df$TAJIMAS[df$MOSTA.1MODEL3!=""]~df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""],
        col=cols,
        ylab="",
        xlab="",las=2)
abline(h=0)
dev.off()
boxplot(df$TAJIMAS[df$MOST.1MODEL2!=""]~df$MOST.1MODEL2[df$MOST.1MODEL2!=""],
        col=cols,
        ylab="",
        xlab="Best Model",las=2)
modelA=aov(df$TAJIMAS[df$MOST.1MODEL1!=""] ~ df$MOST.1MODEL1[df$MOST.1MODEL1!=""])
modelA=aov(df$TAJIMAS[df$MOSTA.1MODEL1!=""] ~ df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""])
modelB=aov(df$TAJIMAS[df$MOST.1MODEL2!=""] ~ df$MOST.1MODEL2[df$MOST.1MODEL2!=""])
modelB=aov(df$TAJIMAS[df$MOSTA.1MODEL2!=""] ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
modelC=aov(df$TAJIMAS[df$MOST.1MODEL3!=""] ~ df$MOST.1MODEL3[df$MOST.1MODEL3!=""])
modelC=aov(df$TAJIMAS[df$MOSTA.1MODEL3!=""] ~ df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""])
summary(modelA) 
summary(modelB) 
summary(modelC) 
TukeyHSD(modelA)  
TukeyHSD(modelB) 
TukeyHSD(modelC) 
## negative tajd = lots of rare alleles, recent sweep, expansion after bottleneck, or linkage
## positive tajd = few rares alleles, balancing selection, sudden contraction
cor((df[,c("TAJIMAS","PERCENT_MISSING","MEAN_FST","MEAN_DXY","MEAN_RECOMB")]),
    use = "pairwise.complete.obs")
## chrom lengths
png(paste("chromlengthtest",date,".png",sep=""))
library(RColorBrewer)
cols=c(brewer.pal(8,"Dark2")[c(1,3,7,4)],"grey")
plotrix::gap.boxplot(df$CHROM_LENGTH[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],col=cols,
                     gap=list(top=c(16000,111000),bottom=c(NA,NA)))
title(ylab="Num Scaffolds")
dev.off()
boxplot(log(df$CHROM_LENGTH[df$MOSTA.1MODEL2!=""]*100000)
        ~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
        col=cols,xlab="Best Model",ylab="Log Chromosome Length")
df_nogen = df[df$DATASET!="GENOME",]
modelA=aov(df$CHROM_LENGTH[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
summary(modelA)
TukeyHSD(modelA) 
modelB=aov(df_nogen$CHROM_LENGTH[df_nogen$MOSTA.1MODEL2!=""]
           ~df_nogen$MOSTA.1MODEL2[df_nogen$MOSTA.1MODEL2!=""])
summary(modelB) 
TukeyHSD(modelB) 
modelC=aov(log(df$CHROM_LENGTH[df$MOSTA.1MODEL2!=""]*100000)
           ~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
summary(modelC)
TukeyHSD(modelC) 
## looking at new stats
pdf(paste("gdm_bivariate_all_sumstats_",date,".pdf",sep=""),height=4,width=6); {
  library(RColorBrewer)
  cols=c(brewer.pal(8,"Dark2")[c(1,3,7,4)],"grey")
  par(mfrow=c(2,2),cex.axis=1)
  par(mar=c(0.1,4,0.3,0.1))
  plotrix::gap.boxplot(df$MEAN_FST[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
                       gap=list(top=c(0.22,0.56),bottom=c(NA,NA)),
                       axes=F,ylim=c(0.01,0.58),ylab="Mean Fst",xlab="Best Model",
                       xaxt="n",col=cols)
  legend("topright",legend=c("IBA","IBD","IBE","IBH","MIX"),
         col=cols,fill=cols)
  title(ylab="Mean Fst")
  axis(2,labels=c(seq(0,0.2,0.05),0.57),
       at=c(seq(0,0.2,0.05),0.235),tick=T,las=1)
  boxplot(df$MEAN_DXY[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
          col=cols,xaxt="n",
          ylab="Mean Dxy",
          xlab="Best Model",las=1)
  boxplot(df$PERCENT_MISSING[df$MOSTA.1MODEL2!=""]*100~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
          col=cols,xaxt="n",
          ylab="Percent Missing",
          xlab="Best Model",las=1)
  boxplot(df$MEAN_RECOMB[df$MOSTA.1MODEL2!=""]/(1e-9)~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
          col=cols,xaxt="n",
          ylab="Mean Recombination Rate",
          xlab="Best Model",las=2)
  boxplot(df$NUMBER.FST.PEAKS[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
          col=cols,xaxt="n",
          ylab="N FST Peaks",
          xlab="Best Model",las=2)
  boxplot(df$NUMBER.FST.PEAKS[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
          col=cols,xaxt="n",ylim=c(0,500),
          ylab="N FST Peaks (Cropped)",
          xlab="Best Model",las=2)
  boxplot(df$NUMBER.DXY.PEAKS[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
          col=cols,xaxt="n",
          ylab="N DXY Peaks",
          xlab="Best Model",las=2)
  boxplot(df$NUMBER.DXY.PEAKS[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
          col=cols,xaxt="n",ylim=c(0,500),
          ylab="N DXY Peaks (Cropped)",
          xlab="Best Model",las=2)
  boxplot(df$NUMBER.ISLANDS[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
          col=cols,xaxt="n",
          ylab="N Islands",
          xlab="Best Model",las=2)
  boxplot(df$TAJIMAS[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
          col=cols,xaxt="n",
          ylab="Tajima's D",
          xlab="Best Model",las=2)
  boxplot(df$NUMBER.SWEEPS[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
          col=cols,xaxt="n",
          ylab="N Sweeps",
          xlab="Best Model",las=2)
  boxplot(df$NUMBER.SWEEPS[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
          col=cols,xaxt="n",ylim=c(0,100),
          ylab="N Sweeps (Cropped)",
          xlab="Best Model",las=2)
  boxplot(df$NUMBER.DXY.LOWS[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
          col=cols,xaxt="n",
          ylab="N Dxy Lows",
          xlab="Best Model",las=2)
  boxplot(df$NUMBER.DXY.LOWS[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
          col=cols,xaxt="n",ylim=c(0,500),
          ylab="N Dxy Lows (Cropped)",
          xlab="Best Model",las=2)
  boxplot(df$CHROM_LENGTH[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
          col=cols,xaxt="n",ylim=c(0,20000),
          ylab="Chrom Length (not Genome)",
          xlab="Best Model",las=2)
  boxplot(df$PROP.FST.PEAK[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
          col=cols,xaxt="n",
          ylab="Proportion FST Peak",
          xlab="Best Model",las=2)
  boxplot(df$PROP.DXY.PEAK[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
          col=cols,xaxt="n",
          ylab="Proportion DXY Peak",
          xlab="Best Model",las=2)
  boxplot(df$PROP.ISL[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
          col=cols,xaxt="n",
          ylab="Proportion Island",
          xlab="Best Model",las=2)
  boxplot(df$PROP.SWP[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
          col=cols,xaxt="n",
          ylab="Proportion Sweep",
          xlab="Best Model",las=2)
  boxplot(df$PROP.SWP[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
          col=cols,xaxt="n",ylim=c(0,0.01),
          ylab="Proportion Sweep (Cropped)",
          xlab="Best Model",las=2)
  boxplot(df$PROP.DXY.LOW[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
          col=cols,xaxt="n",
          ylab="Proportion Low DXY",
          xlab="Best Model",las=2)
  boxplot(df$PROP.DXY.LOW[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
          col=cols,xaxt="n",ylim=c(0,0.35),
          ylab="Proportion Low DXY (Cropped)",
          xlab="Best Model",las=2)
}
dev.off()
small=df[df$MOSTA.1MODEL2!="",c(4:19,33)]
agg=aggregate(cbind(small[,1],small[,2],small[,3],
                    small[,4],small[,5],small[,6],
                    small[,7],small[,8],small[,9],
                    small[,10],small[,11],small[,12],
                    small[,13],small[,14],small[,15],small[,16])
              ~small[,17],FUN=function(x){mean(x,na.rm=T)})
colnames(agg) = colnames(small)[c(17,1:16)]
pdf(paste("gdm_bivariate_all_barplots_",date,".pdf",sep=""),height=4,width=6); {
  library(RColorBrewer)
  cols=c(brewer.pal(8,"Dark2")[c(1,3,7,4)],"grey")
  par(mar=c(0.1,4,0.3,0.1))
  for(i in 2:17){
    print(i)
    barplot(agg[,i],col=cols,ylab=colnames(agg)[i])
  }
}
dev.off()
###
df2 = df[,c("DATASET","SPECIES","MEAN_DXY_50","MEAN_DXY_75","MEAN_DXY_100","MEAN_FST_50","MEAN_FST_75","MEAN_FST_100")]
df2_g = df2[df2$DATASET=="GENOME",]
df2_c = df2[grepl("CHR",df2$DATASET),]
df2_l = df2[grepl("LS",df2$DATASET),]
df2_h = df2[grepl("HIGH",df2$DATASET),]
df2_o = df2[grepl("LOW",df2$DATASET),]
agg_g_m = aggregate(cbind(df2_g$MEAN_DXY_50,df2_g$MEAN_DXY_75,df2_g$MEAN_DXY_100,
                          #df2_g$MEAN_FST_50,df2_g$MEAN_FST_75,
                          df2_g$MEAN_FST_100)~
                      df2_g$SPECIES,data=df2_g,FUN=function(x){mean(x,na.rm=T)})
agg_g_s = aggregate(cbind(df2_g$MEAN_DXY_50,df2_g$MEAN_DXY_75,df2_g$MEAN_DXY_100,
                          #df2_g$MEAN_FST_50,df2_g$MEAN_FST_75,
                          df2_g$MEAN_FST_100)~
                      df2_g$SPECIES,data=df2_g,FUN=function(x){sd(x,na.rm=T)})
barplot(t(as.matrix(agg_g_m[,2:4])),beside=T,col=c("black","darkgrey","lightgrey"),
        names.arg=agg_g_m[,1],las=2,
        main="DXY genome 50-75-100")
barplot(t(as.matrix(agg_g_m[,5])),beside=T,col=c("black","darkgrey","lightgrey"),
        names.arg=agg_g_m[,1],las=2,
        main="fst genome 100")
agg_c_m = aggregate(cbind(df2_c$MEAN_DXY_50,df2_c$MEAN_DXY_75,df2_c$MEAN_DXY_100,
                          df2_c$MEAN_FST_50,df2_c$MEAN_FST_75,
                          df2_c$MEAN_FST_100)~
                      df2_c$SPECIES,data=df2_c,FUN=function(x){mean(x,na.rm=T)})
agg_c_s = aggregate(cbind(df2_c$MEAN_DXY_50,df2_c$MEAN_DXY_75,df2_c$MEAN_DXY_100,
                          df2_c$MEAN_FST_50,df2_c$MEAN_FST_75,
                          df2_c$MEAN_FST_100)~
                      df2_c$SPECIES,data=df2_c,FUN=function(x){sd(x,na.rm=T)})
barplot(t(as.matrix(agg_c_m[,2:4])),beside=T,col=c("black","darkgrey","lightgrey"),
        names.arg=agg_c_m[,1],las=2,
        main="DXY chrom 50-75-100")
barplot(t(as.matrix(agg_c_m[,5:7])),beside=T,col=c("black","darkgrey","lightgrey"),
        names.arg=agg_c_m[,1],las=2,
        main="fst chrom 50-75-100")
agg_l_m = aggregate(cbind(
  #df2_l$MEAN_DXY_50,df2_l$MEAN_DXY_75,
  df2_l$MEAN_DXY_100,
  #df2_l$MEAN_FST_50,df2_l$MEAN_FST_75,
  df2_l$MEAN_FST_100)~
    df2_l$SPECIES,data=df2_l,FUN=function(x){mean(x,na.rm=T)})
agg_l_s = aggregate(cbind(
  #df2_l$MEAN_DXY_50,df2_l$MEAN_DXY_75,
  df2_l$MEAN_DXY_100,
  #df2_l$MEAN_FST_50,df2_l$MEAN_FST_75,
  df2_l$MEAN_FST_100)~
    df2_l$SPECIES,data=df2_l,FUN=function(x){sd(x,na.rm=T)})
barplot(t(as.matrix(agg_l_m[,2])),beside=T,col=c("black","darkgrey","lightgrey"),
        names.arg=agg_l_m[,1],las=2,
        main="DXY lostruct 100")
barplot(t(as.matrix(agg_l_m[,3])),beside=T,col=c("black","darkgrey","lightgrey"),
        names.arg=agg_l_m[,1],las=2,
        main="fst lostruct 100")
## correlation with chromosome size
df_chr = df[!(is.na(df$CHROM_LENGTH)),]
df_chr = df_chr[df_chr$FSTOUTLIER==0,]
df_chr = df_chr[df_chr$LOSTRUCTOUTLIER==0,]
df_chr = df_chr[df_chr$DATASET!="GENOME",]
df_chr = df_chr[df_chr$MOSTA.1MODEL1!="",]
lens = unique(df_chr[,c("DATASET","CHROM_LENGTH")])
lens = lens[order(lens$DATASET),]
tab=table(df_chr[,c("DATASET","MOSTA.1MODEL1")])
plot(log10(lens$CHROM_LENGTH),tab[,1]/rowSums(tab),ylab="PROP IBA")
mod=lm(tab[,1]/rowSums(tab)~log10(lens$CHROM_LENGTH))
abline(mod,col="red"); summary(mod)
plot(log10(lens$CHROM_LENGTH),tab[,2]/rowSums(tab),ylab="PROP IBB")
mod=lm(tab[,2]/rowSums(tab)~log10(lens$CHROM_LENGTH))
abline(mod,col="red"); summary(mod)
plot(log10(lens$CHROM_LENGTH),tab[,3]/rowSums(tab),ylab="PROP IBD")
mod=lm(tab[,3]/rowSums(tab)~log10(lens$CHROM_LENGTH))
abline(mod,col="red"); summary(mod)
plot(log10(lens$CHROM_LENGTH),tab[,4]/rowSums(tab),ylab="PROP IBE")
mod=lm(tab[,4]/rowSums(tab)~log10(lens$CHROM_LENGTH))
abline(mod,col="red"); summary(mod)
plot(log10(lens$CHROM_LENGTH),tab[,5]/rowSums(tab),ylab="PROP IBH")
mod=lm(tab[,5]/rowSums(tab)~log10(lens$CHROM_LENGTH))
abline(mod,col="red"); summary(mod)
plot(log10(lens$CHROM_LENGTH),tab[,6]/rowSums(tab),ylab="PROP MIXED")
mod=lm(tab[,6]/rowSums(tab)~log10(lens$CHROM_LENGTH))
abline(mod,col="red"); summary(mod)
tab=table(df_chr[,c("DATASET","MOSTA.1MODEL2")])
plot(log10(lens$CHROM_LENGTH),tab[,1]/rowSums(tab),ylab="PROP IBA")
mod=lm(tab[,1]/rowSums(tab)~log10(lens$CHROM_LENGTH))
abline(mod,col="red"); summary(mod)
plot(log10(lens$CHROM_LENGTH),tab[,2]/rowSums(tab),ylab="PROP IBD")
mod=lm(tab[,2]/rowSums(tab)~log10(lens$CHROM_LENGTH))
abline(mod,col="red"); summary(mod)
plot(log10(lens$CHROM_LENGTH),tab[,3]/rowSums(tab),ylab="PROP IBE")
mod=lm(tab[,3]/rowSums(tab)~log10(lens$CHROM_LENGTH))
abline(mod,col="red"); summary(mod)
plot(log10(lens$CHROM_LENGTH),tab[,4]/rowSums(tab),ylab="PROP IBH")
mod=lm(tab[,4]/rowSums(tab)~log10(lens$CHROM_LENGTH))
abline(mod,col="red"); summary(mod)
plot(log10(lens$CHROM_LENGTH),tab[,5]/rowSums(tab),ylab="PROP MIXED")
mod=lm(tab[,5]/rowSums(tab)~log10(lens$CHROM_LENGTH))
abline(mod,col="red"); summary(mod)
tab=table(df_chr[,c("DATASET","MOSTA.1MODEL3")])
plot(log10(lens$CHROM_LENGTH),tab[,1]/rowSums(tab),ylab="PROP IBA")
mod=lm(tab[,1]/rowSums(tab)~log10(lens$CHROM_LENGTH))
abline(mod,col="red"); summary(mod)
plot(log10(lens$CHROM_LENGTH),tab[,2]/rowSums(tab),ylab="PROP IBE")
mod=lm(tab[,2]/rowSums(tab)~log10(lens$CHROM_LENGTH))
abline(mod,col="red"); summary(mod)
plot(log10(lens$CHROM_LENGTH),tab[,3]/rowSums(tab),ylab="PROP IBH")
mod=lm(tab[,3]/rowSums(tab)~log10(lens$CHROM_LENGTH))
abline(mod,col="red"); summary(mod)
plot(log10(lens$CHROM_LENGTH),tab[,4]/rowSums(tab),ylab="PROP MIXED")
mod=lm(tab[,4]/rowSums(tab)~log10(lens$CHROM_LENGTH))
abline(mod,col="red"); summary(mod)
model1=aov(log(df_chr$CHROM_LENGTH[df_chr$MOSTA.1MODEL1!=""]) ~ df_chr$MOSTA.1MODEL1[df_chr$MOSTA.1MODEL1!=""])
model2=aov(log(df_chr$CHROM_LENGTH[df_chr$MOSTA.1MODEL2!=""]) ~ df_chr$MOSTA.1MODEL2[df_chr$MOSTA.1MODEL2!=""])
model3=aov(log(df_chr$CHROM_LENGTH[df_chr$MOSTA.1MODEL3!=""]) ~ df_chr$MOSTA.1MODEL3[df_chr$MOSTA.1MODEL3!=""])
summary(model1)
summary(model2)
summary(model3)
## test if % explained variance of IBD models is corr with recomb
df_rec = df[!(is.na(df$MEAN_RECOMB_100)),]
plot(df_rec$MEAN_RECOMB_100[df_rec$MOSTA.1MODEL1=="IBD"],df_rec$MAX1[df_rec$MOSTA.1MODEL1=="IBD"])
mod=lm(df_rec$MAX1[df_rec$MOSTA.1MODEL1=="IBD"]~df_rec$MEAN_RECOMB_100[df_rec$MOSTA.1MODEL1=="IBD"])
abline(mod,col="red")
summary(mod)
plot(df_rec$MEAN_RECOMB_100[df_rec$MOSTA.1MODEL1=="IBA"],df_rec$MAX1[df_rec$MOSTA.1MODEL1=="IBA"])
mod=lm(df_rec$MAX1[df_rec$MOSTA.1MODEL1=="IBA"]~df_rec$MEAN_RECOMB_100[df_rec$MOSTA.1MODEL1=="IBA"])
abline(mod,col="red")
summary(mod)
plot(df_rec$MEAN_RECOMB_100[df_rec$MOSTA.1MODEL2=="IBD"],df_rec$MAX1[df_rec$MOSTA.1MODEL2=="IBD"])
mod=lm(df_rec$MAX1[df_rec$MOSTA.1MODEL2=="IBD"]~df_rec$MEAN_RECOMB_100[df_rec$MOSTA.1MODEL2=="IBD"])
abline(mod,col="red")
summary(mod)
plot(df$MEAN_RECOMB_100,df$MAX1)
plot(df$MEAN_DXY_100,df$MAX1)
plot(df$MEAN_FST_100,df$MAX1)
plot(df$MEAN_FST_100[df$MEAN_FST_100<0.5],df$MAX1[df$MEAN_FST_100<0.5])
plot(df$MEAN_DXY_75,df$MAX1)
plot(df$MEAN_FST_75,df$MAX1)
plot(df$MEAN_DXY_50,df$MAX1)
plot(df$MEAN_FST_50,df$MAX1)
plot(log10(df$CHROM_LENGTH),df$MAX1)
plot(df$PERCENT_MISSING,df$MAX1)
summary(lm(df$MEAN_RECOMB_100~df$MAX1)) 
summary(lm(df$MEAN_DXY_100~df$MAX1))
summary(lm(df$MEAN_FST_100~df$MAX1)) 
summary(lm(df$MEAN_FST_100[df$MEAN_FST_100<0.5]~df$MAX1[df$MEAN_FST_100<0.5])) 
summary(lm(df$MEAN_DXY_75~df$MAX1))
summary(lm(df$MEAN_FST_75~df$MAX1)) 
summary(lm(df$MEAN_DXY_50~df$MAX1))
summary(lm(df$MEAN_FST_50~df$MAX1))
summary(lm(log10(df$CHROM_LENGTH)~df$MAX1))
summary(lm(df$PERCENT_MISSING~df$MAX1)) 
cor(df[,c("MAX1","MEAN_RECOMB_100","MEAN_DXY_100","MEAN_FST_100",
          "MEAN_DXY_75","MEAN_FST_75","MEAN_DXY_50","MEAN_FST_50",
          "CHROM_LENGTH","PERCENT_MISSING")],use="pairwise.complete.obs")[1,]
head(df)
mod1=aov(df$HZAR_width[df$MOSTA.1MODEL1!=""]~
           df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""]
)
summary(mod1)
mod1=glm(df$HZAR_width[df$MOSTA.1MODEL1!=""]~
           df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""]+df$SPECIES[df$MOSTA.1MODEL1!=""]
)
summary(mod1)
car::Anova(mod1)
