files=list.files(path="~/abundance_asciis/",
                 pattern="Idw_abundance_clipped.asc$",full.names = T)
for (file in files){
  print(file)
  asc = raster(file)
  print(mean(values(asc),na.rm=T))
  print(sd(values(asc),na.rm=T))
}
abun_stack = stack(files)
abun_sum = calc(abun_stack,fun=sum,na.rm=T)
abun_sd = calc(abun_stack,fun=sd,na.rm=T)
abun_norm = abun_stack
for(i in 1:10){
  this = abun_norm[[i]]
  this[this<1]=0
  this[this>=1]=1
  abun_norm[[i]]=this
}
abun_norm_sum = calc(abun_norm,fun=sum,na.rm=T)
yellow=colorRampPalette(brewer.pal(9,"YlOrRd")[3:9])
abun_sum[is.na(abun_stack[[1]])] = NA
abun_norm_sum[is.na(abun_stack[[1]])] = NA
writeRaster(abun_norm_sum,"AbundanceStackRaster.asc",format="ascii")
pdf("abundance_graphs_across_species.pdf",height=10,width=10)
par(mfrow=c(2,2))
plot(abun_sum,main="abundance summed",col=yellow(11))
plot(log10(abun_sum),main="log abundance summed",col=yellow(11))
plot(abun_sd,main="abundance stdev",col=yellow(11))
plot(abun_norm_sum,main="binary abundance summed",col=yellow(11))
dev.off()
red=colorRampPalette(brewer.pal(12,"RdYlBu"))
png("binary_abundance_sum.png")
plot(abun_norm_sum,main="binary abundance summed",col=red(11))
dev.off()
overall="means_and_stdevs_across_species.txt"
df = read.table(overall,header=T,sep="\t",row.names = 1)
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
pdf("all_species_barplot_mean_sd_variables.pdf",height=10,width=10)
par(mfrow=c(3,3),mar=c(4,4,1,0.5))
colnames(df)[seq(1,ncol(df),2)] = c("Mean Normalized RF Distance","Mean Cline Width (degrees)",
                                    "Mean Cline Center","Mean Recombination Rate (x 10^-9)",
                                    "Mean FST","Mean DXY","Mean Missing Data","Mean Abundance")
colnames(df)[seq(2,ncol(df),2)] = c("SD Normalized RF Distance","SD Cline Width (degrees)",
                                    "SD Cline Center","SD Recombination Rate (x 10^-9)",
                                    "SD FST","SD DXY","SD Missing Data","SD Abundance")
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
colors=c("black","red","orange","goldenrod","green",
         "blue","cyan","magenta","brown","grey")
short = sapply(rownames(df),FUN=function(x){strsplit(x," ")[[1]][2]})
short = substr(short,1,3)
for(i in seq(1,ncol(df),2)) {
  maximum=max(df[,i])+max(df[,i+1])
  x=barplot(df[,i],names=short,las=2,
            ylim=c(0,maximum*1.01),
            ylab=colnames(df)[i],
            col=col2(10))
  box()
  error.bar(x,df[,i],df[,(i+1)])
}
corrplot(cor(df),method="color",order="hclust",diag=T,type="lower",tl.pos="ld",
         tl.cex=0.7,tl.col="black")
dev.off()
png("cline_width_vs_cline_center_with_dots.png")
plot(df[,3],df[,5],col="black",pch=22,bg=colors,
     ylim=c(1,18),xlim=c(3,17),
     ylab=colnames(df)[5],
     xlab=colnames(df)[3])
segments(x0=df[,3]-df[,4],y0=df[,5],
         x1=df[,3]+df[,4],y1=df[,5],col="black",lty=3)
segments(x0=df[,3],y0=df[,5]-df[,6],
         x1=df[,3],y1=df[,5]+df[,6],col="black",lty=3)
legend("topleft",fill=colors,
       legend=short)
points(df[,3],df[,5],col="black",pch=22,bg=colors,
       ylim=c(1,18),xlim=c(3,17),
       ylab=colnames(df)[5],
       xlab=colnames(df)[3],cex=2)
dev.off()
