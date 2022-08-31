data="bivariate_gdm_results.csv"
df = read.csv(data,header=T)
df = df[df$GDMSUBSET==1 & !(is.na(df$GDMSUBSET)),]
a1=aggregate(df$proportion_fst_outliers_under_5sd~df$SPECIES+df$GDMCOLOR,FUN=mean)
a2=aggregate(df$proportion_fst_outliers_over_5sd~df$SPECIES+df$GDMCOLOR,FUN=mean)
a3=merge(a1,a2)
plot((a3$`df$proportion_fst_outliers_under_5sd`*100),(a3$`df$proportion_fst_outliers_over_5sd`*100),
xlim=c(0,2),ylim=c(0,2))
df$MOSTA.1MODEL1[df$MOSTA.1MODEL1==""] = NA
df$MOSTA.1MODEL2[df$MOSTA.1MODEL2==""] = NA
df$MOSTA.1MODEL3[df$MOSTA.1MODEL3==""] = NA
colors = RColorBrewer::brewer.pal(12,"Dark2")
model1=c(colors[c(1,2,3,7,4)],"grey")
model2=c(colors[c(1,3,7,4)],"grey")
model3=c(colors[c(1,7,4)],"grey")
pdf("fst_outliers_5sd_vs_gdm_model.pdf")
par(mfrow=c(2,2))
## Univariate
boxplot(df$number_fst_outliers_under_5sd[!(is.na(df$MOSTA.1MODEL1))]~df$MOSTA.1MODEL1[!(is.na(df$MOSTA.1MODEL1))],col=model1,ylab="Number FST Outliers Under 5 SD",xlab="Univariate Models",main="IBH-IBA, IBH-IBB")
box(col="red",lwd=2)
boxplot(df$number_fst_outliers_over_5sd[!(is.na(df$MOSTA.1MODEL1))]~df$MOSTA.1MODEL1[!(is.na(df$MOSTA.1MODEL1))],col=model1,ylab="Number FST Outliers Over 5 SD",xlab="Univariate Models",main="n.s.")
boxplot(df$proportion_fst_outliers_under_5sd[!(is.na(df$MOSTA.1MODEL1))]~df$MOSTA.1MODEL1[!(is.na(df$MOSTA.1MODEL1))],col=model1,ylab="Proportion FST Outliers Under 5 SD",xlab="Univariate Models",main="IBH-IBB, IBH-IBE")
box(col="red",lwd=2,lty=2)
boxplot(df$proportion_fst_outliers_over_5sd[!(is.na(df$MOSTA.1MODEL1))]~df$MOSTA.1MODEL1[!(is.na(df$MOSTA.1MODEL1))],col=model1,ylab="Proportion FST Outliers Over 5 SD",xlab="Univariate Models",main="n.s.")
boxplot(df$number_fst_outliers_under_5sd[!(is.na(df$MOSTA.1MODEL2))]~df$MOSTA.1MODEL2[!(is.na(df$MOSTA.1MODEL2))],col=model2,ylab="Number FST Outliers Under 5 SD",xlab="Bivariate Models",main="n.s.")
boxplot(df$number_fst_outliers_over_5sd[!(is.na(df$MOSTA.1MODEL2))]~df$MOSTA.1MODEL2[!(is.na(df$MOSTA.1MODEL2))],col=model2,ylab="Number FST Outliers Over 5 SD",xlab="Bivariate Models",main="n.s.")
boxplot(df$proportion_fst_outliers_under_5sd[!(is.na(df$MOSTA.1MODEL2))]~df$MOSTA.1MODEL2[!(is.na(df$MOSTA.1MODEL2))],col=model2,ylab="Proportion FST Outliers Under 5 SD",xlab="Bivariate Models",main="n.s.")
boxplot(df$proportion_fst_outliers_over_5sd[!(is.na(df$MOSTA.1MODEL2))]~df$MOSTA.1MODEL2[!(is.na(df$MOSTA.1MODEL2))],col=model2,ylab="Proportion FST Outliers Over 5 SD",xlab="Bivariate Models",main="n.s.")
boxplot(df$number_fst_outliers_under_5sd[!(is.na(df$MOSTA.1MODEL3))]~df$MOSTA.1MODEL3[!(is.na(df$MOSTA.1MODEL3))],col=model3,ylab="Number FST Outliers Under 5 SD",xlab="Multivariate Models",main="n.s.")
boxplot(df$number_fst_outliers_over_5sd[!(is.na(df$MOSTA.1MODEL3))]~df$MOSTA.1MODEL3[!(is.na(df$MOSTA.1MODEL3))],col=model3,ylab="Number FST Outliers Over 5 SD",xlab="Multivariate Models",main="n.s.")
boxplot(df$proportion_fst_outliers_under_5sd[!(is.na(df$MOSTA.1MODEL3))]~df$MOSTA.1MODEL3[!(is.na(df$MOSTA.1MODEL3))],col=model3,ylab="Proportion FST Outliers Under 5 SD",xlab="Multivariate Models",main="n.s.")
boxplot(df$proportion_fst_outliers_over_5sd[!(is.na(df$MOSTA.1MODEL3))]~df$MOSTA.1MODEL3[!(is.na(df$MOSTA.1MODEL3))],col=model3,ylab="Proportion FST Outliers Over 5 SD",xlab="Multivariate Models",main="n.s.")
dev.off()
boxplot(df$number_fst_outliers_under_5sd~df$SPECIES)
boxplot(df$number_fst_outliers_over_5sd~df$SPECIES)
boxplot(df$proportion_fst_outliers_under_5sd~df$SPECIES)
boxplot(df$proportion_fst_outliers_over_5sd~df$SPECIES)
## Multivariate -- 
## number -5 SD
modelA=aov(df$number_fst_outliers_under_5sd[!(is.na(df$MOSTA.1MODEL3))]~df$MOSTA.1MODEL3[!(is.na(df$MOSTA.1MODEL3))])
summary(modelA) 
TukeyHSD(modelA) 
modelA=aov(df$number_fst_outliers_under_5sd[!(is.na(df$MOSTA.1MODEL3)) & df$SPECIES!="NITENS"]~df$MOSTA.1MODEL3[!(is.na(df$MOSTA.1MODEL3)) & df$SPECIES!="NITENS"])
summary(modelA) 
TukeyHSD(modelA) 
## number +5 SD
modelA=aov(df$number_fst_outliers_over_5sd[!(is.na(df$MOSTA.1MODEL3))]~df$MOSTA.1MODEL3[!(is.na(df$MOSTA.1MODEL3))])
summary(modelA)
TukeyHSD(modelA) 
## prop -5 SD
modelA=aov(df$proportion_fst_outliers_under_5sd[!(is.na(df$MOSTA.1MODEL3))]~df$MOSTA.1MODEL3[!(is.na(df$MOSTA.1MODEL3))])
summary(modelA)
TukeyHSD(modelA) 
## prop +5 SD
modelA=aov(df$proportion_fst_outliers_over_5sd[!(is.na(df$MOSTA.1MODEL3))]~df$MOSTA.1MODEL3[!(is.na(df$MOSTA.1MODEL3))])
summary(modelA)
TukeyHSD(modelA) 
