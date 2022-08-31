df=read.table("/Users/kprovost/Dropbox (AMNH)/locations_morphology_genome_dissertation.csv",sep=",",header=T)
plot(as.numeric(df$LONG),as.numeric(df$LAT),
pch=as.numeric(as.factor(df$TYPE))-1,
col=c("black","red","blue")[as.numeric(as.factor(df$TYPE))])
library(raster)
shp1=shapefile("/Users/kprovost/Documents/OneDrive - The Ohio State University/Environment/Environmental_Layers_Dissertation/cb_2016_us_state_500k/cb_2016_us_state_500k.shp")
shp2=shapefile("/Users/kprovost/Documents/OneDrive - The Ohio State University/Environment/Environmental_Layers_Dissertation/mexstates/mexstates.shp")
for(spp in sort(unique(df$SPECIES))){
print(spp)
df_s = df[df$SPECIES==spp,]
png(paste(spp,"_localities_morph_genome.png",sep=""))
plot(df$LONG,df$LAT,type="n",main=spp)
plot(shp1,add=T)
plot(shp2,add=T)
points(df_s$LONG,df_s$LAT,
pch=as.numeric(as.factor(df_s$TYPE))-1,
col=c("black","red","blue")[as.numeric(as.factor(df_s$TYPE))])
legend("bottomleft",pch=0:2,col=c("black","red","blue"),
legend=c("Morph","Both","Seq"))
dev.off()
}
