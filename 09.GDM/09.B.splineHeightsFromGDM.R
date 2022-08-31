path="~/Distances/"
setwd(path)
outfile = "extracting_best_column_ngsdist_lostruct_28July2022.alltogether"
files = list.files(path=path,
                   pattern="modelparameters",recursive=F,full.names = T)
files = files[!(grepl(".txt.gz$",files))]
newdf=NULL
for (i in 1:length(files)) {
  print(paste(i,"/",length(files)))
  file = files[i]
  base = tools::file_path_sans_ext(basename(file))
  types=strsplit(base,"_")[[1]]
  spp=types[1]
  variables=types[2]
  dataset=paste(types[4:length(types)],sep="-",collapse="-")
  txt=readLines(file,n=20)
  null=txt[grepl("NULL Deviance",txt)]
  dev= txt[grepl("GDM Deviance",txt)]
  exp= txt[grepl("Percent Deviance Explained",txt)]
  null=as.numeric(strsplit(null,":")[[1]][2])
  dev= as.numeric(strsplit(dev, ":")[[1]][2])
  exp= as.numeric(strsplit(exp, ":")[[1]][2])
  linetoadd = cbind("SPECIES"=spp,"DATASET"=dataset,"MODEL"=variables,"NULL"=null,"DEV"=dev,"EXP"=exp)
  if(is.null(newdf)){
    newdf = linetoadd
    write.table(linetoadd,outfile,append=T,row.names=F,col.names = T,sep="\t",quote=F)
  } else {
    newdf = rbind(newdf,linetoadd)
    write.table(linetoadd,outfile,append=T,row.names=F,col.names = F,sep="\t",quote=F)
  }
  R.utils::gzip(file,overwrite=T)
}
newdf=unique(newdf)
write.table(newdf,outfile,row.names = F,append = T,sep="\t",quote=F)
newdf = as.data.frame(newdf)
models_list = sort(unique(newdf$MODEL))
dataset_list = sort(unique(newdf$DATASET))
species_list = sort(unique(newdf$SPECIES))
dataset_species = expand.grid(dataset_list, species_list, stringsAsFactors = FALSE)
pivot_df = NULL
for(model_i in models_list){
  print(model_i)
  newdf_mod = newdf[newdf$MODEL==model_i,]
  newdf_mod = newdf_mod[,c("SPECIES","DATASET","EXP")]
  colnames(newdf_mod) = c("SPECIES","DATASET",model_i)
  if(is.null(pivot_df)){
    pivot_df = newdf_mod
  } else {
    pivot_df = merge(pivot_df,newdf_mod,by=c("SPECIES","DATASET"),all=T)
  }
}
write.table(pivot_df,paste(outfile,".PIVOT.txt",sep=""),row.names = F,append = T,sep="\t",quote=F)
## get spline heights 
files = list.files(path=path,
                   pattern="splineheights.*txt",recursive=F,full.names = T)
files=files[grepl("SUBSET",files)]
fulldf = data.frame()
for(i in 1:length(files)) {
  print(paste(i,"/",length(files)))
  file = files[i]
  df = read.table(file,header=T,sep=",",colClasses = c("character","numeric"),row.names = 1)
  df = t(df)
  df=as.data.frame(df)
  filesplit=strsplit(basename(file),"_")[[1]]
  species=filesplit[1]
  params=filesplit[2]
  chrom=tools::file_path_sans_ext(filesplit[4])
  df$SPP = species
  df$PARAMS = params
  df$CHROM = chrom
  write.table(df,"~/splineheights_28July2022.txt",append=T)
  fulldf=plyr::rbind.fill(fulldf,df)
}
okay=c("ABUN","ENV","IBD","LGM","PRES","STR",
       "STR-ABUN","STR-ENV","STR-IBD","STR-LGM","STR-PRES",
       "STR-IBD-ABUN","STR-IBD-ENV","STR-IBD-LGM","STR-IBD-PRES")
fulldf=fulldf[,c("SPP","PARAMS","CHROM","ABUN","ENV","IBD","LGM","PRES","STR")]
fulldf=unique(fulldf)
fulldf=fulldf[rowSums(is.na(fulldf))!=6,]
fulldf=fulldf[!(grepl("PC",fulldf$PARAMS)),]
fulldf=fulldf[fulldf$PARAMS %in% okay,]
write.table(fulldf,"splineheights.txt",append=F)
fulldfna = fulldf
fulldfna[is.na(fulldfna)] = 0
fulldfna$sum = rowSums(fulldfna[,c(4:9)])
fulldf_relative = fulldfna
for(row in 1:nrow(fulldfna)) {
  print(row)
  fulldf_relative[row,4:9] = (fulldfna[row,4:9] / fulldfna[row,10])
}
abun = aggregate(fulldf_relative$ABUN~fulldf_relative$PARAMS,FUN=function(x){mean(x,na.rm=T)}); colnames(abun) = c("PARAMS","ABUN")
env = aggregate(fulldf_relative$ENV~fulldf_relative$PARAMS,FUN=function(x){mean(x,na.rm=T)}); colnames(env) = c("PARAMS","ENV")
ibd = aggregate(fulldf_relative$IBD~fulldf_relative$PARAMS,FUN=function(x){mean(x,na.rm=T)}); colnames(ibd) = c("PARAMS","IBD")
lgm = aggregate(fulldf_relative$LGM~fulldf_relative$PARAMS,FUN=function(x){mean(x,na.rm=T)}); colnames(lgm) = c("PARAMS","LGM")
pres = aggregate(fulldf_relative$PRES~fulldf_relative$PARAMS,FUN=function(x){mean(x,na.rm=T)}); colnames(pres) = c("PARAMS","PRES")
str = aggregate(fulldf_relative$STR~fulldf_relative$PARAMS,FUN=function(x){mean(x,na.rm=T)}); colnames(str) = c("PARAMS","STR")
forbarplot = merge(merge(abun,env,all=T),merge(merge(ibd,lgm,all=T),merge(pres,str,all=T),all=T),all=T)
barplot(forbarplot)
par(mar=c(8,6,2,2))
boxplot(fulldf$ABUN~fulldf$PARAMS,las=2,xlab="")
