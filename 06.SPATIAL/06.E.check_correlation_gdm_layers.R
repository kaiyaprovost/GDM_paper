library(raster)
path="~/Distances/"
x_1 = list.files(path=path,
                 pattern="AGGFACT-1.MORPH-AND-GENE-DISTANCES.txt.converted.txtREDUCED.csv",recursive = T,full.names = T)
x_2 = list.files(path=path,
                 pattern="_distancematrix_env.csv",recursive = T,full.names = T)
x_3 = list.files(path=path,
                 pattern="_distancematrix_geog.csv",recursive = T,full.names = T)
x_4 = list.files(path=path,
                 pattern="_distancematrix_str.csv",recursive = T,full.names = T)
x=c(x_1,x_2,x_3,x_4)
spp_list=c("BELLII","BILINEATA","BRUNNEICAPILLUS","CRISSALE","CURVIROSTRE",
           "FLAVICEPS","FUSCA","MELANURA","NITENS","SINUATUS")
library(spatialEco)
library(gstat)                                         
library(sp)
library(RColorBrewer)
for(spp in sample(spp_list[c(1:10)])){
  x_spp = x[grepl(spp,x)]
  x_spp_short = sapply(x_spp,FUN=function(x){
    splits=strsplit(x,split="/")[[1]]
    lastsplit = splits[length(splits)]
  }
  )
  x_spp_short=make.unique(x_spp_short)
  names(x_spp_short)=NULL
  x_spp_short  = gsub("distancematrix_","",x_spp_short)
  x_spp_short  = gsub(".csv","",x_spp_short)
  x_spp_short  = gsub("AddedPredThin_BestModel_.*addedL.*yers","Env",x_spp_short)
  x_spp_short  = gsub(".rescaled.ascAGGFACT-1.MORPH-AND-GENE-DISTANCES.txt.converted.txtREDUCED","",x_spp_short)
  x_spp_short  = gsub("_clipped.asc","",x_spp_short)
  x_spp_short  = gsub(".asc","",x_spp_short)
  x_spp_short  = gsub("AGGFACT-1.MORPH-AND-GENE-DISTANCES.txt.converted.txtREDUCED","",x_spp_short)
  x_list = lapply(x_spp,FUN=function(x){
    print(x)
    y=read.table(x,sep=",",header=T,row.names = NULL,check.names=T)
    head(y)
    if(ncol(y)<=1){
      y=read.table(x,sep="\t",header=T,row.names = NULL,check.names=T)
    }
    if(ncol(y)<=1){
      y=read.table(x,sep=" ",header=T,row.names = NULL,check.names=T)
    }
    y[,1] = make.unique(y[,1])
    rownames(y) = y[,1]
    y=y[,-1]
    y=y[sort(rownames(y)),sort(colnames(y))]
    return(y)
  })
  names(x_list) = x_spp_short
  names=lapply(x_list,FUN=function(x){
    y=unique(c(colnames(x),rownames(x)))
    return(y)
  })
  list_names=sort(unique(unlist(names)))
  new_x_list=lapply(x_list,FUN=function(x){
    miss_names_row= list_names[which(!(list_names %in% unique(c(rownames(x)))))]
    new_df_row = as.data.frame(matrix(data=NA,nrow=length(miss_names_row),ncol=ncol(x)))
    rownames(new_df_row) = miss_names_row
    colnames(new_df_row) = colnames(x)
    new_x = rbind(x,new_df_row)
    miss_names_col= list_names[which(!(list_names %in% unique(c(colnames(new_x)))))]
    new_df_col = as.data.frame(matrix(data=NA,ncol=length(miss_names_col),nrow=nrow(new_x)))
    rownames(new_df_col) = rownames(new_x)
    colnames(new_df_col) = miss_names_col
    new_x = cbind(new_x,new_df_col)
    new_x=new_x[sort(rownames(new_x)),sort(colnames(new_x))]
    new_x=as.matrix(new_x)
    if(sum(colnames(new_x)!=rownames(new_x)) <= 0){
      diag(new_x)[is.na(diag(new_x))] = 0
    } 
    x_scale=scales::rescale(new_x)
    print(dim(x_scale))
    return(x_scale)
  })
  names(new_x_list) = x_spp_short
  new_x_raster = lapply(new_x_list,FUN=function(x){
    x_ras = raster(as.matrix(x))
  })
  GEO = new_x_raster[grepl("geog",x_spp_short)]
  ENV = new_x_raster[grepl("env",x_spp_short)]
  ABUN = new_x_raster[grepl("Idw_abundance",x_spp_short)]
  STR = new_x_raster[grepl("str",x_spp_short)]
  PRES = new_x_raster[grepl("worldclim",x_spp_short)]
  MID = new_x_raster[grepl("MID",x_spp_short)]
  LGM = new_x_raster[grepl("LGM",x_spp_short)]
  x_stack = stack(new_x_raster)
  GEO_stack = mean(stack(GEO),na.rm=T)
  ENV_stack = mean(stack(ENV),na.rm=T)
  ABUN_stack = mean(stack(ABUN),na.rm=T)
  STR_stack = mean(stack(STR),na.rm=T)
  PRES_stack = mean(stack(PRES),na.rm=T)
  MID_stack = mean(stack(MID),na.rm=T)
  LGM_stack = mean(stack(LGM),na.rm=T)
  x_stack_mean = stack(GEO_stack,ABUN_stack,STR_stack,
                       PRES_stack,LGM_stack,
                       ENV_stack)
  names(x_stack_mean) = c("IBD","IBA","IBB","IBE-res","IBH","IBE-dist")
  not_na = (!(is.na(values(ENV_stack))) & !(is.na(values(GEO_stack))) & !(is.na(values(ABUN_stack)))
            & !(is.na(values(STR_stack))) & !(is.na(values(PRES_stack))) & !(is.na(values(MID_stack))) 
            & !(is.na(values(LGM_stack))))
  RAS_na = GEO_stack
  RAS_na[not_na] = 1
  RAS_na[!(not_na)] = NA
  ENV_res = ENV_stack; values(ENV_res)=NA; values(ENV_res)[!(is.na(values(ENV_stack))) & !(is.na(values(GEO_stack)))] = lm(values(ENV_stack)~values(GEO_stack))$residuals
  ABUN_res = ABUN_stack; values(ABUN_res)=NA; values(ABUN_res)[!(is.na(values(ABUN_stack)))& !(is.na(values(GEO_stack)))] = lm(values(ABUN_stack)~values(GEO_stack))$residuals
  STR_res = STR_stack; values(STR_res)=NA; values(STR_res)[!(is.na(values(STR_stack)))& !(is.na(values(GEO_stack)))] = lm(values(STR_stack)~values(GEO_stack))$residuals
  PRES_res = PRES_stack; values(PRES_res)=NA; values(PRES_res)[!(is.na(values(PRES_stack)))& !(is.na(values(GEO_stack)))] = lm(values(PRES_stack)~values(GEO_stack))$residuals
  MID_res = MID_stack; values(MID_res)=NA; values(MID_res)[!(is.na(values(MID_stack)))& !(is.na(values(GEO_stack)))] = lm(values(MID_stack)~values(GEO_stack))$residuals
  LGM_res = LGM_stack; values(LGM_res)=NA; values(LGM_res)[!(is.na(values(LGM_stack)))& !(is.na(values(GEO_stack)))] = lm(values(LGM_stack)~values(GEO_stack))$residuals
  ENV_res[is.na(RAS_na)] =NA
  ABUN_res[is.na(RAS_na)] =NA
  STR_res[is.na(RAS_na)] =NA
  PRES_res[is.na(RAS_na)] =NA
  MID_res[is.na(RAS_na)] =NA
  LGM_res[is.na(RAS_na)] =NA
  ENV_res_m = as.matrix(ENV_res)
  PRES_res_m = as.matrix(PRES_res)
  ABUN_res_m = as.matrix(ABUN_res)
  STR_res_m = as.matrix(STR_res)
  MID_res_m = as.matrix(MID_res)
  LGM_res_m = as.matrix(LGM_res)
  keep=which(colSums(!(is.na(ENV_res_m)))>1)
  ENV_res_m=ENV_res_m[keep,keep]
  LGM_res_m=LGM_res_m[keep,keep]
  PRES_res_m=PRES_res_m[keep,keep]
  ABUN_res_m=ABUN_res_m[keep,keep]
  STR_res_m=STR_res_m[keep,keep]
  MID_res_m=MID_res_m[keep,keep]
  ENV_res=raster(ENV_res_m)
  MID_res=raster(MID_res_m)
  ABUN_res=raster(ABUN_res_m)
  PRES_res=raster(PRES_res_m)
  STR_res=raster(STR_res_m)
  LGM_res=raster(LGM_res_m)
  x_stack_res = stack(ENV_res,ABUN_res,STR_res,PRES_res,LGM_res)
  names(x_stack_res) = c("IBE-dist","IBA","IBB","IBE-res","IBH")
  x_stats_res=layerStats(x_stack_res,"pearson",na.rm=T,as.sample=F)
  cor = x_stats_res$`pearson correlation coefficient`
  png(paste(spp,"_gdm_corr_residuals.png",sep=""))
  corrplot::corrplot(cor,method="color")
  corrplot::corrplot(cor,method="number",add=T,col=c("black"),bg=rgb(0,0,0,0),
                     cl.pos="n")
  dev.off()
  x_stats_means = layerStats(x_stack_mean,"pearson",na.rm=T,as.sample=F)
  png(paste(spp,"_cormat_plots_per_thing.png",sep=""),width=12,height=12,units="in",res=1200)
  par(mfrow=c(4,4),mar=c(2,2,0,0))
  cormat = matrix(nrow=nlayers(x_stack_mean),ncol=nlayers(x_stack_mean))
  for(i in 1: nlayers(x_stack_mean)){
    for(j in 1:nlayers(x_stack_mean)){
      r_i=(x_stack_mean[[i]])
      r_j=(x_stack_mean[[j]])
      plot(values(r_j),values(r_i),ylab=names(x_stack_mean)[j],xlab=names(x_stack_mean)[i])
      mod=lm(values(r_i)~values(r_j))
      cor_val=sqrt(summary(mod)$adj.r.squared)
      if(mod$coefficients[2]<0) {cor_val=cor_val*-1}
      cormat[i,j] = cor_val
    }
  }
  dev.off()
  colnames(cormat) = names(x_stack_mean)
  rownames(cormat) = names(x_stack_mean)
  png(paste(spp,"_gdm_raw_calc_corr_means.png",sep=""))
  corrplot::corrplot(cormat,method="color")
  corrplot::corrplot(cormat,method="number",add=T,col=c("black"),bg=rgb(0,0,0,0),
                     cl.pos="n")
  dev.off()
  cor = x_stats_means$`pearson correlation coefficient`
  cor_flat = cor
  cor_flat[cor_flat>=1]=1
  cor_flat[cor_flat<=-1]=-1
  png(paste(spp,"_gdm_corr_means_flat_pearsons.png",sep=""))
  corrplot::corrplot(cor_flat,method="color")
  corrplot::corrplot(cor_flat,method="number",add=T,col=c("black"),bg=rgb(0,0,0,0),
                     cl.pos="n")
  dev.off()
  names(x_stack) = x_spp_short
  x_stats=layerStats(x_stack,"pearson",na.rm=T,as.sample=F)
  cor = x_stats$`pearson correlation coefficient`
  cor_flat = cor
  cor_flat[cor_flat>=1]=1
  cor_flat[cor_flat<=-1]=-1
  colnames(cor)
  png(paste(spp,"_gdm_corr.png",sep=""))
  corrplot::corrplot(cor_flat,
                     is.corr=T,method="color",
                     col=colorRampPalette(brewer.pal(11,"RdBu"))(200),
                     order="hclust")
  dev.off()
}
