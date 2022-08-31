## DROP THE MIX IN THE UNI AND BIV
rm(list=ls())
setwd("~")
date=format(Sys.time(), "%d%b%Y")
do_process=F
do_combine=F
do_ngs=F
do_mix=F
do_chi=F
## process from the raw files
if (do_process==T) {
  outfile = paste("extracting_best_column_ngsdist_lostruct_",date,".alltogether",sep="")
  files = list.files(path="~/Distances/",
                     pattern="variable_deviance",recursive=F,full.names = T)
  files = files[grepl("SUBSET",files)]
  newdf = c()
  for (file in files) {
    print(file)
    df = read.csv(file,row.names = 1)
    if(sum(colSums(!is.na(df))) != 0) {
      rightcol=max(which(df[1,] == min(df[1,],na.rm=T)))
      values=(df[,names(df)[rightcol]])
      names(values) = (rownames(df))
      if(rightcol==ncol(df)) {
        toadd = "_UNI"
        rightcol=ncol(df)-1
      } else {
        toadd = ""
      }
      rnam=c("PC1T","PC2T","PC3T","PC1P","PC2P","PC3P",
             "IBD","LGM","PRES","ENV","MID","STR","ABUN")
      otherdf1 = read.csv(sub("deviance","importance",file),row.names=1)
      otherdf2 = read.csv(sub("deviance","significance",file),row.names=1)
      otherdf3 = read.csv(sub("deviance","permutations",file),row.names=1)
      if(nrow(otherdf1) == length(rnam)) {
        rownames(otherdf1) = rnam
        rownames(otherdf2) = rnam
        rownames(otherdf2) = rnam
      } else {
        rownames(otherdf1) = make.unique(as.character(rownames(otherdf1)))
        rownames(otherdf2) = make.unique(as.character(rownames(otherdf2)))
        rownames(otherdf2) = make.unique(as.character(rownames(otherdf3)))
      }
      val1 = as.data.frame(rbind(otherdf1[,colnames(otherdf1)[min(rightcol,ncol(otherdf1),na.rm=T)]]))
      names(val1) = paste(rownames(otherdf1),"imp")
      val2 = as.data.frame(rbind(otherdf1[,colnames(otherdf2)[min(rightcol,ncol(otherdf2),na.rm=T)]]))
      names(val2) = paste(rownames(otherdf2),"sig")
      val3 = as.data.frame(rbind(otherdf3[,colnames(otherdf3)[min(rightcol,ncol(otherdf3),na.rm=T)]]))
      names(val3) = paste(rownames(otherdf3),"per")
      toname = strsplit(file,"/")[[1]]
      folder=toname[10]
      if (is.null(newdf)) {
        newdf = t(as.data.frame(c(paste(toname[length(toname)],toadd,sep=""),
                                  values,val1,val2,val3)))
      } else {
        newrow = t(as.data.frame(c(paste(toname[length(toname)],toadd,sep=""),
                                   values,val1,val2,val3)))
        newdf = rbind(newdf,newrow)
      }
      rownames(newdf) = NULL
    }
  }
  write.table(newdf,outfile,append=F,row.names = F,sep="\t")
  newdf=read.table("~/extracting_best_column.alltogether",
                   sep="\t",header=T)
  ## imp cols = 6:18, sig cols == 19:30, per cols = 31:44
  ## metadata 45:48
  imps = as.data.frame(newdf[,c(1:5,6:18,45:48)])
  sigs = as.data.frame(newdf[,c(1:5,19:30,45:48)])
  imps$VARIABLE = imps$CHROM
  sigs$VARIABLE = sigs$CHROM
  imps = imps[order(imps$VARIABLE),]
  impsnums = imps[,8:18]
  impsnums = as.matrix(impsnums)
  for(i in 1:nrow(impsnums)) {
    for (j in 1:ncol(impsnums)) {
      if(is.na(impsnums[i,j])) {
        impsnums[i,j] = 0
      }
      if(is.infinite(impsnums[i,j])) {
        impsnums[i,j] = 0
      }
      if(impsnums[i,j] < 0) {
        impsnums[i,j] = 0
      }
      impsnums[i,j] = as.numeric(impsnums[i,j])
    }
  }
  for(i in 1:nrow(impsnums)) {
    explain = as.numeric(as.character(imps$`Percent deviance explained`[i]))
    impsnums[i,] = scales::rescale(as.numeric(impsnums[i,]))
    impsnums[i,] = ((as.numeric(impsnums[i,]) / sum(as.numeric(impsnums[i,]))))
  }
  ## get relative explanatory power
  for(i in 1:nrow(impsnums)) {
    impsnums[i,] = as.numeric(as.character(impsnums[i,]))*as.numeric(as.character(imps$`Percent deviance explained`[i]))
  }
  colnames(impsnums) = colnames(imps[,8:18])
  rownames(impsnums) = imps$SPECIES
  impsnums=data.matrix(impsnums)
  rownames(impsnums)=NULL
  colnames(impsnums)=NULL
  dimmy=dim(impsnums)
  impsnums=as.numeric(impsnums)
  impsnums = matrix(impsnums,ncol=11)
  ifremove = rowSums(impsnums,na.rm = T)
  ifremove[is.na(ifremove)] = 0
  ifremove[is.infinite(ifremove)] = 0
  toremove = which(ifremove==0)
  if(length(toremove)!= 0) {
    impsnums = impsnums[-toremove,]
    ## remove if entire row is 0
    env = rowSums(impsnums[,c(1,2,3,4,8)])
    his = rowSums(impsnums[,c(6,7,9,10)])
    ibd = rowSums(impsnums[,c(5,11)])
    combo = cbind(env,his,ibd)
    combo=as.data.frame(combo)
    combo$SPECIES = imps$SPECIES[-toremove]
    combo$VARIABLE = imps$VARIABLE[-toremove]
  } else {
    env = rowSums(impsnums[,c(1,2,3,4,8)])
    his = rowSums(impsnums[,c(6,7,9,10)])
    ibd = rowSums(impsnums[,c(5,11)])
    combo = cbind(env,his,ibd)
    combo=as.data.frame(combo)
    combo$SPECIES = imps$SPECIES
    combo$VARIABLE = imps$VARIABLE
  }
  short = substr(combo$SPECIES,1,3)
  shorta = gsub(".csv","",combo$VARIABLE)
  shorta[shorta=="test"] = "allmorph"
  ids = paste(shorta,short,sep=".")
  for(spp in unique(combo$SPECIES)) {
    print(spp)
    subcombo = combo[combo$SPECIES==spp,]
    plotmax = round((max(rowSums(subcombo[,1:3]))),-1)
    png(paste(spp,"_all_gdm_results_bk.png",sep=""),width=900)
    barplot(t(subcombo[,1:3]),col=c("#d7191c","#ffffbf","#2b83ba"),
            names.arg = subcombo$VARIABLE,ylab="% Variation Explained",
            las=2,cex.names = 1,main=spp,
            ylim=c(0,plotmax+10))
    legend("top",legend=c("ENV","HIS","IBD"),ncol=3,
           fill=c("#d7191c","#ffffbf","#2b83ba"),col=c("#d7191c","#ffffbf","#2b83ba"),
           bty="n",box.col="white",border="white")
    dev.off()
    png(paste(spp,"_all_gdm_results.png",sep=""),width=900)
    par(bg = rgb(0,0,0,0), fg = 'white',
        col="white",col.axis="white",col.lab="white",col.main="white",col.sub="white",
        font=2,xpd=T)
    barplot(t(subcombo[,1:3]),col=c("#d7191c","#ffffbf","#2b83ba"),
            names.arg = subcombo$VARIABLE,ylab="% Variation Explained",
            las=2,cex.names = 1,main=spp,
            ylim=c(0,plotmax+10))
    legend("top",legend=c("ENV","HIS","IBD"),ncol=3,
           fill=c("#d7191c","#ffffbf","#2b83ba"),col=c("#d7191c","#ffffbf","#2b83ba"),
           bty="n",box.col="white",border="white")
    dev.off()
  }
  png("alltogether_gdm_results.png",width=900)
  par(bg = rgb(0,0,0,0), fg = 'white',
      col="white",col.axis="white",col.lab="white",col.main="white",col.sub="white",
      font=2,xpd=T)
  barplot(t(combo[,1:3]),col=c("#d7191c","#ffffbf","#2b83ba"),
          names.arg = ids,ylab="% Variation Explained",las=2,cex.names = 0.5)
  legend("top",legend=c("ENV","HIS","IBD"),ncol=3,
         fill=c("#d7191c","#ffffbf","#2b83ba"),col=c("#d7191c","#ffffbf","#2b83ba"),
         bty="n",box.col="white",border="white")
  dev.off()
  png("morphology_gdm_results_bk.png",width=900)
  barplot(t(combo[,1:3]),
          col=c("#d7191c","#ffffbf","#2b83ba"),
          names.arg = combo$SPECIES,
          las=2,cex.names=0.5,
          ylab="% Variation Explained"#,
  )
  legend("top",legend=c("ENV","HIS","IBD"),ncol=3,
         fill=c("#d7191c","#ffffbf","#2b83ba"),col=c("#d7191c","#ffffbf","#2b83ba"),
         bty="n",box.col="white",border="white")
  dev.off()
  par(bg = rgb(0,0,0,0), fg = 'white',
      col="white",col.axis="white",col.lab="white",col.main="white",col.sub="white",
      font=2,xpd=T)
  barplot(t(impsnums),las=2,col=as.numeric(as.factor(colnames(imps)[8:18])))
  par(bg = rgb(0,0,0,0), fg = 'white',
      col="white",col.axis="white",col.lab="white",col.main="white",col.sub="white",
      font=2,xpd=T)
  barplot(as.numeric(as.character(imps$`Percent deviance explained`)),
          col=as.factor(imps$SPECIES)
  )
}
## generate the ngs data plots
if(do_ngs == T){
  bi = read.csv("~/bivariate_gdm_results_NGSDIST_20July2022.csv",
                sep=",",skip=0,stringsAsFactors = F)
  bi = unique(bi)
  bi_trimmed = bi[,c("DATASET","SPECIES",
                     "MOSTA.1MODEL1","MOSTA.1MODEL2","MOSTA.1MODEL3",
                     "MOSTTOT.1MODEL","MOSTATOT.1MODEL",
                     "MAX1","MAX2","MAX3",
                     "MAXTOT")]
  bi_trimmed = unique(bi_trimmed)
  bi_trimmed$MOSTA.1MODEL1=as.character(bi_trimmed$MOSTA.1MODEL1)
  bi_trimmed$MOSTA.1MODEL1[(bi_trimmed$MOSTA.1MODEL1=="")] = "MIXED"
  bi_trimmed$MOSTA.1MODEL2=as.character(bi_trimmed$MOSTA.1MODEL2)
  bi_trimmed$MOSTA.1MODEL2[(bi_trimmed$MOSTA.1MODEL2=="")] = "MIXED"
  bi_trimmed$MOSTA.1MODEL3=as.character(bi_trimmed$MOSTA.1MODEL3)
  bi_trimmed$MOSTA.1MODEL3[(bi_trimmed$MOSTA.1MODEL3=="")] = "MIXED"
  bi_trimmed$MOSTTOT.1MODEL=as.character(bi_trimmed$MOSTTOT.1MODEL)
  bi_trimmed$MOSTTOT.1MODEL[(bi_trimmed$MOSTTOT.1MODEL=="")] = "MIXED"
  brewcols = RColorBrewer::brewer.pal(8,"Dark2")
  blue=brewcols[3] 
  yellow=brewcols[7]
  red=brewcols[4]
  purple=brewcols[1]
  orange=brewcols[2]
  green=brewcols[5]
  gold=brewcols[6]
  grey=brewcols[8]
  black="#000000"
  for (datanum in 1:4) {
    if(datanum==1){
      dx = bi_trimmed[,c("DATASET","SPECIES","MAX1","MOSTA.1MODEL1")]
    } else if (datanum==2) {
      dx = bi_trimmed[,c("DATASET","SPECIES","MAX2","MOSTA.1MODEL2")]
    } else if (datanum==3) {
      dx = bi_trimmed[,c("DATASET","SPECIES","MAX3","MOSTA.1MODEL3")]
    } else {
      dx = bi_trimmed[,c("DATASET","SPECIES","MAXTOT","MOSTTOT.1MODEL")]
    }
    colnames(dx) = c("DATASET","SPECIES","MAX","MOSTA.1MODEL")
    dx$SCALEDMAX = scales::rescale(dx$MAX,to=c(100,255))
    dx$MAX[is.na(dx$MAX)] = 0
    dx$SCALEDMAX[is.na(dx$SCALEDMAX)] = 0
    maximumvalue=85
    dtt = table(bi_trimmed[,c("DATASET","SPECIES")])
    dtt=dtt[rownames(dtt)!="",colnames(dtt)!=""]
    dttn = dtt
    dtte = dtt
    dx$MOSTA.1MODEL.F = as.factor(dx$MOSTA.1MODEL)
    if(datanum==1){
      colorlist = c(purple,orange,blue,yellow,red,black)
      dx$MOSTA.1MODEL.COLOR = as.numeric(as.factor(dx$MOSTA.1MODEL.F))
      dx$MOSTA.1MODEL.COLOR[dx$MOSTA.1MODEL.F=="IBA"] = 1 ## purple
      dx$MOSTA.1MODEL.COLOR[dx$MOSTA.1MODEL.F=="IBB"] = 2 ## orange
      dx$MOSTA.1MODEL.COLOR[dx$MOSTA.1MODEL.F=="IBD"] = 3 ## blue
      dx$MOSTA.1MODEL.COLOR[dx$MOSTA.1MODEL.F=="IBE"] = 4 ## yellow
      dx$MOSTA.1MODEL.COLOR[dx$MOSTA.1MODEL.F=="IBH"] = 5 ## red
      dx$MOSTA.1MODEL.COLOR[dx$MOSTA.1MODEL.F=="MIXED"] = 6 ## black
    } else if (datanum==2) {
      colorlist = c(purple,blue,yellow,red,black)
      dx$MOSTA.1MODEL.COLOR = as.numeric(as.factor(dx$MOSTA.1MODEL.F))
      dx$MOSTA.1MODEL.COLOR[dx$MOSTA.1MODEL.F=="IBA"] = 1 ## purple
      dx$MOSTA.1MODEL.COLOR[dx$MOSTA.1MODEL.F=="IBD"] = 2 ## blue
      dx$MOSTA.1MODEL.COLOR[dx$MOSTA.1MODEL.F=="IBE"] = 3 ## yellow
      dx$MOSTA.1MODEL.COLOR[dx$MOSTA.1MODEL.F=="IBH"] = 4 ## red
      dx$MOSTA.1MODEL.COLOR[dx$MOSTA.1MODEL.F=="MIXED"] = 5 ## black
    } else if (datanum==3) {
      colorlist = c(purple,yellow,red,black)
      dx$MOSTA.1MODEL.COLOR = as.numeric(as.factor(dx$MOSTA.1MODEL.F))
      dx$MOSTA.1MODEL.COLOR[dx$MOSTA.1MODEL.F=="IBA"] = 1 ## purple
      dx$MOSTA.1MODEL.COLOR[dx$MOSTA.1MODEL.F=="IBE"] = 2 ## yellow
      dx$MOSTA.1MODEL.COLOR[dx$MOSTA.1MODEL.F=="IBH"] = 3 ## red
      dx$MOSTA.1MODEL.COLOR[dx$MOSTA.1MODEL.F=="MIXED"] = 4 ## black
    } else {
      colorlist = c(purple,orange,blue,yellow,red,black)
      dx$MOSTA.1MODEL.COLOR = as.numeric(as.factor(dx$MOSTA.1MODEL.F))
      dx$MOSTA.1MODEL.COLOR[dx$MOSTA.1MODEL.F=="IBA"] = 1 ## purple
      dx$MOSTA.1MODEL.COLOR[dx$MOSTA.1MODEL.F=="IBB"] = 2 ## orange
      dx$MOSTA.1MODEL.COLOR[dx$MOSTA.1MODEL.F=="IBD"] = 3 ## blue
      dx$MOSTA.1MODEL.COLOR[dx$MOSTA.1MODEL.F=="IBE"] = 4 ## yellow
      dx$MOSTA.1MODEL.COLOR[dx$MOSTA.1MODEL.F=="IBH"] = 5 ## red
      dx$MOSTA.1MODEL.COLOR[dx$MOSTA.1MODEL.F=="MIXED"] = 6 ## black
    }
    for(ri in 1:nrow(dtt)) {
      for(ci in 1:ncol(dtt)) {
        print(paste(ri,"-",ci))
        variable = rownames(dtt)[ri]
        species = colnames(dtt)[ci]
        print(paste(variable,species))
        if(variable!="" & species !=""){
          entry = dx[dx$DATASET==variable,]
          entry = entry[entry$SPECIES==species,]
          entry = entry[complete.cases(entry),]
          entry_1 = as.character(as.factor(entry$MOSTA.1MODEL.F))
          entry_2 = as.numeric(entry$MOSTA.1MODEL.COLOR)
          entry_3 = as.numeric(as.character(entry$SCALEDMAX))
          if(length(entry_3) == 0){
            entry_3 = 0
          } else if(is.na(entry_3)){
            entry_3 = 0
          }
          if(length(entry_2) == 0){
            entry_2 = 0
          } else if(is.na(entry_2)){
            entry_2 = 0
          }
          if(length(entry_1) == 0){
            entry_1 = 0
          } else if(is.na(entry_1)){
            entry_1 = 0
          }
          dtt[ri,ci] = entry_1
          dttn[ri,ci] = as.numeric(entry_2)
          dtte[ri,ci] = entry_3
        }
      }
    }
    toordercols_spp = c("NITENS","BILINEATA","BRUNNEICAPILLUS","SINUATUS","CRISSALE","FUSCA","MELANURA","FLAVICEPS","BELLII","CURVIROSTRE")
    toordercols = sapply(toordercols_spp, function(x) {which(colnames(dtt) == x)})
    toorderrows_gen = c("GENOME")
    toorderrows_chr = c(1,"1A","1B",2:4,"4A",5:28,
                        "LG2","LG5","LGE22","MTDNA","Z")
    toorderrows_chr=paste("CHR",toorderrows_chr,sep="")
    toorderrows_mor = c("MORPH","PC1","PC2","PC3")
    toorderrows_fst = c("HIGH","LOW","NORM",
                        "HIGH75","LOW75","NORM75",
                        "HIGH50","LOW50","NORM50")
    toorderrows_lostruct = c("LS0","LS1","LS2","LS3")
    toorderrows = sapply(c(toorderrows_gen,toorderrows_chr,toorderrows_lostruct,toorderrows_fst,toorderrows_mor), function(x) {which(rownames(dtt) == x)})
    dims=dim(dttn)
    dttn=matrix(as.numeric(dttn),nrow=dims[1],ncol=dims[2])
    dtt=dtt[as.numeric(toorderrows),as.numeric(toordercols)]
    dttn=dttn[as.numeric(toorderrows),as.numeric(toordercols)]
    dtte=dtte[as.numeric(toorderrows),as.numeric(toordercols)]
    dttn[dttn==0] = NA
    narows = is.na(rownames(dtt))
    dttn = dttn[!narows,]
    dtte = dtte[!narows,]
    dtt = dtt[!narows,]
    ff = factor(as.matrix(dttn))
    fft = factor(as.matrix(dtt))
    fx=matrix(dttn,ncol=ncol(dtt))
    colnames(fx) = colnames(dtt)
    rownames(fx) = rownames(dtt)
    fxd = as.data.frame(fx)
    #use labels to assign colors
    col<-c("IBA"="purple","IBE"="yellow","IBD"="blue",
           "IBH"="red","MIXED"="black","IBB"="orange")
    imgflip<-function(x) {t(x[nrow(x):1,])}
    png(paste("~/heatmap_gdm_results_",datanum,"_bivariate_ngsdist_",date,".png",sep=""),width=900)
    par(mfrow=c(1,1),xpd=TRUE,mar=c(5,10,6,1))
    image((dttn),col=colorlist,
          xaxt="n",yaxt="n")
    axis(3,at=seq(0,1,length.out=nrow(dtt)),
         labels=(rownames(fx)),las=2)
    axis(2,at=seq(0,1,length.out=ncol(dtt)),
         labels=(colnames(fx)),las=2)
    legend(x = 0,y=-0.1,legend=levels(ff),ncol = 5,
           bty="n",
           fill=colorlist)
    dev.off()
    class(dttn) = "matrix"
    class(dtte) = "matrix"
    fills=colorlist
    col2rgb(fills,alpha=T)
    fillsf = fills[dttn]
    fillsfc = col2rgb(fillsf,alpha=T)
    dttef = matrix(dtte)
    fillsfc[4,] = floor(dttef[,1])
    transparents = rgb(red=fillsfc[1,],green=fillsfc[2,],
                       blue=fillsfc[3,],alpha=fillsfc[4,],
                       maxColorValue = 255)
    transp_f = factor(as.matrix(transparents))
    imgcolors = as.numeric(transp_f)
    fx<-matrix(as.numeric(imgcolors), ncol=ncol(dtt))
    colnames(fx) = colnames(dtt)
    rownames(fx) = rownames(dtt)
    colmat = sapply(colorlist,FUN=function(ramp){
      trans = paste(ramp,"64",sep="") ## 32 for 50, 64 for 100
      solid = paste(ramp,"FF",sep="")
      colfunc <- colorRampPalette(c(solid, trans),alpha=T)
      (matrix(colfunc(maximumvalue), ncol=1))
    })
    ## for each column make a raster and then plot that 
    library(grid)
    colmatf = factor(colmat)
    colcol = as.numeric(colmatf)
    colff<-matrix(as.numeric(colcol), ncol=ncol(colmat))
    colnames(colff) = colnames(colmat)
    rownames(colff) = rownames(colmat)
    if(datanum==1){
      png(paste("heatmap_legend_",date,".png",sep="_"),width=250)
      par(mar=c(4,4,4,4))
      image(colff,col=levels(colmatf),
            xaxt="n",yaxt="n",main="% Explained\n")
      axis(3,at=seq(0,1,length.out=6),
           labels=round(seq(maximumvalue,0,length.out=6)),las=1,
           cex=0.8)
      axis(2,at=seq(0,1,length.out=6),
           labels=c("IBA","IBB","IBD","IBE","IBH","MIX"),las=2)
      dev.off()
    }
    png(paste("~/heatmap_gdm_results_",datanum,"_bivariate_ngs_transparent_",date,".png",sep=""),width=900)
    layout(t(matrix(c(1,3,1,3,1,2,1,2,1,3,1,3),nrow=2)), width = c(6,1),height = c(1,1))
    par(xpd=TRUE,mar=c(1,3,6,1))
    image(fx,col=levels(transp_f),
          xaxt="n",yaxt="n")
    axis(3,at=seq(0,1,length.out=nrow(dtt)),
         labels=(rownames(fx)),las=2)
    axis(2,at=seq(0,1,length.out=ncol(dtt)),
         labels=substr(colnames(fx),1,3),las=2)
    legend(x = 0,y=-0.1,legend=levels(ff),ncol = 5,
           bty="n",
           fill=colorlist)
    image(colff,col=levels(colmatf),
          xaxt="n",yaxt="n",main="% Explained")
    axis(3,at=seq(0,1,length.out=6),
         labels=round(seq(maximumvalue,0,length.out=6)),las=2)
    axis(2,at=seq(0,1,length.out=6),
         labels=c("IBA","IBB","IBD","IBE","IBH","MIX"),las=2)
    dev.off()
    ##64-FF
  }
  if(do_mix==T){
    datasetOrder = c("GENOME",1,"1A","1B",2:4,"4A",5:28,"LG2","LG5","LGE22","mtDNA","Z",
                     "MORPH","PC1","PC2","PC3","LS0","LS1","LS2","LS3",
                     "HIGH","LOW","NORM","HIGH75","LOW75","NORM75","HIGH50","LOW50","NORM50")
    pdf(paste("~/heatmaps_total_ngsdist_separate_mixed_",date,".pdf",sep=""),width=13,height=10)
    par(mfrow=c(5,2))
    for(spp in sort(unique(bi$SPECIES))){
      if(spp != ""){
        print(spp)
        bi_spp = bi[bi$SPECIES==spp,]
        bi_spp = bi_spp[,c("DATASET",
                           "MOSTTOT.1MODEL","MOSTTOT",
                           "MOSTATOT.1MODEL","MOSTATOT")]
        bi_spp$DATASET=sub("CHR","",bi_spp$DATASET)
        ## find and add missing datasets
        missing_chrom = datasetOrder[which(!(datasetOrder %in% unique(bi_spp$DATASET)))]
        if(length(missing_chrom)>1){
          missing_chrom = cbind(DATASET=missing_chrom,
                                MOSTTOT.1MODEL=6,MOSTTOT="",
                                MOSTATOT.1MODEL=6,MOSTATOT="")
          bi_spp = rbind(bi_spp,missing_chrom)
        }
        bi_spp = bi_spp[ order(match(bi_spp$DATASET, datasetOrder)), ]
        bi_spp$ABUN = grepl("ABUN",bi_spp$MOSTTOT)
        bi_spp$STR = grepl("STR",bi_spp$MOSTTOT)
        bi_spp$LGM = grepl("LGM",bi_spp$MOSTTOT)
        bi_spp$PRES = grepl("PRES",bi_spp$MOSTTOT)
        bi_spp$ENV = grepl("ENV",bi_spp$MOSTTOT)
        bi_spp$IBD = grepl("IBD",bi_spp$MOSTTOT)
        bi_spp$ABUN[bi_spp$ABUN == TRUE] = 1
        bi_spp$STR[bi_spp$STR == TRUE] = 2
        bi_spp$IBD[bi_spp$IBD == TRUE] = 3
        bi_spp$ENV[bi_spp$ENV == TRUE] = 4
        bi_spp$PRES[bi_spp$PRES == TRUE] = 4
        bi_spp$LGM[bi_spp$LGM == TRUE] = 5
        bi_spp[bi_spp == FALSE] = 7
        bi_spp$MOSTATOT.1MODEL[bi_spp$MOSTATOT.1MODEL=="IBA"]=1
        bi_spp$MOSTATOT.1MODEL[bi_spp$MOSTATOT.1MODEL=="IBB"]=2
        bi_spp$MOSTATOT.1MODEL[bi_spp$MOSTATOT.1MODEL=="IBD"]=3
        bi_spp$MOSTATOT.1MODEL[bi_spp$MOSTATOT.1MODEL=="IBE"]=4
        bi_spp$MOSTATOT.1MODEL[bi_spp$MOSTATOT.1MODEL=="IBH"]=5
        bi_spp$MOSTATOT.1MODEL[bi_spp$MOSTATOT.1MODEL=="MIXED"]=6
        bi_spp$MOSTATOT.1MODEL = as.numeric(bi_spp$MOSTATOT.1MODEL)
        colorlist = c(purple,orange,blue,yellow,yellow,red,"black","white")
        par(mar=c(5,5,2,0))
        image(as.matrix(bi_spp[,c("ABUN","STR","IBD","ENV","PRES","LGM","MOSTATOT.1MODEL")]),
              col=colorlist,yaxt="n",zlim=c(1,7),xaxt="n",main=spp)
        axis(2,at=seq(0,1,length.out=7),
             labels=c("IBA","IBB","IBD","IBE-E","IBE-P","IBH","OVERALL"),las=2)
        axis(1,at=seq(0,1,length.out=54),
             labels=datasetOrder,las=2)
      }
    }
    dev.off()
  }
}
## calculate chi-squared tests
if(do_chi == T){
  bi=read.table("bivariate_gdm_results_NGSDIST_20July2022.csv",sep=",",header=T)
  bi_trimmed=bi[,c("DATASET","SPECIES","STRTYPE","MOSTA.1MODEL1","MOSTA.1MODEL2","MOSTA.1MODEL3")]
  bi_trimmed[bi_trimmed==""] = NA
  bi_trimmed$MOSTA.1MODEL1[grepl(" ",bi_trimmed$MOSTA.1MODEL1)] = "MIXED"
  bi_trimmed$MOSTA.1MODEL2[grepl(" ",bi_trimmed$MOSTA.1MODEL2)] = "MIXED"
  bi_trimmed$MOSTA.1MODEL3[grepl(" ",bi_trimmed$MOSTA.1MODEL3)] = "MIXED"
  species_tab_1 = table(bi_trimmed[,c("SPECIES","MOSTA.1MODEL1")])
  species_tab_2 = table(bi_trimmed[,c("SPECIES","MOSTA.1MODEL2")])
  species_tab_3 = table(bi_trimmed[,c("SPECIES","MOSTA.1MODEL3")])
  structure_tab_1 = table(bi_trimmed[,c("STRTYPE","MOSTA.1MODEL1")])
  structure_tab_2 = table(bi_trimmed[,c("STRTYPE","MOSTA.1MODEL2")])
  structure_tab_3 = table(bi_trimmed[,c("STRTYPE","MOSTA.1MODEL3")])
  dataset_tab_1 = table(bi_trimmed[,c("DATASET","MOSTA.1MODEL1")])
  dataset_tab_2 = table(bi_trimmed[,c("DATASET","MOSTA.1MODEL2")])
  dataset_tab_3 = table(bi_trimmed[,c("DATASET","MOSTA.1MODEL3")])
  spp1_sim = chisq.test(species_tab_1,simulate.p.value = T)
  spp1_raw = chisq.test(species_tab_1,simulate.p.value = F)
  spp2_sim = chisq.test(species_tab_2,simulate.p.value = T)
  spp2_raw = chisq.test(species_tab_2,simulate.p.value = F)
  spp3_sim = chisq.test(species_tab_3,simulate.p.value = T)
  spp3_raw = chisq.test(species_tab_3,simulate.p.value = F)
  str1_sim = chisq.test(structure_tab_1,simulate.p.value = T)
  str1_raw = chisq.test(structure_tab_1,simulate.p.value = F)
  str2_sim = chisq.test(structure_tab_2,simulate.p.value = T)
  str2_raw = chisq.test(structure_tab_2,simulate.p.value = F)
  str3_sim = chisq.test(structure_tab_3,simulate.p.value = T)
  str3_raw = chisq.test(structure_tab_3,simulate.p.value = F)
  dat1_sim = chisq.test(dataset_tab_1,simulate.p.value = T)
  dat1_raw = chisq.test(dataset_tab_1,simulate.p.value = F)
  dat2_sim = chisq.test(dataset_tab_2,simulate.p.value = T)
  dat2_raw = chisq.test(dataset_tab_2,simulate.p.value = F)
  dat3_sim = chisq.test(dataset_tab_3,simulate.p.value = T)
  dat3_raw = chisq.test(dataset_tab_3,simulate.p.value = F)
}
