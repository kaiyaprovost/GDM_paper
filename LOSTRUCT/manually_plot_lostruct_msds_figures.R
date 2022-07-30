path = "/Users/kprovost/Documents/MDSCOORDS/BELLII/TGUT/Vireo-bellii_PseudoNC"
setwd(path)
files = list.files(path=path,pattern="subset",full.names = T,recursive=T)
#files = files[(grepl("ALL",files))]
corner.cols <- RColorBrewer::brewer.pal(3,"Dark2")

for(file in files[1:length(files)]){
  print(file)
  try({
  df = read.table(file,header=T,sep=",",fill=T)
  df = df[!(df$chrom %in% c("black","empty",corner.cols)),]
  df$ccols = NA
  df$ccols4 = NA
  df$cpch = NA
  df$cpch4 = NA
   
  for(spp in sort(unique(df$species))){
    dfspp = df[df$species==spp,]
    #spp = strsplit(basename(file),"\\.")[[1]][2]
    
    mds.corners <- lostruct::corners( dfspp[,c("MDS1","MDS2")], prop=.05 ,k=3)
    corner.cols <- RColorBrewer::brewer.pal(3,"Dark2")
    corner.pch <- c(15,17,19)
    ccols <- rep("black",nrow(dfspp))
    cpch <- rep(20,nrow(dfspp))
    for (k in 1:ncol(mds.corners)) {
      ccols[ mds.corners[,k] ] <- corner.cols[k]
      cpch[ mds.corners[,k] ] <- corner.pch[k]
    }

    mds.corners4 <- lostruct::corners( dfspp[,c("MDS1","MDS2")], prop=.05 ,k = 4)
    corner.cols4 <- RColorBrewer::brewer.pal(4,"Dark2")
    corner.pch4 <- c(15,17,19,21)
    ccols4 <- rep("black",nrow(dfspp))
    cpch4 <- rep(20,nrow(dfspp))
    for (k in 1:ncol(mds.corners4)) {
      ccols4[ mds.corners4[,k] ] <- corner.cols4[k]
      cpch4[ mds.corners4[,k] ] <- corner.pch4[k]
    }
    
    png(paste(file,"_",spp,"_MSDS_plot_manual_lostruct_3.png",sep=""),height=3,width=3,units="in",res=600)
    par(mar=c(4,4,1,0))
    plot(dfspp$MDS1,dfspp$MDS2,
         col=ccols,main="",pch=cpch,
         ylab="MDS2",xlab="MDS1",las=2)
    dev.off()
    png(paste(file,"_",spp,"_MSDS_plot_manual_lostruct_4.png",sep=""),height=3,width=3,units="in",res=600)
    par(mar=c(4,4,1,0))
    plot(dfspp$MDS1,dfspp$MDS2,
         col=ccols4,main="",pch=cpch4,
         ylab="MDS2",xlab="MDS1",las=2)
    dev.off()
    
    png(paste(file,"_",spp,"_MSDS_boxplot_manual_lostruct.png",sep=""))
    par(mar=c(4,4,1,0),mfrow=(c(2,1)))
    boxplot(dfspp$MDS1~dfspp$chrom,main=spp,ylab="MDS1",xlab="Chrom",las=2)
    boxplot(dfspp$MDS2~dfspp$chrom,main=spp,ylab="MDS2",xlab="Chrom",las=2)
    dev.off()
    
    ## 3: bil, bru
    ## 4: cri
    ## not enough: sin, fus, nit, mel, cur
    ## weird: fla
    

    
    if(spp=="Vireo-bellii") {
      #dfsm = dfspp[dfspp$chrom!="mtDNA",]
      png(paste(file,"_",spp,"_MSDS_plot_manual_lostruct_3_NOMTDNA.png",sep=""),height=3,width=3,units="in",res=600)
      par(mar=c(4,4,1,0))
      plot(dfspp$MDS1[dfspp$chrom!="mtDNA"],dfspp$MDS2[dfspp$chrom!="mtDNA"],
           col=ccols[dfspp$chrom!="mtDNA"],main="",pch=cpch[dfspp$chrom!="mtDNA"],
           ylab="MDS2",xlab="MDS1",las=2)
      dev.off()
      
      png(paste(file,"_",spp,"_MSDS_plot_manual_lostruct_4_NOMTDNA.png",sep=""),height=3,width=3,units="in",res=600)
      par(mar=c(4,4,1,0))
      plot(dfspp$MDS1[dfspp$chrom!="mtDNA"],dfspp$MDS2[dfspp$chrom!="mtDNA"],
           col=ccols4[dfspp$chrom!="mtDNA"],main="",pch=cpch4[dfspp$chrom!="mtDNA"],
           ylab="MDS2",xlab="MDS1",las=2)
      dev.off()
      
      png(paste(file,"_",spp,"_MSDS_boxplot_manual_lostruct_NOMTDNA.png",sep=""))
      par(mar=c(4,4,1,0),mfrow=(c(2,1)))
      boxplot(dfsm$MDS1~dfsm$chrom,main=spp,ylab="MDS1",xlab="Chrom",las=2)
      boxplot(dfsm$MDS2~dfsm$chrom,main=spp,ylab="MDS2",xlab="Chrom",las=2)
      dev.off()
      
      
      
      
      
    }
    
    df$ccols[df$species==spp] = ccols
    df$ccols4[df$species==spp] = ccols4
    df$cpch[df$species==spp] = cpch
    df$cpch4[df$species==spp] = cpch4
    write.table(df,paste(file,"_COLORLOSTRUCT.csv",sep=""),sep=",",row.names = F,quote=F)
  }
  })
  
  
  
}


## final models 
filelist = list.files(path="/Users/kprovost/Documents/MDSCOORDS",
                      pattern="COLORLOSTRUCT.csv$",full.names = T)

## make the files
for(fi in filelist){
  df = data.table::fread(fi,header=T,sep=",",data.table=F)
  #png(paste(fi,"_lostruct_4.png"))
  #plot(df$MDS1,df$MDS2,col=df$ccols4,pch=df$cpch4-15)
  #dev.off()
  print(fi)
  print(table(df$ccols,df$ccols4))
  
  
  tab=table(df$chrom,df$ccols4)
  size = rowSums(tab)
  tab2 = tab/size
  png(paste(fi,"_lostruct_4_BARPLOT_sizescale.png"),width=1400)
  barplot(t(tab2),col=colnames(tab2),width=size+500)
  dev.off()
}

## each window is 50,000 wide and do not overlap
## so window 0 runs from starting 0 to starting 50000

for(fi in filelist[6:10]){
  df = data.table::fread(fi,header=T,sep=",",data.table=F)
  df_fix = data.frame()
  for(k in 0:4) {
    df_k = df
    df_k$windowstart = (as.numeric(df_k$window)*50000)+(k*10000)
    df_k$windowend = df_k$windowstart + 100000
    df_k$midPos = df_k$windowstart + 50000
    df_fix = gtools::smartbind(df_fix,df_k)
  }
  df_fix = unique(df_fix)
  df_fix = df_fix[order(df_fix$chrom,df_fix$windowstart),]
  write.table(df_fix,paste(fi,"_FixWindows.csv",sep=""),sep=",",row.names = F,quote = F)
}

filelist2 = list.files(path="/Users/kprovost/Documents/MDSCOORDS",
                      pattern="_FixWindows.csv$",full.names = T)
## merge into one big window thing
df_list = lapply(filelist2,FUN=function(x){data.table::fread(x,sep=",",header=T,data.table=F)})
df_new = do.call(gtools::smartbind,df_list)
df_new = df_new[,c("species","chrom","ccols","ccols4","windowstart","windowend","midPos")]
colnames(df_new) = c("species","chr","color","color4","windowstarts","windowstops","midPos")
write.table(df_new,"/Users/kprovost/Documents/MDSCOORDS/colors.txt",sep=",",row.names = F,quote=F)
