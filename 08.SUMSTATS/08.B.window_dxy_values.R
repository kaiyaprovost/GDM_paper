specieslist = sort(c("bel","bil","bru",
                     "cri","cur","fla","fus",
                     "mel","nit","sin"))
chromlist=sort(c("1","1A","1B","2","3","4","4A","5",
                 "6","7","8","9","10","11","12","13",
                 "14","15","16","17","18","19","20",
                 "21","22","23","24","25","26","27",
                 "28","LG2","LG5","LGE22","mtDNA","Z"))
orig_num_to_keep = 40
windowsize = 100000
movesize = 10000
printstatus = 10
listfiles = list.files(pattern="Dxy_persite.txt",full.names = T,recursive = T)
num_to_keep = orig_num_to_keep
for(chrom in chromlist) { 
  chromfiles = listfiles[grep(chrom,listfiles)]
  if(length(chromfiles)>0){
    for(spp in specieslist) {
      sppfiles = chromfiles[grep(spp,chromfiles)]
      if(length(sppfiles)>0){
        for(file in sppfiles) {
          print(file)
          if(file.exists(file)){
            print("READING CSV")
            df=data.table::fread(file,sep="\t",header=T,stringsAsFactors=F,
                                 colClasses=c("character","numeric","numeric"),
                                 showProgress=T)
            chrom=strsplit(basename(file),"_")[[1]][2]
            dataset=strsplit(basename(file),"_")[[1]][4]
            df$position = as.numeric(as.character(df$position))
            df$dxy = as.numeric(as.character(df$dxy))
            df = df[complete.cases(df),]
            if(nrow(df) <=0){
              print("NO DATA")
              next
            } else {
              ## how many rows are just zero? 
              justzero=sum(df$dxy==0)
              percentzero=justzero/nrow(df)
              print(paste("Proportion that is zero:",round(percentzero,2)))
              only_zero_data = df[df$dxy==0,]
              only_full_data = df[df$dxy!=0,]
              print("READING IN SCAFFOLDS")
              num_on_scaff = rev(sort(table(df$chromo)))
              totalscaffs = length(num_on_scaff)
              print(paste("TOTAL SCAFFS READ:",totalscaffs))
              if (totalscaffs != num_to_keep) {
                num_to_keep=totalscaffs
              }
              outfile = paste(spp,"_",dataset,"_Dxy_WINDOWS_",num_to_keep,"-",chrom,".txt",sep="")
              outfile_TEMP = paste(spp,"_",dataset,"_Dxy_WINDOWS_",num_to_keep,"-",chrom,".temp",sep="")
              pngfile = paste(spp,"_",dataset,"_DxyMEAN_WINDOWS_",num_to_keep,"-",chrom,".png",sep="")
              sumfile = paste(spp,"_",dataset,"_DxySUMS_WINDOWS_",num_to_keep,"-",chrom,".png",sep="")
              sdvfile = paste(spp,"_",dataset,"_DxySDVS_WINDOWS_",num_to_keep,"-",chrom,".png",sep="")
              rawfile = paste(spp,"_",dataset,"_DxyRAW_",num_to_keep,"-",chrom,".png",sep="")
              if(file.exists(pngfile) && overwrite==F) {
                print("SKIPPING")
              } else {
                if(num_to_keep != 1){
                  print("SUBSETTING SCAFFOLDS")
                  scaff_subset = names(num_on_scaff)[1:num_to_keep]
                  df_z = only_zero_data[only_zero_data$chromo %in% scaff_subset,]
                  df_f = only_full_data[only_full_data$chromo %in% scaff_subset,]
                  matches = match(df_f$chromo,scaff_subset)
                  df_f$chromnum = matches
                  df_f = df_f[order(df_f[,"chromnum"], df_f[,"position"] ),]
                  head(df_f)
                  summary(df_f$dxy)
                  summary(df_f$position)
                  nscaf = (length(scaff_subset))
                } else {
                  nscaf=1
                  df_z = only_zero_data
                  df_f = only_full_data
                  scaff_subset = df_f$chromo[1]
                }
                closest_snp = min(min(df_f$position),min(df_z$position))
                furthest_snp = max(max(df_f$position),max(df_z$position))
                print(paste("NUMBER OF SCAFFS TO OUTPUT:",nscaf))
                nums = as.numeric(num_on_scaff)
                scaf_list = c()
                start_list = c()
                mean_list = c()
                dev_list = c()
                snp_list = c()
                sum_list = c()
                print("writing rawfile")
                png(rawfile,width=700,height=350)
                palette(c("red","orange","goldenrod","green","blue","purple",
                          "cyan","black","brown","magenta"))
                plot(x=as.numeric(df_f$position),y=as.numeric(df_f$dxy),
                     col=as.numeric(as.factor(df_f$chromo)),cex=0.2,
                     main=spp,xlab="Position",ylab="Raw DXY NO ZERO",
                     ylim=c(0,0.5))
                text(x=mean(as.numeric(df_f$position),na.rm=T),
                     y=0.9,labels=as.character(round(as.numeric(mean(df_f$dxy,na.rm=T)),2)))
                dev.off()
                if(rawonly==F){
                  print("BEGINNING SCAFFOLDS")
                  for (i in 1:length(scaff_subset)) {
                    print(paste(i," of ",nscaf,spp))
                    scaf = (scaff_subset[i])
                    number = nums[i]
                    if(nscaf==1){
                      this_scaf_f = df_f
                      this_scaf_z = df_z
                    } else {
                      this_scaf_f = df_f[df_f$chromo==scaf,]
                      this_scaf_z = df_z[df_z$chromo==scaf,]
                    }
                    furthest = max(max(this_scaf_f$position),max(this_scaf_z$position))
                    ## get windows
                    window_ends = seq(windowsize,max(furthest,windowsize),movesize)
                    window_starts = seq(1,furthest,movesize)
                    window_starts = window_starts[1:length(window_ends)]
                    windows = as.data.frame(cbind("starts"=window_starts,"ends"=window_ends))
                    totwindows=(nrow(windows))
                    print(totwindows)
                    for (j in 1:nrow(windows)){
                      window = windows[j,]
                      if(j %% printstatus == 0) {
                        print(paste(j,"/",totwindows))
                      }
                      start = window$starts
                      end = window$ends
                      nonzero = this_scaf_f[this_scaf_f$position>=start & this_scaf_f$position<=end,]
                      numzero = sum(this_scaf_z$position>=start & this_scaf_z$position<=end)
                      num_snps = numzero + nrow(nonzero)
                      sum_dxy = sum(nonzero$dxy,na.rm=T)
                      average_dxy = sum_dxy/num_snps
                      deviation_dxy = sqrt((sum((c(nonzero$dxy,rep(0,numzero))-average_dxy)^2))/num_snps)
                      if(j==1 && i==1){
                        write.table(cbind("scafs","starts","means","stdvs","snps","sums"),outfile_TEMP,
                                    row.names = F,col.names = F,append=T)
                      }
                      write.table(cbind(scaf,start,average_dxy,deviation_dxy,num_snps,sum_dxy),
                                  outfile_TEMP,
                                  append=T,row.names = F,col.names = F)
                    }
                  }
                  output_dataframe = read.table(outfile_TEMP,header = T,stringsAsFactors = F)
                  output_dataframe=unique(output_dataframe)
                  if(plotMeans==T){
                    png(pngfile,width=700,height=700)
                    par(mfrow=c(2,1))
                    palette(c("red","orange","goldenrod","green","blue","purple",
                              "cyan","black","brown","magenta"))
                    plot(as.numeric(output_dataframe$means),
                         col=as.numeric(as.factor(output_dataframe$scafs)),cex=0.2,
                         main=spp,xlab="Window (Scaffold)",ylab="Mean DXY",
                         ylim=c(0,0.5))
                    plot(as.numeric(output_dataframe$means),
                         col=as.numeric(as.factor(output_dataframe$scafs)),cex=0.2,
                         main=spp,xlab="Window (Scaffold)",ylab="Mean DXY")
                    dev.off()
                    png(sumfile,width=700,height=350)
                    palette(c("red","orange","goldenrod","green","blue","purple",
                              "cyan","black","brown","magenta"))
                    plot(as.numeric(output_dataframe$sums),
                         col=as.numeric(as.factor(output_dataframe$scafs)),cex=0.2,
                         main=spp,xlab="Window (Scaffold)",ylab="Sum DXY")
                    dev.off()
                    png(sdvfile,width=700,height=350)
                    palette(c("red","orange","goldenrod","green","blue","purple",
                              "cyan","black","brown","magenta"))
                    plot(as.numeric(output_dataframe$stdvs),
                         col=as.numeric(as.factor(output_dataframe$scafs)),cex=0.2,
                         main=spp,xlab="Window (Scaffold)",ylab="STDV DXY")
                    dev.off()
                  }
                  write.table(output_dataframe,outfile,append=T)
                }
              }
            }
            file.rename(file,paste("DONE/",basename(file),sep=""))
          } else {print("file does not exist, skipping")}
        }
      }
    }
  }
}