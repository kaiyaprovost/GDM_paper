## lostruct on genotype likelihoods from PCAngsd
library(R.utils)
library(tools)
library(lostruct)
library(data.table)
path="~/COV/" ## Change this as needed
setwd(path)
spplist = c("bellii","bilineata","brunneicapillus","crissale","curvirostre",
            "flaviceps","fusca","melanura","nitens","sinuatus"
)
dirs = list.dirs(path,full.names = T,recursive = T)
for(dir in dirs) {
  print(dir)
  list_cov_files = list.files(path=dir,pattern="PCAngsd.1.cov$",
                              recursive=F,full.names = T)
  list_cov_files= list_cov_files[file.exists(list_cov_files) && !dir.exists(list_cov_files)]
  for(species in spplist){
    print(species)
    numPCS=2
    spp_cov_files = list_cov_files[grepl(species,list_cov_files)]
    if(length(spp_cov_files)>0){
      ## generate chromosomes etc
      try({
        print("GENERATING SPLIT DATA")
        splitdata = lapply(spp_cov_files,FUN=function(x){
          y=basename(x)
          if(grepl("NOWEIRD",y,fixed=T)){
            thisspp=species
          } else {
            thisspp = paste(species,"-NOWEIRD",sep="")
          }
          splits=strsplit(y,"1_Tgut_")[[1]]
          if(length(splits)>1){
            chromsplit = strsplit(splits[2],"\\.")[[1]]
            chrom=chromsplit[1]
            if(chromsplit[3]=="F"){window="F"} else {window=chromsplit[2]}
          } else {
            chrom="GENOME"
            window="F"
          }
          z=cbind(file=y,
                  species=thisspp,
                  chrom=chrom,
                  window=window)
        })
        splitdata_df = as.data.frame(do.call(rbind,splitdata))
      })
      try({
        print("CALCULATING EIGENVECTORS")
        cov_df_list = lapply(spp_cov_files,FUN=function(x){
          as.matrix(data.table::fread(x,header=F,data.table=F)) ## supports gz and normal files
        })
        window_eigs <- lapply(cov_df_list, function (cm) lostruct::cov_pca(k=numPCS, covmat=cm))
        window_eigs_df = do.call(rbind,window_eigs)
        eig.file <- file.path(dir,paste(species,"eigenvectors.txt",sep=""))
        write.table(window_eigs_df,eig.file)
      })
      try({
        for(spp in unique(splitdata_df$species)){
          print(spp)
          mds.file <- file.path(dir,paste(species,"_subset",spp,"_mds_coords.csv",sep=""))
          this_spp_nums = which(splitdata_df$species==spp)
          splitdata_df_spp = splitdata_df[this_spp_nums,]
          window_eigs_df_spp = window_eigs_df[this_spp_nums,]
          print("DOING PCA DISTANCE MATRIX")
          pc.distmat <- lostruct::pc_dist(x=window_eigs_df_spp, npc=numPCS, w = 1, normalize = "L1", mc.cores = 1)
          na.inds <- is.na( window_eigs_df_spp[,1] ) # there may be windows with missing data
          print("CALCULATING MDS PER SPECIES")
          mds.coords <- cbind( data.frame(splitdata_df_spp),
                               cmdscale( pc.distmat[!na.inds,!na.inds], k=numPCS )[ ifelse( na.inds, NA, cumsum(!na.inds) ), ]
          )
          colnames(mds.coords)[-(1:4)] <- paste0("MDS",seq_len(numPCS))
          write.csv( mds.coords, mds.file, row.names=FALSE )
        }
      })
    }
  }
}
