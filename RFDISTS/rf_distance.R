rm(list=ls())
Sys.setenv('R_MAX_VSIZE'=64000000000)
spp="bellii"

## convert distance matrices to newicks

dists = list.files(pattern="ngsdist",recursive = T,full.names = T)
dists = dists[grepl(spp,dists)]

for(dist in dists) {
  df = read.table(dist,sep="\t",header=T,row.names = 1)
  df = as.dist(df)
  tree=ape::nj(df)
  write.tree(tree,paste(dist,".newick",sep=""))
}

newickfiles1 = list.files(pattern="newick$",recursive = T,full.names = T)
newickfiles1 = newickfiles1[grepl(spp,newickfiles1)]

trees = lapply(newickfiles,FUN=read.NEWICK)
class(trees) = "multiPhylo"
phydist=phytools::multiRF(trees)
colnames(phydist) = names
rownames(phydist) = names
phydist=as.matrix(phydist)

if(nrow(phydist)!=length(allnewicks)){
rowtoadd=rep(NA,ncol(phydist))
for(missing in missingcolnames){
  phydist=rbind(phydist,rowtoadd)
}
coltoadd=rep(NA,nrow(phydist))
for(missing in missingcolnames){
  phydist=cbind(phydist,coltoadd)
}

names=c(names,missingcolnames)
colnames(phydist) = names
rownames(phydist) = names

}

phydist = phydist[,gtools::mixedorder(colnames(phydist))]
phydist = phydist[gtools::mixedorder(rownames(phydist)),]

allcol=which(colnames(phydist)=="all")
neworder=c(1:(allcol-1),(allcol+1):length(colnames(phydist)),allcol)
phydist=phydist[neworder,neworder]
write.table(phydist,paste("phydist_rfdist_",spp,".txt",sep=""),sep="\t",row.names = T)



