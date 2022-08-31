spplist=(c("BELLII"))
chromlist=(c(1,"1A","1B",2:4,"4A",5:28,
"LG2","LG5","LGE22","mtDNA","Z"
))
dynamic_require <- function(package,lib=NULL) {
if (eval(parse(text = paste("require(", package, ")"))))
return(TRUE)
if(is.null(lib)) {
install.packages(package,repos='http://cran.us.r-project.org')
} else {
install.packages(package,lib=lib,repos='http://cran.us.r-project.org')
}
return(eval(parse(text = paste(
"require(", package,  ")"))))
}
packages = c("vegan","sp","rgdal","raster","parallel","iterators","gdm","foreach","ecodist","doParallel")
for (p in packages) {
dynamic_require(p)
}
print(version)
outersect <- function(x, y) {unique(sort(c(setdiff(x, y),setdiff(y, x))))}
gdm.varImp.edit = function (spTable, geo, splines = NULL, knots = NULL, fullModelOnly = FALSE, 
nPerm = 50, parallel = FALSE, cores = 2, sampleSites = 1, 
sampleSitePairs = 1, outFile = NULL) {
k <- NULL
if (class(spTable)[1] != "gdmData") {
warning("spTable class does not include type 'gdmData'. Make sure your data is in site-pair format or the gdm model will not fit.")
}
if (!(class(spTable)[1] == "gdmData" | class(spTable)[1] == 
"matrix" | class(spTable)[1] == "data.frame")) {
stop("spTable argument needs to be gdmData, a matrix, or a data frame")
}
if (ncol(spTable) < 6) {
stop("Not enough columns in data. (Minimum need: Observed, weights, X0, Y0, X1, Y1)")
}
if (nrow(spTable) < 1) {
stop("Not enough rows in data")
}
if (!(geo == TRUE | geo == FALSE)) {
stop("geo argument must be either TRUE or FALSE")
}
if (is.null(splines) == FALSE & class(splines) != "numeric") {
stop("argument splines needs to be a numeric data type")
}
if (is.null(knots) == FALSE & class(knots) != "numeric") {
stop("argument knots needs to be a numeric data type")
}
if (!(fullModelOnly == TRUE | fullModelOnly == FALSE)) {
stop("fullModelOnly argument must be either TRUE or FALSE")
}
if ((is.null(nPerm) == FALSE & is.numeric(nPerm) == FALSE) | 
nPerm < 1) {
stop("argument nPerm needs to be a positive integer")
}
if (!(parallel == TRUE | parallel == FALSE)) {
stop("parallel argument must be either TRUE or FALSE")
}
if (parallel == TRUE & is.null(cores) == TRUE) {
stop("If parallel==TRUE, the number of cores must be specified")
}
if ((is.null(cores) == FALSE & is.numeric(cores) == FALSE) | 
cores < 1) {
stop("argument cores needs to be a positive integer")
}
if (is.numeric(sampleSites) == FALSE | sampleSites < 0 | 
sampleSites > 1) {
stop("argument sampleSites needs to be a positive number between 0 and 1")
}
if (is.numeric(sampleSitePairs) == FALSE | sampleSitePairs < 
0 | sampleSitePairs > 1) {
stop("argument sampleSitePairs needs to be a positive number between 0 and 1")
}
if (sampleSites == 0) {
stop("a sampleSites value of 0 will remove all sites from the analysis")
}
if (sampleSitePairs == 0) {
stop("a sampleSitePairs value of 0 will remove all sites from the analysis")
}
if (is.null(outFile) == FALSE) {
if (is.character(outFile) == FALSE) {
stop("argument outFile needs to be a character string of the directory and file name you wish the tables to be written to")
}
outFileChar <- nchar(outFile)
if (substr(outFile, outFileChar - 5, outFileChar) != 
".RData") {
outFile <- paste(outFile, ".RData", sep = "")
}
if (length(strsplit(outFile, "/")[[1]]) > 1) {
splitOutFile <- strsplit(outFile, "/")[[1]][-length(strsplit(outFile, 
				   "/")[[1]])]
dir.create(paste(splitOutFile, collapse = "/"))
}
else {
outFile <- paste("./", outFile, sep = "")
}
}
nPerm <- as.integer(nPerm)
cores <- as.integer(cores)
if (sampleSites < 1) {
spTable <- removeSitesFromSitePair(spTable, sampleSites = sampleSites)
if (sampleSitePairs < 1) {
warning("You have selected to randomly remove sites and/or site-pairs.")
}
}
if (sampleSitePairs < 1) {
numRm <- sample(1:nrow(spTable), round(nrow(spTable) * 
(1 - sampleSitePairs)))
spTable <- spTable[-c(numRm), ]
}
rtmp <- spTable[, 1]
if (length(rtmp[rtmp < 0]) > 0) {
stop("Response spTable has negative values. Must be between 0 - 1.")
}
if (length(rtmp[rtmp > 1]) > 0) {
stop("Response spTable has values greater than 1. Must be between 0 - 1.")
}
nVars <- (ncol(spTable) - 6)/2
varNames <- colnames(spTable[c(7:(6 + nVars))])
varNames <- sapply(strsplit(varNames, "s1."), "[[", 2)
if (geo == TRUE) {
nVars <- nVars + 1
varNames <- c("Geographic", varNames)
}
sortMatX <- sapply(1:nrow(spTable), function(i, spTab) {
c(spTab[i, 3], spTab[i, 5])
}, spTab = spTable)
sortMatY <- sapply(1:nrow(spTable), function(i, spTab) {
c(spTab[i, 4], spTab[i, 6])
}, spTab = spTable)
sortMatNum <- sapply(1:nrow(spTable), function(i) {
c(1, 2)
})
sortMatRow <- sapply(1:nrow(spTable), function(i) {
c(i, i)
})
fullSortMat <- cbind(as.vector(sortMatX), as.vector(sortMatY), 
as.vector(sortMatNum), as.vector(sortMatRow), rep(NA, 
						 length(sortMatX)))
siteByCoords <- as.data.frame(unique(fullSortMat[, 1:2]))
numSites <- nrow(siteByCoords)
for (i in 1:numSites) {
fullSortMat[which(fullSortMat[, 1] == siteByCoords[i, 
	   1] & fullSortMat[, 2] == siteByCoords[i, 2]), 5] <- i
}
indexTab <- matrix(NA, nrow(spTable), 2)
for (iRow in 1:nrow(fullSortMat)) {
indexTab[fullSortMat[iRow, 4], fullSortMat[iRow, 3]] <- fullSortMat[iRow, 
						5]
}
rm(fullSortMat)
rm(sortMatX)
rm(sortMatY)
rm(sortMatNum)
rm(sortMatRow)
rm(siteByCoords)
exBySite <- lapply(1:numSites, function(i, index, tab) {
rowSites <- which(index[, 1] %in% i)
if (length(rowSites) < 1) {
rowSites <- which(index[, 2] %in% i)
}
exSiteData <- tab[rowSites[1], ]
return(exSiteData)
}, index = indexTab, tab = spTable)
outSite <- which(!(1:numSites %in% indexTab[, 1]))
for (i in 1:length(exBySite)) {
siteRow <- exBySite[[i]]
if (i %in% outSite) {
siteRow <- siteRow[grep("s2.", colnames(siteRow))]
colnames(siteRow) <- sapply(strsplit(colnames(siteRow), 
"s2."), "[[", 2)
}
else {
siteRow <- siteRow[grep("s1.", colnames(siteRow))]
colnames(siteRow) <- sapply(strsplit(colnames(siteRow), 
"s1."), "[[", 2)
}
exBySite[[i]] <- siteRow
}
siteData <- do.call("rbind", exBySite)
modelTestValues <- matrix(NA, 4, nVars, dimnames = list(c("Model deviance", 
			"Percent deviance explained", "Model p-value", "Fitted permutations"), 
		  c("fullModel", paste("fullModel-", seq(1, nVars - 1), 
							   sep = ""))))
devReductVars <- matrix(NA, nVars, nVars - 1)
rownames(devReductVars) <- varNames
colnames(devReductVars) <- c("fullModel", paste("fullModel-", 
  seq(1, nVars - 2), sep = ""))
pValues <- numPermsFit <- devReductVars
currSitePair <- spTable
nullGDMFullFit <- 0
for (v in 1:nVars) {
fullGDM <- gdm::gdm(currSitePair, geo = geo, splines = splines, 
knots = knots)
if (is.null(fullGDM) == TRUE) {
warning(paste("The model did not fit when testing variable number ", 
v, ". Terminating analysis. Returning output objects filled up to the point of termination.", 
sep = ""))
break
}
if (parallel == TRUE) {
cl <- makeCluster(cores, outfile = "")
registerDoParallel(cl)
permSitePairs <- foreach(k = 1:nPerm, .verbose = F, 
.packages = c("gdm"), .export = c("permutateSitePair")) %dopar% 
permutateSitePair(currSitePair, siteData, indexTab, 
varNames)
permGDM <- try(foreach(k = 1:length(permSitePairs), 
.verbose = F, .packages = c("gdm")) %dopar% gdm::gdm(permSitePairs[[k]], 
								  geo = geo, splines = NULL, knots = NULL))
stopCluster(cl)
}
else {
permSitePairs <- lapply(1:nPerm, function(i, csp, 
sd, it, vn) {
permutateSitePair(csp, sd, it, vn)
}, csp = currSitePair, sd = siteData, it = indexTab, 
vn = varNames)
permGDM <- lapply(permSitePairs, gdm, geo = geo, 
splines = NULL, knots = NULL)
}
permModelDev <- sapply(permGDM, function(mod) {
mod$gdmdeviance
})
modPerms <- length(which(sapply(permModelDev, is.null) == 
TRUE))
if (modPerms > 0) {
permModelDev <- unlist(permModelDev[-(which(sapply(permModelDev, 
		 is.null) == T))])
}
modelTestValues[1, v] <- fullGDM$gdmdeviance
modelTestValues[2, v] <- fullGDM$explained
modelTestValues[3, v] <- sum(permModelDev <= fullGDM$gdmdeviance)/(nPerm - 
						 modPerms)
modelTestValues[4, v] <- nPerm - modPerms
if (length(varNames) < 2) {
break
}
if (geo == TRUE) {
noGeoGDM <- gdm::gdm(currSitePair, geo = FALSE, splines = NULL, 
knots = NULL)
if (parallel == TRUE) {
cl <- makeCluster(cores, outfile = "")
registerDoParallel(cl)
permSitePairs <- foreach(k = 1:nPerm, .verbose = F, 
.packages = c("gdm"), .export = c("permutateSitePair")) %dopar% 
permutateSitePair(currSitePair, siteData, indexTab, 
varNames)
permGDM <- try(foreach(k = 1:length(permSitePairs), 
.verbose = F, .packages = c("gdm")) %dopar% 
gdm::gdm(permSitePairs[[k]], geo = geo, splines = NULL, 
knots = NULL))
stopCluster(cl)
}
else {
permSitePairs <- lapply(1:nPerm, function(i, 
  csp, sd, it, vn) {
permutateSitePair(csp, sd, it, vn)
}, csp = currSitePair, sd = siteData, it = indexTab, 
vn = varNames)
permGDM <- lapply(permSitePairs, gdm, geo = geo, 
splines = NULL, knots = NULL)
}
permModelDev <- sapply(permGDM, function(mod) {
mod$gdmdeviance
})
modPerms <- length(which(sapply(permModelDev, is.null) == 
TRUE))
if (modPerms > 0) {
permModelDev <- unlist(permModelDev[-(which(sapply(permModelDev, 
		   is.null) == T))])
}
if (is.null(noGeoGDM$gdmdeviance) == TRUE) {
permDevReduct <- -9999
devReductVars[1, v] <- -9999
pValues[1, v] <- -9999
}
else {
permDevReduct <- noGeoGDM$gdmdeviance - permModelDev
devReductVars[1, v] <- 100 * abs((noGeoGDM$explained - 
fullGDM$explained)/fullGDM$explained)
pValues[1, v] <- sum(permDevReduct >= (noGeoGDM$gdmdeviance - 
 fullGDM$gdmdeviance))/(nPerm - modPerms)
}
numPermsFit[1, v] <- nPerm - modPerms
}
for (varChar in varNames) {
if (varChar != "Geographic") {
testVarCols1 <- grep(paste("^s1.", varChar, "$", 
sep = ""), colnames(currSitePair))
testVarCols2 <- grep(paste("^s2.", varChar, "$", 
sep = ""), colnames(currSitePair))
testSitePair <- currSitePair[, -c(testVarCols1, 
testVarCols2)]
noVarGDM <- gdm::gdm(testSitePair, geo = geo, splines = NULL, 
knots = NULL)
if (parallel == TRUE) {
cl <- makeCluster(cores, outfile = "")
registerDoParallel(cl)
noVarSitePairs <- foreach(k = 1:nPerm, .verbose = F, 
.packages = c("gdm"), .export = c("permutateVarSitePair")) %dopar% 
permutateVarSitePair(currSitePair, siteData, 
indexTab, varChar)
permGDM <- try(foreach(k = 1:length(noVarSitePairs), 
.verbose = F, .packages = c("gdm")) %dopar% 
gdm::gdm(noVarSitePairs[[k]], geo = geo, splines = NULL, 
knots = NULL))
stopCluster(cl)
}
else {
noVarSitePairs <- lapply(1:nPerm, function(i, 
	 csp, sd, it, vn) {
permutateVarSitePair(csp, sd, it, vn)
}, csp = currSitePair, sd = siteData, it = indexTab, 
vn = varChar)
permGDM <- lapply(noVarSitePairs, gdm, geo = geo, 
splines = NULL, knots = NULL)
}
permModelDev <- sapply(permGDM, function(mod) {
mod$gdmdeviance
})
modPerms <- length(which(sapply(permModelDev, 
is.null) == TRUE))
if (modPerms > 0) {
permModelDev <- unlist(permModelDev[-(which(sapply(permModelDev, 
			 is.null) == T))])
}
if (is.null(noVarGDM$gdmdeviance) == TRUE) {
permDevReduct <- -9999
ggg <- which(rownames(devReductVars) %in% varChar)
devReductVars[ggg, v] <- rep(-9999, times = length(ggg))
pValues[ggg, v] <- rep(-9999, times = length(ggg))
}
else {
permDevReduct <- noVarGDM$gdmdeviance - permModelDev
devReductVars[which(rownames(devReductVars) %in% 
varChar), v] <- 100 * abs((noVarGDM$explained - 
			 fullGDM$explained)/fullGDM$explained)
pValues[which(rownames(pValues) %in% varChar), 
v] <- sum(permDevReduct >= (noVarGDM$gdmdeviance - 
fullGDM$gdmdeviance))/(nPerm - modPerms)
}
numPermsFit[which(rownames(numPermsFit) %in% 
varChar), v] <- nPerm - modPerms
}
}
if (fullModelOnly == TRUE) {
break
}
tempPVals <- as.numeric(pValues[c(1:nVars), v])
tempDevs <- as.numeric(devReductVars[c(1:nVars), v])
tempPVals <- tempPVals[!is.na(tempPVals)]
tempDevs <- tempDevs[!is.na(tempDevs)]
varToOmit <- which.max(tempPVals)
for (iCheck in 1:length(varNames)) {
if (tempPVals[iCheck] == tempPVals[varToOmit]) {
if (tempDevs[iCheck] < tempDevs[varToOmit]) {
varToOmit <- iCheck
}
}
}
if (varToOmit == 1 & geo == TRUE) {
geo <- FALSE
varNames <- varNames[-1]
}
else {
nameToRemove <- varNames[varToOmit]
varNames <- varNames[-varToOmit]
removeFromSitePs1 <- grep(paste("^s1.", nameToRemove, 
"$", sep = ""), colnames(currSitePair))
removeFromSitePs2 <- grep(paste("^s2.", nameToRemove, 
"$", sep = ""), colnames(currSitePair))
currSitePair <- currSitePair[, -c(removeFromSitePs1, 
removeFromSitePs2)]
}
## manually edit the splines and the knots if they are set so that the variables will run 
if(!(is.null(splines))){
splines = splines[-varToOmit]
}
if(!(is.null(knots))){
knots = knots[-varToOmit]
}
}
outObject <- list(modelTestValues, devReductVars, pValues, 
numPermsFit)
if (is.null(outFile) == FALSE) {
save(outObject, file = outFile)
}
return(outObject)
}
loadMorphData = function(spp,doGENE=F) {
if (doGENE == F) {
morph = "KLP_Master_Spreadsheet_Morphology_AGGREGATED_19_July_2019.csv"
morphdf = read.csv(morph)
morphdf = morphdf[,2:40]
} else {
morph = "CATALOG_BY_SPECIES.txt"
morphdf = read.csv(morph,sep = "\t")
colnames(morphdf)[1] = ("CATALOG.NUMBER")
}
morphdf = unique(morphdf)
morphdf = morphdf[morphdf$CATALOG.NUMBER != "",]
morphlocs = morphdf[,c("LAT","LONG")]
morphspp = morphdf$SPP
return(list(morphdf,morphlocs,morphspp))
}
loadEnvData = function(morphlocs,morphspp,morphdf,doGENE=F,spp,verbose=F,doEnvpcs=T,numenvpcs) {
## pure environments 
envfile = "ENMS_multilayer_onlyworldclim.tif"
env = raster::stack(envfile)
morphlocs[,1] = as.numeric(as.character(morphlocs[,1]))
morphlocs[,2] = as.numeric(as.character(morphlocs[,2]))
envdata = as.data.frame(raster::extract(env,morphlocs[,2:1]))
envdata$SPP = morphspp
envdata$CATALOG.NUMBER = morphdf$CATALOG.NUMBER
inputdata = merge(morphdf,envdata)
if(verbose==T){print("fin input data")}
data = unique(inputdata)
data = data[data$CATALOG.NUMBER != "",]
rownames(data) = make.names(data$CATALOG.NUMBER,unique=T)
if(verbose==T){print("fin rownames")}
if (doGENE == F) {
bioclim = data[,c(40:58)] 
} else {
bioclim = data[,c(5:23)] 
}
colnames(bioclim) = c("BIO1","BIO2","BIO3","BIO4","BIO5","BIO6","BIO7","BIO8","BIO9","BIO10","BIO11","BIO12","BIO13","BIO14","BIO15","BIO16","BIO17","BIO18","BIO19")
bioclim = bioclim[complete.cases(bioclim),]
outenvfile="raw_bioclim_values.txt"
write.csv(x=bioclim,file=outenvfile)
if(verbose==T){print("checking pcas")}
if (doGENE == F) {
m = cor(bioclim,use = "pairwise.complete.obs")
xy <- t(combn(colnames(m),2))
z = data.frame(xy,dist = m[xy])
z[abs(z$dist) >= 0.7,]
## to remove because over 0.7: "BIO2","BIO3","BIO4","BIO8","BIO9","BIO10","BIO15","BIO16","BIO17"
bioclim = bioclim[,c("BIO1","BIO5","BIO6","BIO7","BIO11","BIO12","BIO13","BIO14","BIO18","BIO19")]
if(doEnvpcs==F) {
all_pca = prcomp(bioclim,center=T,scale.=T)
all_scores <- predict(all_pca)[,1:numenvpcs]
colnames(all_scores) <- paste("PC",1:numenvpcs,"A",sep="")
toprint_all = rbind(all_pca$rotation,summary(all_pca)$importance)
write.csv(toprint_all,paste("all_pca.csv",sep = ""))
}
else {
temp_pca <-
prcomp(bioclim[,c("BIO1","BIO5","BIO6","BIO7","BIO11")],center = TRUE,scale. = TRUE)
temp_scores <- predict(temp_pca)[,1:numenvpcs]
colnames(temp_scores) <- paste("PC",1:numenvpcs,"T",sep="")
prec_pca <-
prcomp(bioclim[,c("BIO12","BIO13","BIO14","BIO18","BIO19")],center = TRUE,scale. = TRUE)
summary(prec_pca)
prec_scores <- predict(prec_pca)[,1:numenvpcs]
colnames(prec_scores) <- paste("PC",1:numenvpcs,"P",sep="")
toprint_temp = rbind(temp_pca$rotation,summary(temp_pca)$importance)
toprint_prec = rbind(prec_pca$rotation,summary(prec_pca)$importance)
write.csv(toprint_temp,paste("temperature_pca.csv",sep = ""))
write.csv(toprint_prec,paste("precipitation_pca.csv",sep = ""))
}
} else {
if(doEnvpcs==F) {
all_pca = read.csv(paste("all_pca.csv",sep = ""),row.names = 1)
all_pca = all_pca[1:5,]
} else {
temp_pca = read.csv(paste("temperature_pca.csv",sep = ""),row.names = 1)
temp_pca = temp_pca[1:5,]
prec_pca = read.csv(paste("precipitation_pca.csv",sep = ""),row.names = 1)
prec_pca = prec_pca[1:5,]
temp_scores = bioclim[,c(1:5)]
prec_scores = bioclim[,c(6:10)]
for (i in 1:nrow(prec_scores)) {
row = t(prec_scores[i,])
newrow = colSums(row * prec_pca)
prec_scores[i,] = newrow
row = t(temp_scores[i,])
newrow = colSums(row * temp_pca)
temp_scores[i,] = newrow
}
temp_scores <- temp_scores[,1:numenvpcs]
colnames(temp_scores) <- paste("PC",1:numenvpcs,"T",sep="")
prec_scores = prec_scores[,1:numenvpcs]
colnames(prec_scores) <- paste("PC",1:numenvpcs,"P",sep="")
}
}
if(verbose==T){print("fin pca check")}
temp_pca = read.csv(paste("temperature_pca.csv",sep = ""),row.names = 1)
temp_pca = temp_pca[1:5,]
prec_pca = read.csv(paste("precipitation_pca.csv",sep = ""),row.names = 1)
prec_pca = prec_pca[1:5,]
temp_scores = bioclim[,c(1:5)]
prec_scores = bioclim[,c(6:10)]
for (i in 1:nrow(prec_scores)) {
row = t(prec_scores[i,])
newrow = colSums(row * prec_pca)
prec_scores[i,] = newrow
row = t(temp_scores[i,])
newrow = colSums(row * temp_pca)
temp_scores[i,] = newrow
}
temp_scores <- temp_scores[,1:numenvpcs]
colnames(temp_scores) <- paste("PC",1:numenvpcs,"T",sep="")
prec_scores = prec_scores[,1:numenvpcs]
colnames(prec_scores) <- paste("PC",1:numenvpcs,"P",sep="")
bioclim_final <-
cbind(bioclim[,1:2],temp_scores,prec_scores)
bioclim_final = bioclim_final[,3:ncol(bioclim_final)]
for (x in 2:ncol(bioclim_final)) {
bioclim_final[,x] = scales::rescale(bioclim_final[,x])
}
if(verbose==T){print("cleaning up gene data")}
data = data[data$SPP == spp,]
gene = data[,c("CATALOG.NUMBER","LAT","LONG")]
colnames(gene) = c("Sample","Latitude","Longitude")
if(class(gene$Sample) == "factor"){gene$Sample = droplevels(gene$Sample)}
gene$Sample = as.character(gene$Sample)
rownames(gene) <- make.unique(gene$Sample)
samples = gene$Sample
outbiofile="processed_bioclim_values.txt"
write.csv(x=bioclim_final,file=outbiofile)
return(list(bioclim_final,data,gene,samples))
}
generateMatrixCsvString = function(matVariable,doPCA=F,doGENE=F,doChroms=F,spp,whichpca,whichchrom,verbose=F,
basestr="~/CSVS/",
baseabun="~/CONVERTED/",
basepres=baseabun,basemid=baseabun,baselgm=baseabun,
baseibd=basestr,baseenv=basestr,basegene=basestr) {
if(matVariable=="STR") {
csvstring=paste(basestr,spp,"/",spp,"_distancematrix_STR.csv",sep = "")
} else if (matVariable=="PRES") {
csvstring=list.files(path = paste(basepres,spp,sep = ""),pattern = glob2rx("*worldclim.asc*AGGFACT-1.*MORPH-AND-GENE-DISTANCES.txt.converted.txtREDUCED.csv$"),full.names = T)[1]
} else if (matVariable=="MID") {
csvstring=list.files(path = paste(basemid,spp,sep = ""),pattern = glob2rx("*MID.asc*AGGFACT-1.*MORPH-AND-GENE-DISTANCES.txt.converted.txtREDUCED.csv$"),full.names = T)[1]
} else if (matVariable=="LGM") {
csvstring=list.files(path = paste(baselgm,spp,sep = ""),pattern = glob2rx("*LGM.asc*AGGFACT-1.*MORPH-AND-GENE-DISTANCES.txt.converted.txtREDUCED.csv$"),full.names = T)[1]
} else if (matVariable=="IBD") {
csvstring=paste(baseibd,spp,"/",spp,"_distancematrix_geog.csv",sep = "")
} else if (matVariable=="ENV") {
csvstring=paste(baseenv,spp,"/",spp,"_distancematrix_env.csv",sep = "")
} else if (matVariable=="ABUN") {
csvstring=paste(baseabun,spp,"/",spp,"_Idw_abundance_clipped.ascAGGFACT-1.MORPH-AND-GENE-DISTANCES.txt.converted.txtREDUCED.csv",sep = "")
} else if (matVariable=="GENE") {
if (doPCA == F) {
if (doGENE == F) {
csvstring = paste(basegene,spp,"/",spp,"_distancematrix_morph.csv",sep = "")
} else {
## do genes is true
if (doChroms==F) {
csvstring = paste(basegene,spp,"/",spp,"_distancematrix_FULLGENOME.csv.converted",sep = "")
} else {
csvstring = list.files(path=paste(basegene,spp,"/",sep=""),pattern=paste("Tgut_",whichchrom,".csv.converted$",sep=""),full.names = T)
}
}
} else {
csvstring = paste(basegene,spp,"/",spp,"_distancematrix_pc",whichpca,"morph.csv",sep = "")
}
} else {
csvstring=NULL
}
if(verbose==T){print(paste("LOADED FILE:",csvstring))}
return(csvstring)
}
readDataChangeHeaderGDM = function(csvstring,makeUnique=F,verbose=F) {
if(file.exists(csvstring)) {
df = (read.csv(csvstring,header = F,sep = "\t",stringsAsFactors = F))
if (ncol(df) == 1) {
df = read.csv(csvstring,header = F,sep = ",",stringsAsFactors = F)
}
## make sure only using the first one
colnames(df) = NULL
colnames(df) = as.character(df[1,])
splits = strsplit(as.character(df[1,]),"/")
firstsplits = sapply(splits,FUN=function(x){return(x[1])})
firstsplits[is.na(firstsplits)] = ""
colnames(df) = firstsplits
df[1,] = firstsplits
colnames(df)[1] = "Sample"
df[,1] = factor(df[,1])
df = unique(df)
df = df[which(df[,1] %in% colnames(df)),]
colnames(df) = make.unique(colnames(df))
rownames(df) = make.unique(as.character(df[,1]))
if(makeUnique==T) {
df = as.data.frame(unique(df))
df = as.data.frame(t(unique(df)))
df$Sample = rownames(df)
df = df[which(rownames(df) %in% c("Sample",colnames(df))),]
df = as.data.frame(t(unique(df))) ## bad 
df = as.data.frame(unique(df))
df = df[which(rownames(df) %in% c("Sample",colnames(df))),]
df = df[order(rownames(df)),order(colnames(df))]
df = df[-c(which(rownames(df)=="Sample")),-c(which(colnames(df)=="Sample"))]
df = cbind(rownames(df),df)
colnames(df)[1] = "Sample"
rownames(df) = colnames(df)[2:ncol(df)]
} else {
df = df[-1,]
df = df[order(rownames(df)),order(colnames(df))]
df = df[,-c(which(colnames(df)=="Sample"))]
df = cbind(rownames(df),df)
colnames(df)[1] = "Sample"
}
rownames(df) = make.unique(stringr::str_replace_all(rownames(df),"[^[:alnum:]]",""))
colnames(df) = make.unique(stringr::str_replace_all(colnames(df),"[^[:alnum:]]",""))
df$Sample = stringr::str_replace_all(df$Sample,"[^[:alnum:]]","")
return(df)
} else {
print("FILE DOES NOT EXIST")
return(NULL)
}
}
generateFormatted = function(univariate=F,columnnumber=NULL,bioclim_final,distPredList,GENE,dogeo=F,IBD,verbose=F) {
if(verbose==T){
print("Colnames")
print(colnames(bioclim_final))
print("Check dimensions of data")
for (dist in distPredList) {
print(dim(dist))
print(colnames(dist))
}
print("gene")
print(dim(GENE))
print(colnames(GENE))
}
if (univariate==F) {
if(verbose==T){
print("MULTI")
}
print(duplicated(GENE$Sample))
print("----")
print(duplicated(bioclim_final$Sample))
formatted <- gdm::formatsitepair(bioData = GENE,bioFormat = 3,siteColumn = "Sample",XColumn = "Longitude",YColumn = "Latitude",
predData = bioclim_final,distPreds = distPredList)
if(verbose==T){print("DIST")}
formatted$distance = scales::rescale(formatted$distance)
selected_col = "FULL"
} else {
i=columnnumber
if (i <= 6) {
if(verbose==T){print(i)}
selected_col = colnames(bioclim_final)[i+1]
if(verbose==T){print(paste("Selected Column:",selected_col))}
bioclim_subset = bioclim_final[,c("Sample",selected_col,"Latitude","Longitude")]
if(verbose==T){print("Start geo")}
if (dogeo==F) {
if(verbose==T){print("no geo")}
formatted <- gdm::formatsitepair(bioData = GENE,bioFormat = 3,siteColumn = "Sample",XColumn = "Longitude",YColumn = "Latitude",
predData = bioclim_subset)
} else {
if(verbose==T){print("yes geo")}
formatted <- gdm::formatsitepair(bioData = GENE,bioFormat = 3,siteColumn = "Sample",XColumn = "Longitude",YColumn = "Latitude",
predData = bioclim_subset,distPreds = list(as.matrix(IBD)))
}
if(verbose==T){print("format")}
formatted$distance = scales::rescale(formatted$distance)
} else {
todolist = distPredList
selected_col = names(todolist)[i-6]
if(dogeo==T) {
ibdnum = which(names(todolist)=="IBD")
todolist = todolist[unique(c(ibdnum,i-6))]
} else {
todolist = todolist[unique(c(i-6))]
}
bioclim_subset = bioclim_final[,c("Sample","Latitude","Longitude")]
formatted <- gdm::formatsitepair(bioData = GENE,bioFormat = 3,siteColumn = "Sample",XColumn = "Longitude",YColumn = "Latitude",
predData = bioclim_subset,
distPreds = todolist)
formatted$distance = scales::rescale(formatted$distance)
}
}
return(list(formatted,selected_col))
}
createModels = function(formatted,univariate,dogeo,verbose=F,prefix,suffix) {
numvars=((ncol(formatted)-6)/2)
if(verbose==T){print(numvars)}
if (numvars > 2) {
geotag = F; fulltag="full"
} else {
geotag = dogeo;fulltag="selected column"
}
if(verbose==T){
print("formatted head")
print(head(formatted))
print(formatted$distance)
print("before dimensions:"); print(dim(formatted))
}
formatted = formatted[complete.cases(formatted),]
if(verbose==T){print("after dimensions:"); print(dim(formatted))}
model <- gdm::gdm(formatted,geo = geotag,splines = NULL,knots = NULL)
if(is.null(model)) {
model <- gdm::gdm(formatted,geo = F,splines = rep(10,numvars),knots = NULL) 
counter=-1
} else {
counter=-2 
}
while(is.null(model)) {
counter = counter + 1
model <- gdm::gdm(formatted,geo = F,knots = rep(counter,3*numvars))
if(counter >= 100) {
break 
}
}
if(verbose==T){print(paste("Final Counter (-2 = default, -1 = 10 splines):",counter))}
print(paste("Explained by:",fulltag,model$explained)) 
print(paste("Null deviance: ",model$nulldeviance))
print(paste("GDM deviance: ",model$gdmdeviance)) 
if(verbose==T){print("Model splines:")
print(model$splines)}
return(model)
}
doImportanceTesting = function(doImp,model,formatted,distPredList,doPCA,doGENE,spp,doChroms,
whichchrom,verbose=F,prefix,fullmod=T,nPerm=50,doEnvpcs,numenvpcs,dogeo) {
start_time <- Sys.time()
if(verbose==T){print(start_time)}
if(model$splines[1] == 3) {
modTest <- gdm.varImp(formatted,geo = F,fullModelOnly = fullmod,nPerm=nPerm)
} else {
modTest <- gdm.varImp(formatted,geo = F,fullModelOnly = T,nPerm=nPerm)
}
end_time <- Sys.time()
if(verbose==T){print(end_time - start_time)}
## check if it worked
print(head(modTest))
if(is.null(modTest)){null_or_sum=T} else if(sum(modTest[[2]][,1],na.rm=T) <= 0) {null_or_sum=T} else {null_or_sum=F}
if (null_or_sum) {
if(verbose==T){print("Running Fix 1")}
modTest <- gdm.varImp(formatted,geo = F,fullModelOnly = F,nPerm=nPerm)
if(is.null(modTest)){null_or_sum=T} else if(sum(modTest[[2]][,1],na.rm=T) <= 0) {null_or_sum=T} else {null_or_sum=F}
}
if (null_or_sum) {
if(verbose==T){print("Running Fix 2")}
modTest <- gdm.varImp(formatted,geo = F,fullModelOnly = T,nPerm=nPerm)
if(is.null(modTest)){null_or_sum=T} else if(sum(modTest[[2]][,1],na.rm=T) <= 0) {null_or_sum=T} else {null_or_sum=F}
}
if (null_or_sum) {
print("MODEL FAILED TO CONVERGE")
numvars=((ncol(formatted)-6)/2)
modTest <- gdm.varImp.edit(formatted,geo = F,fullModelOnly = T,nPerm=nPerm,
splines = rep(10,numvars),knots = NULL) 
counter=-1 
if(verbose==T){print("TESTING COUNTERS")}
if(is.null(modTest)){null_or_sum=T} else if(sum(modTest[[2]][,1],na.rm=T) <= 0) {null_or_sum=T} else {null_or_sum=F}
if(verbose==T){print("END COUNTERS")}
print(null_or_sum)
while(null_or_sum) {
counter = counter + 1
if(verbose==T){print(counter)}
modTest <- gdm.varImp.edit(formatted,geo = F,fullModelOnly = T,nPerm=nPerm,
knots = rep(counter,3*numvars))
if(is.null(modTest)){null_or_sum=T} else if(sum(modTest[[2]][,1],na.rm=T) <= 0) {null_or_sum=T} else {null_or_sum=F}
if(counter >= 100) {
print("MODEL FAILED TO CONVERGE AGAIN")
break 
}
if(is.null(modTest)){null_or_sum=T} else if(sum(modTest[[2]][,1],na.rm=T) <= 0) {null_or_sum=T} else {null_or_sum=F}
}
if(verbose==T){print("SEEING COUNTER")}
if(verbose==T){print(counter)}
if(counter == -1) {
modTest <- gdm.varImp.edit(formatted,geo = F,fullModelOnly = F,nPerm=nPerm,
splines = rep(10,numvars),knots = NULL) 
} else if (counter < 100) {
modTest <- gdm.varImp.edit(formatted,geo = F,fullModelOnly = F,nPerm=nPerm,
knots = rep(counter,3*numvars))
} else {
return()
}
if(verbose==T){print("DONE COUNTERS")}
}
if(doEnvpcs==T) {
start=(numenvpcs*2)+1
} else {
start=2
}
if(verbose==T){print("done envpcs")}
rownames(modTest[[2]])[(start):length(model$predictors)] = names(distPredList)
rownames(modTest[[3]])[(start):length(model$predictors)] = names(distPredList)
rownames(modTest[[4]])[(start):length(model$predictors)] = names(distPredList)
if (doPCA == F) {
if (doGENE == F) {
suffix="_allmorph" 
} else {
if(doChroms==F) {
suffix="_gene" } else {suffix=paste("_chr",whichchrom,sep="") }
}
} else {
suffix=paste("_pc",whichpca,sep = "") 
}
if(dogeo==T) {
suffix=paste("_geo",suffix,sep="")
}
if(verbose==T){print("done suffixes")}
if(!(is.null(modTest))) {
if(verbose==T){print("start deviances")}
#Deviance
print(head(modTest[[1]]))
print(paste(spp,"_",prefix,"_variable_deviance",suffix,".csv",sep = ""))
write.csv(modTest[[1]],paste(spp,"_",prefix,"_variable_deviance",suffix,".csv",sep = ""))
if(verbose==T){print("start importance")}
#Importance
modTest[[2]]
write.csv(modTest[[2]],paste(spp,"_",prefix,"_variable_importance",suffix,".csv",sep = ""))
if(verbose==T){print("start signifs")}
#Significance
modTest[[3]]
write.csv(modTest[[3]],paste(spp,"_",prefix,"_variable_significance",suffix,".csv",sep = ""))
#Permuts Calculated
if(verbose==T){print("start permutes")}
modTest[[4]]
write.csv(modTest[[4]],paste(spp,"_",prefix,"_variable_permutations",suffix,".csv",sep = ""))
if(verbose==T){print("IMPORTANCE PLOT")}
png(paste(spp,"_",prefix,"_importance",suffix,".png",sep = ""),height = 500,width = 1200)
if(model$splines[1] == 3) {
par(mfrow = c(3,4)) 
}
barplot(modTest[[2]][,1] / colSums(modTest[[2]],na.rm = T)[1],main = "full Model",las = 2,ylim=c(0,1))
if(model$splines[1] == 3) {
lapply(2:ncol(modTest[[2]]),FUN = function(x) {barplot(modTest[[2]][,x] / colSums(modTest[[2]],na.rm = T)[x], main = paste("full Model -",x - 1), las = 2 ,ylim=c(0,1))})
}
dev.off()
}
if(verbose==T){print("done importance testing")}
}
runGDMs = function(doimp = F, doPCA = F,doGENE = F,whichpca = 1,univariate = T,dogeo = T,doChroms = F,whichchrom="Z",
matList=c("STR","PRES","MID","LGM","IBD","ENV","ABUN","GENE"),
spplist=c("BELLII","BILINEATA","BRUNNEICAPILLUS","CRISSALE","CURVIROSTRE","FLAVICEPS","FUSCA","MELANURA","NITENS","SINUATUS"),
doEnvpcs=T,verbose=F,numenvpcs,overwrite=F,genedistmatrix=NULL) {
print(paste("Do Importance?",doimp,"Do Gene?",doGENE,"Do Pca?",doPCA,"Which PCA?",whichpca))
print(paste("Do univariate?",univariate,"Do IBD?",dogeo,"Do chromosomes?",doChroms,"Which chroms?"))
print(whichchrom)
prefix = paste(matList[2:length(matList)],collapse="-")
if(doEnvpcs==T){
prefix=paste("PC",numenvpcs,"-",prefix,sep="")
}
if (doPCA == F) {
if (doGENE == F) {
suffix="_allmorph" 
} else {
if(doChroms==F) {
suffix="_gene" } else {suffix=paste("_chr",whichchrom,sep="") }
}
} else {
suffix=paste("_pc",whichpca,sep = "") 
}
if(dogeo==T) {
suffix=paste("_geo",suffix,sep="")
}
for (spp in spplist) {
print(spp)
modelstring2 = paste(spp,"_",prefix,"_modelparameters",suffix,".txt",sep = "") 
if(file.exists(modelstring2) && overwrite==F) {
print("SKIPPING")
} else {
if(verbose==T){print("loading morph")}
morphlist = loadMorphData(spp,doGENE)
morphdf = morphlist[[1]]
morphlocs = morphlist[[2]]
morphspp = morphlist[[3]]
setwd("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/")
path1="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/"
path2=paste(path1,"GDM_results/",sep="")
path3=paste(path2,"univariate/",sep="")
path4=paste(path2,"biivariate/",sep="")
path5=paste(path2,"multivariate/",sep="")
if(verbose==T){print("loading environment")
print(head(morphlocs))
print(head(morphspp))
print(head(morphdf))
print(doGENE)
print(spp)
print(verbose)
print(numenvpcs)
}
envlist = loadEnvData(morphlocs,morphspp,morphdf,doGENE,spp,verbose=verbose,numenvpcs=numenvpcs,doEnvpcs=doEnvpcs)
if(verbose==T){print("done env")}
bioclim_final = envlist[[1]]
data = envlist[[2]]
gene = envlist[[3]]
samples = envlist[[4]]
gsplit = strsplit(as.character(gene$Sample),"/")
gene$Sample = sapply(gsplit,FUN=function(x){return(x[1])})
gene=unique(gene)
if(verbose==T){print("relabeling")}
colnames(bioclim_final) = (stringr::str_replace_all(colnames(bioclim_final),"[^[:alnum:]]",""))
bioclim_final$Sample = (stringr::str_replace_all(rownames(bioclim_final),"[^[:alnum:]]",""))
rownames(bioclim_final) = (stringr::str_replace_all(rownames(bioclim_final),"[^[:alnum:]]",""))
bsplit = strsplit(as.character(bioclim_final$Sample),"/")
bioclim_final$Sample = sapply(bsplit,FUN=function(x){return(x[1])})
bioclim_final$Sample = stringr::str_replace_all(bioclim_final$Sample,"[^[:alnum:]]","")
bioclim_final=unique(bioclim_final)
rownames(gene) = stringr::str_replace_all(rownames(gene),"[^[:alnum:]]","")
colnames(gene) = stringr::str_replace_all(colnames(gene),"[^[:alnum:]]","")
gene$Sample = stringr::str_replace_all(gene$Sample,"[^[:alnum:]]","")
bioclim_final = merge(bioclim_final,gene)
rownames(bioclim_final) = make.unique(bioclim_final$Sample)
bioclim_final$Sample = stringr::str_replace_all(bioclim_final$Sample,"[^[:alnum:]]","")
if(doEnvpcs==F) {
bioclim_final$ZERO = rep(0,nrow(bioclim_final))
bioclim_final = bioclim_final[,c("Sample","ZERO","Latitude","Longitude")]
}
## load matList 
if(verbose==T){print("loading matlist")}
toexit=F
for (matVariable in matList) {
print(matVariable)
if(is.null(genedistmatrix) | matVariable != "GENE") {
csvstring = generateMatrixCsvString(matVariable,doPCA,doGENE,doChroms,spp,whichpca,whichchrom,verbose=verbose)
} else {
if(verbose==T){print(paste("IMPORTING GENEDISTMATRIX"))}
csvstring = genedistmatrix
} 
if(verbose==T){print(paste("LOADED FILE:",csvstring))}
if(sum(!file.exists(csvstring))!=0){
toexit=T
print("FILE DOES NOT EXIST, SKIPPING")
print(csvstring)
} else {
if(length(csvstring) > 1) {
df = lapply(csvstring,FUN=function(x){readDataChangeHeaderGDM(x,verbose=verbose)})
df = df[order(rownames(df)),order(colnames(df))]
assign(matVariable,df)
} else if (length(csvstring) == 1) {
df = readDataChangeHeaderGDM(csvstring,makeUnique = T,verbose=verbose) 
df = df[order(rownames(df)),order(colnames(df))]
assign(matVariable,df)
} else {
print("FILE DOES NOT EXIST, SKIPPING")
toexit = T
}
}
}
if(toexit==T) {
next
}
if(verbose==T){print("cross-referencing matlists")}
# generate the ones to keep
tokeep = c()
missing = c()
originals=c()
allpossibleones=c()
for (matVariable in c(matList,"bioclim_final")) {
if(verbose==T){print(matVariable)
print("Samples here:")}
tempsamp = sort(get(matVariable)$Sample)
if(verbose==T){print(tempsamp)}
allpossibleones =c(allpossibleones,get(matVariable)$Sample)
if (is.null(tokeep)) { 
tokeep = sort(intersect(get(matVariable)$Sample,get(matVariable)$Sample))
originals=tokeep
} else { 
tempmiss = tokeep
tokeep = intersect(sort(tokeep),get(matVariable)$Sample) 
missing = c(missing,tempmiss[which(!(tempmiss %in% tokeep))])
}
if(verbose==T){print(paste("temp keep:",length(tokeep)))
print(sort(tokeep))}
}
tokeep = c("Sample",tokeep)
allpossibleones = sort(unique(allpossibleones))
if(verbose==T){print(paste("Keeping:",length(tokeep),"of",length(originals)+1))}
if(verbose==T){print(tokeep)}
if(verbose==T){
if(length(missing) > 0) {
print("ORIGNALS:")
print(sort(originals))
print("MISSING ONES:")
print(sort(missing))
print("CHECK IF MISSING ARE IN ALL POSSIBLE:")
print(allpossibleones[which((allpossibleones %in% missing))])
}
}
## remove the ones that don't match 
## MAKE SURE YOU DON'T GET RID OF SAMPLE
for (matVariable in matList) {
if(verbose==T){print(matVariable)}
toset = get(matVariable)
toset = toset[order(rownames(toset)),order(colnames(toset))]
samplecolumn = which(colnames(toset) == "Sample")
if(setequal(rownames(toset),colnames(toset)[colnames(toset)!="Sample"])) {
if(setequal(toset$Sample,colnames(toset)[colnames(toset)!="Sample"])) {
toset = toset[which(toset$Sample %in% tokeep),
which(colnames(toset) %in% tokeep)]
} else {
if(verbose==T){print("SAMPLE AND COLS DON'T MATCH")}
sam_keep = length(intersect(tokeep,toset$Sample))
col_keep = length(intersect(tokeep,colnames(toset)[colnames(toset)!="Sample"]))
if (sam_keep >= col_keep) {
toset = toset[which(toset$Sample %in% c(tokeep,"Sample")),
c(samplecolumn,which(toset$Sample %in% c(tokeep,"Sample")))]
} else {
toset = toset[which(colnames(toset) %in% c(tokeep,"Sample")),
c("Sample",which(colnames(toset) %in% c(tokeep,"Sample")))]
}
}
} else {
if(verbose==T){print("ROWS AND COLS DON'T MATCH")}
}
assign(matVariable,toset)
}
if(verbose==T){print("done removing matches")}
bioclim_final = bioclim_final[which(bioclim_final$Sample %in% tokeep),]
distPredList = lapply(matList,FUN = function(x) {df = get(matVariable); df = as.matrix(df); return(df) })
numallnames = lapply(distPredList,FUN = function(x) {return(length(x[,1])) })
numunqnames = lapply(distPredList,FUN = function(x) {return(length(unique(x[,1]))) })
duplistname = lapply(distPredList,FUN = function(x) {return(duplicated(unique(x[,1]))) })
names(distPredList) = matList
distPredList = distPredList[names(distPredList)!="GENE"]
if(verbose==T){print("done making list")}
print(matList)
print(numallnames)
print(numunqnames)
print(duplistname)
if(univariate==F) {
if(verbose==T){print("start multivar")}
formattedlist = generateFormatted(univariate,columnnumber=NULL,bioclim_final,distPredList,GENE,dogeo,IBD,verbose=verbose)
if(verbose==T){print("end format")}
formatted = formattedlist[[1]]
formatted$distance = scales::rescale(formatted$distance)
oldrows = nrow(formatted)
selected_col = formattedlist[[2]]
## run models 
if(verbose==T){print("start models")}
model = createModels(formatted,univariate,dogeo,verbose=verbose,prefix=prefix,suffix=suffix)
if(verbose==T){print("end models")}
if(verbose==T){print("start null check")}
if(!is.null(model)) {
isplines = gdm::isplineExtract(model)
splineheight = isplines$y[nrow(isplines$y),]
if(verbose==T){print("not null")}
} else {
print("FINAL MODEL IS BAD, STOP")
break
}
if(verbose==T){print("end null check")}
if(verbose==T){print("output params file")}
capture.output(summary(model),file=modelstring2)
if(verbose==T){print("calculate start")}
if(doEnvpcs==T) {
start=(numenvpcs*2)+1
} else {
if(dogeo==F) {start=2} else {start=3}
}
print(splineheight)
print(names(distPredList))
if(length(names(distPredList))==1 && names(splineheight)[length(splineheight)] == "matrix_1") {
names(splineheight)[length(splineheight)] = names(distPredList)
} else {
names(splineheight)[start:length(splineheight)] = names(distPredList)
}
print("Spline heights:")
print(splineheight)
heightstring = paste(spp,"_",prefix,"_splineheights",suffix,".png",sep = "") 
png(heightstring)
barplot(splineheight,las=2,main=spp)
dev.off()
heightstring2 = paste(spp,"_",prefix,"_splineheights",suffix,".txt",sep = "") 
write.csv(splineheight,heightstring2)
if (!(is.null(model))) {
if(doEnvpcs==T) {modeloff=(numenvpcs*2)+1} else {modeloff=2}
model$predictors[modeloff:length(model$predictors)] = names(distPredList)
splinestring = paste(spp,"_",prefix,"_splines",suffix,".png",sep = "") 
png(splinestring,height = 900,width = 1500) 
if(modeloff==7) {
plot(model,plot.layout = c(2,6))
} else {
plot(model,plot.layout = c(2,4))
}
dev.off()
}
if(doimp==T) {
doImportanceTesting(doImp,model,formatted,distPredList,doPCA,doGENE,spp,doChroms,whichchrom,prefix=prefix,verbose=verbose,doEnvpcs=doEnvpcs,numenvpcs=numenvpcs,dogeo=dogeo)
}
} else {
if(doEnvpcs==T) {colsep=(numenvpcs*2)} else {colsep=0}
for (columnnumber in 1:(colsep+length(distPredList))) {
formattedlist = generateFormatted(univariate,columnnumber,bioclim_final,distPredList,GENE,dogeo,IBD,verbose=verbose)
formatted = formattedlist[[1]]
selected_col = formattedlist[[2]]
if(verbose==T){print("make uni model")}
model = createModels(formatted,univariate,dogeo,verbose=verbose,prefix=prefix,suffix=suffix)
if (!(is.null(model))) {
length(model$predictors)
if(dogeo==T) {geotext="IBD"} else {geotext=""}
pngtext = paste(substr(spp,1,3),"_",prefix,"_",selected_col,"_splines",geotext,suffix,".png",sep = "")
png(pngtext,height = 450,width = 600) 
plot(model,plot.layout = c(1,5))
mtext(text=paste(spp,selected_col,"Expl:",(model$explained)))
dev.off()
}
}
}
}
}
}
matList=c("GENE","STR","IBD","PRES","ABUN","ENV","LGM")
combinations = c(combn(2:length(matList),1,simplify = F),
combn(2:length(matList),2,simplify = F),
combn(2:length(matList),3,simplify = F)
)
genemat_list_full = list.files(pattern="ngsdist",full.names = T)
for(spp in spplist){
genemat_list_spp = genemat_list_full[grepl(spp,genemat_list_full)]
for(chr in chromlist){
genemat_list = genemat_list_spp[grepl(chr,genemat_list_spp)]
for(genedistmatrix_i in 1:length(genemat_list)){
genedistmatrix = genemat_list[genedistmatrix_i]
spp = genemat_spp[genedistmatrix_i]
runGDMs(doimp = F,doPCA = F,doGENE = T,whichpca = 1,
univariate = F,dogeo = F,doChroms = T,
whichchrom=chr,matList=c("GENE","ABUN"),spplist=spp,
doEnvpcs=F,
verbose=F,numenvpcs=3,genedistmatrix=genedistmatrix,overwrite=F)
runGDMs(doimp = F,doPCA = F,doGENE = T,whichpca = 1,
univariate = F,dogeo = F,doChroms = T,
whichchrom=chr,matList=c("GENE","ENV"),spplist=spp,
doEnvpcs=F,
verbose=F,numenvpcs=3,genedistmatrix=genedistmatrix,overwrite=F)
runGDMs(doimp = F,doPCA = F,doGENE = T,whichpca = 1,
univariate = F,dogeo = F,doChroms = T,
whichchrom=chr,matList=c("GENE","IBD"),spplist=spp,
doEnvpcs=F,
verbose=F,numenvpcs=3,genedistmatrix=genedistmatrix,overwrite=F)
runGDMs(doimp = F,doPCA = F,doGENE = T,whichpca = 1,
univariate = F,dogeo = F,doChroms = T,
whichchrom=chr,matList=c("GENE","LGM"),spplist=spp,
doEnvpcs=F,
verbose=F,numenvpcs=3,genedistmatrix=genedistmatrix,overwrite=F)
runGDMs(doimp = F,doPCA = F,doGENE = T,whichpca = 1,
univariate = F,dogeo = F,doChroms = T,
whichchrom=chr,matList=c("GENE","PRES"),spplist=spp,
doEnvpcs=F,
verbose=F,numenvpcs=3,genedistmatrix=genedistmatrix,overwrite=F)
runGDMs(doimp = F,doPCA = F,doGENE = T,whichpca = 1,
univariate = F,dogeo = F,doChroms = T,
whichchrom=chr,matList=c("GENE","STR"),spplist=spp,
doEnvpcs=F,
verbose=F,numenvpcs=3,genedistmatrix=genedistmatrix,overwrite=F)
runGDMs(doimp = F,doPCA = F,doGENE = T,whichpca = 1,
univariate = F,dogeo = F,doChroms = T,
whichchrom=chr,matList=c("GENE","STR","ABUN"),spplist=spp,
doEnvpcs=F,
verbose=F,numenvpcs=3,genedistmatrix=genedistmatrix,overwrite=F)
runGDMs(doimp = F,doPCA = F,doGENE = T,whichpca = 1,
univariate = F,dogeo = F,doChroms = T,
whichchrom=chr,matList=c("GENE","STR","ENV"),spplist=spp,
doEnvpcs=F,
verbose=F,numenvpcs=3,genedistmatrix=genedistmatrix,overwrite=F)
runGDMs(doimp = F,doPCA = F,doGENE = T,whichpca = 1,
univariate = F,dogeo = F,doChroms = T,
whichchrom=chr,matList=c("GENE","STR","IBD"),spplist=spp,
doEnvpcs=F,
verbose=F,numenvpcs=3,genedistmatrix=genedistmatrix,overwrite=F)
runGDMs(doimp = F,doPCA = F,doGENE = T,whichpca = 1,
univariate = F,dogeo = F,doChroms = T,
whichchrom=chr,matList=c("GENE","STR","LGM"),spplist=spp,
doEnvpcs=F,
verbose=F,numenvpcs=3,genedistmatrix=genedistmatrix,overwrite=F)
runGDMs(doimp = F,doPCA = F,doGENE = T,whichpca = 1,
univariate = F,dogeo = F,doChroms = T,
whichchrom=chr,matList=c("GENE","STR","PRES"),spplist=spp,
doEnvpcs=F,
verbose=F,numenvpcs=3,genedistmatrix=genedistmatrix,overwrite=F)
runGDMs(doimp = F,doPCA = F,doGENE = T,whichpca = 1,
univariate = F,dogeo = F,doChroms = T,
whichchrom=chr,matList=c("GENE","STR","IBD","ABUN"),spplist=spp,
doEnvpcs=F,
verbose=F,numenvpcs=3,genedistmatrix=genedistmatrix,overwrite=F)
runGDMs(doimp = F,doPCA = F,doGENE = T,whichpca = 1,
univariate = F,dogeo = F,doChroms = T,
whichchrom=chr,matList=c("GENE","STR","IBD","ENV"),spplist=spp,
doEnvpcs=F,
verbose=F,numenvpcs=3,genedistmatrix=genedistmatrix,overwrite=F)
runGDMs(doimp = F,doPCA = F,doGENE = T,whichpca = 1,
univariate = F,dogeo = F,doChroms = T,
whichchrom=chr,matList=c("GENE","STR","IBD","LGM"),spplist=spp,
doEnvpcs=F,
verbose=F,numenvpcs=3,genedistmatrix=genedistmatrix,overwrite=F)
runGDMs(doimp = F,doPCA = F,doGENE = T,whichpca = 1,
univariate = F,dogeo = F,doChroms = T,
whichchrom=chr,matList=c("GENE","STR","IBD","PRES"),spplist=spp,
doEnvpcs=F,
verbose=F,numenvpcs=3,genedistmatrix=genedistmatrix,overwrite=F)
}
}
}