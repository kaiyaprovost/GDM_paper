library(subsppLabelR)
library(ecospat)
library(ENMTools)
library(rasterExtras)
library(landsim)
library(sys)
library(scales)
dynamic_require <- function(package) {
  if (eval(parse(text = paste("require(", package, ")"))))
    return(TRUE)
  install.packages(package)
  return(eval(parse(text = paste(
    "require(", package, ")"
  ))))
}
packages = c(
  "dismo",
  "GISTools",
  "raster",
  "rgdal",
  "ENMeval",
  "phyloclim",
  "data.table",
  "dplyr",
  "EMCluster",
  "knor",
  "maps",
  "MASS",
  "parallel",
  "plotly",
  "rgeos",
  "roxygen2",
  "rworldmap",
  "sp",
  "spThin",
  "spocc",
  "spThin",
  "viridis",
  "auk",
  "rebird"
)
for (p in packages) {
  dynamic_require(p)
}
path = "~/WORLDCLIM/"
setwd(path)
specieslist=c(
  "Amphispiza bilineata",
  "Auriparus flaviceps",
  "Campylorhynchus brunneicapillus",
  "Cardinalis sinuatus",
  "Melozone fusca",
  "Phainopepla nitens",
  "Polioptila melanura",
  "Toxostoma crissale",
  "Toxostoma curvirostre",
  "Vireo bellii"
)
## stack time slices
for (spp in rev(specieslist)) {
  print(spp)
  current_f = list.files(path=path,pattern=paste("Thresh_EqualSensSpec_",spp,"_worldclim.asc$",sep=".+"),recursive = T,full.names = T)[1]
  current = raster (current_f)
  LGM_f = list.files(path=path,pattern=paste("Thresh_EqualSensSpec_",spp,"_LGM.asc$",sep=".+"),recursive = T,full.names = T)[1]
  LGM = raster (LGM_f)
  MID_f = list.files(path=path,pattern=paste("Thresh_EqualSensSpec_",spp,"_MID.asc$",sep=".+"),recursive = T,full.names = T)[1]
  MID = raster (MID_f)
  MID_2 = MID
  values(MID_2)[values(MID_2)==1 & !(is.na(values(MID_2)))] = 2
  LGM_4 = LGM
  values(LGM_4)[values(LGM_4)==1 & !(is.na(values(LGM_4)))] = 4
  summed = current+MID_2+LGM_4
  png(paste("Thresh_Timestack_",spp,"_difcolors.png",sep=""))
  plot(summed,col=c(
    "black", #//
    "#FF7F00", #C
    "#FB9A99", #M
    "#E31A1C", #CM
    "#1F78B4", #L
    "#6A3D9A", #CL
    "#CAB2D6", #ML
    "grey" #CLM
  ))
  legend(legend=c("None","Now-Only","Mid-Only","Now-Mid",
                  "LGM-Only","Now-LGM","Mid-LGM","All"),
         x="right",
         title=spp,
         bty="n",
         fill=c("black", #//
                "#FF7F00", #C
                "#FB9A99", #M
                "#E31A1C", #CM
                "#1F78B4", #L
                "#6A3D9A", #CL
                "#CAB2D6", #ML
                "grey" #CLM
         ))
  dev.off()
}
