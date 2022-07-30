## CHANGE THIS TO DETERMINE WHERE TO START
FIRSTSTEP = 8

#detach("package:subsppLabelR", unload = TRUE)
#library(devtools)
#devtools::install_github('kaiyaprovost/subsppLabelR', force = F)
# devtools::install_github("cran/ecospat")
# devtools::install_github("danlwarren/ENMTools")
# devtools::install_github("rsh249/rasterExtras")
# devtools::install_github("petrelharp/landsim")
library(subsppLabelR)
library(ecospat)
library(ENMTools)
library(rasterExtras)
library(landsim)
library(sys)

#sudo R CMD javareconf JAVA_HOME=/System/Library/Frameworks/JavaVM.framework/Versions/A/Headers/

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
## rJava is causing issues

for (p in packages) {
  dynamic_require(p)
}

#####

if (FIRSTSTEP <= 1) {

bandwidth = 2
#a value in kilometers
nclus = 24
options(java.parameters = "-Xmx1g")
resam = 100
samsize = 50000

Env = raster::stack(
  list.files(
    path = '/Users/kprovost/Dropbox (AMNH)/Classes/Spatial Bioinformatics/spatial_bioinformatics-master/ENM/wc2-5/',
    #path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/bio_2-5m_bil/",
    pattern = "\\.bil$",
    full.names = T
  )
)

LGM = raster::stack(
  list.files(
    path = '/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/past_climates/cclgmbi_2-5m/',
    pattern = "\\.tif$",
    full.names = T
  )
)

MID = raster::stack(
  list.files(
    path = '/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/past_climates/ccmidbi_2-5m/',
    pattern = "\\.tif$",
    full.names = T
  )
)

ext = raster::extent(c(-125,-60, 10, 50)) ## for cropping the layers for backproject
#ext = raster::extent(c(-115,-97,26,37)) ## for making the models final layers
Env = raster::crop(Env, ext)
LGM = raster::crop(LGM, ext)
MID = raster::crop(MID, ext)


bg = Env[[1]] ## just for plotting
## in case you need to redo the worldclim:
# writeRaster(Env,filename="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/ENMS_multilayer_onlyworldclim.tif",
#             options="INTERLEAVE=BAND", overwrite=TRUE)
# writeRaster(LGM,filename="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/ENMS_multilayer_LGM.tif",
#             options="INTERLEAVE=BAND", overwrite=TRUE)
# writeRaster(MID,filename="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/ENMS_multilayer_MID.tif",
#             options="INTERLEAVE=BAND", overwrite=TRUE)

## step 1: import data -- with subspecies this time


  print("STEP1")
  list_of_taxa = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/southwestern_subspecies.txt"
  taxa = read.csv(list_of_taxa, sep = "\t", header = F)
  colnames(taxa) = c("GEN", "SPP", "SUB")
  unique_species = unique(taxa[, 1:2])
  
  ## just extract Phainopepla
  #uni_phaino = unique_species[unique_species$GEN=="Phainopepla",]
  
  outputDir = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/OUTPUTS/"
  setwd(outputDir)
  
  for (rownum in 6:nrow(unique_species)) {
    temp = taxa[taxa$GEN == unique_species$GEN[rownum],]
    temp = temp[temp$SPP == unique_species$SPP[rownum],]
    subspp = droplevels.factor(unlist(temp$SUB))
    spp = paste((unlist(temp[1, 1:2])), sep = " ", collapse = " ")
    print(spp)
    print(subspp)
    
    print("Getting Subspecies Information")
    pointLimit = 1000
    dbToQuery = c("gbif", "bison", "ecoengine", "vertnet") ## removed everything besides gbif and vertnet?
    
    outputDir_spp = paste(outputDir, spp, sep = "")
    dir.create(outputDir_spp)
    setwd(outputDir_spp)
    
    listFromSubspeciesOcc = subspeciesOccQuery(
      spp = spp,
      subsppList = subspp,
      pointLimit = pointLimit,
      dbToQuery = dbToQuery
    )
    
    print("Labeling data by subspecies")
    labeledLoc = labelSubspecies(subsppOccList = listFromSubspeciesOcc)
    
    ## EXPORT THE OCCURRENCE DATA!
    ## NOTE: you don't need to do this. it is identical to merging goodsubspp and suspectsubspp and only taking first four cols
    ## (name	longitude	latitude	subspecies)
    print("Exporting Occurrence Data")
    write.table(
      labeledLoc,
      paste(
        paste(
          "OccurrenceDatabase",
          spp,
          paste(subspp, collapse = " "),
          sep = "_"
        ),
        ".occ",
        sep = ""
      ),
      quote = FALSE,
      sep = "\t",
      row.names = FALSE
    )
    
  }
}

#####

## step 2: also import the data currently on my computer

if (FIRSTSTEP <= 2) {
  print("STEP2")
  localities_ebird = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/locality information/EBIRD/ebd_trim_spplatlong_unique.txt"
  ebd = read.csv(localities_ebird, sep = "\t", header = F)
  colnames(ebd) = c("name", "latitude", "longitude")
  ebd = ebd[, c(1, 3, 2)]
  ## split by species
  out <- split(ebd , f = ebd$name)
  
  path = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/locality information/EBIRD/by_species/"
  setwd(path)
  
  for (name in names(out)) {
    temp = out[[name]]
    
    name = gsub("(", "", name, fixed = T)
    name = gsub(")", "", name, fixed = T)
    name = gsub("/", " or ", name, fixed = T)
    name = gsub(".", "", name, fixed = T)
    print(name)
    
    write.table(
      temp,
      file = paste("Ebird ", name, ".txt", sep = ""),
      col.names = colnames(ebd),
      row.names = F
    )
    
  }
  
  ## from ebird:
  ebird_focal = read.csv(
    "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/locality information/EBIRD/Ebird_focal_species.txt",
    row.names = NULL,
    sep = " "
  )
  #ebird_focal$subspecies = rep("unknown",nrow(ebird_focal))
  write.table(
    ebird_focal,
    "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/locality information/EBIRD/Ebird_focal_species.txt",
    row.names = F
  )
  
  ## gbif
  gbif = read.csv(
    "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/locality information/enms_GBIF/Localities/AllSpecies_NotThinned_NoUncertains.txt",
    sep = "\t",
    as.is = T
  )
  gbif = gbif[, c(3, 10, 9, 4)]
  colnames(gbif) = c("name", "longitude", "latitude", "subspecies")
  gbif$subspecies[gbif$subspecies == ""] = "unknown"
  
  write.table(
    gbif,
    "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/locality information/enms_GBIF/localities_compiled_for_enm.txt",
    row.names = F
  )
  
  ## tissues and songs
  
  file1 = read.csv(
    "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/locality information/tissues and points/localities/allSpecies_11July2017.csv",
    header = T,
    sep = ",",
    as.is = T
  )
  file1 = file1[, c(56, 4, 3, 61)]
  colnames(file1) = c("name", "longitude", "latitude", "subspecies")
  file1 = unique(file1)
  
  file2 = read.csv(
    "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/locality information/tissues and points/localities/testTable_unique.txt",
    sep = ",",
    header = T,
    as.is = T
  )
  file2 = file2[, c(3, 2, 1, 4)]
  colnames(file2) = c("name", "longitude", "latitude", "subspecies")
  file2$subspecies = rep("unknown", nrow(file2))
  file2 = unique(file2)
  
  file3 = read.csv(
    "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/locality information/tissues and points/localities/towhees_unique.txt",
    sep = "\t",
    as.is = T
  )
  file3 = file3[, c(1, 4, 3, 2)]
  colnames(file3) = c("name", "longitude", "latitude", "subspecies")
  file3$subspecies = rep("unknown", nrow(file3))
  file3 = unique(file3)
  
  file4 = read.csv(
    "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/locality information/tissues and points/localities/Vertnet_unique.txt",
    sep = "\t",
    header = T,
    as.is = T
  )
  file4 = file4[, c(84, 85, 141, 143, 144)]
  file4 = unique(file4)
  file4 = file4[file4$decimallatitude != "",]
  file4$name = paste(file4$genus, " ", file4$specificepithet, sep = "")
  file4 = file4[, c(6, 2, 1, 5)]
  colnames(file4) = c("name", "longitude", "latitude", "subspecies")
  file4 = file4[file4$subspecies != "taxonrank",]
  file4$subspecies[file4$subspecies == ""] = "unknown"
  
  
  merged = rbind(file1, file2, file3, file4)
  merged = unique(merged)
  merged = merged[merged$longitude != "",]
  write.table(
    merged,
    "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/locality information/tissues and points/EverythingMergedTogether.txt",
    row.names = F
  )
  
  
  every_loc = rbind(ebird_focal, gbif, merged)
  every_loc = unique(every_loc)
  every_loc = every_loc[every_loc$latitude != "",]
  
  write.table(
    every_loc,
    "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/locality information/all_source_localities_merged_together.txt",
    row.names = F
  )
  
  ## now need to split the everything by species
  out_loc <- split(every_loc , f = every_loc$name)
  for (name in sort(names(out_loc))) {
    if (name != "") {
      temp = out_loc[[name]]
      
      name = gsub("(", "", name, fixed = T)
      name = gsub(")", "", name, fixed = T)
      name = gsub("/", " or ", name, fixed = T)
      name = gsub(".", "", name, fixed = T)
      name = gsub(" ", "_", name, fixed = T)
      print(name)
      name_for_file = paste(strsplit(name, "_")[[1]][1:2],
                            sep = "_",
                            collapse = "_")
      print(name_for_file)
      
      write.table(
        temp,
        file = paste(
          "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/locality information/merged_by_species/All Locs Combined ",
          name_for_file,
          ".txt",
          sep = ""
        ),
        col.names = colnames(every_loc),
        row.names = F,
        append = T,
        sep = "\t"
      )
    }
    
  }
  
  ## now import all of the .occ files with the merged file
  occ = c()
  setwd(outputDir)
  paths = Sys.glob(file.path(outputDir, "*", "/", "*.occ"))
  for (path in paths) {
    print(path)
    temp = read.csv(path, sep = "\t")
    print(colnames(temp))
    if (is.null(occ)) {
      occ = temp
      
    } else {
      occ = gtools::smartbind(occ, temp)
    }
    print(nrow(occ))
  }
  
  occ = unique(occ)
  nrow(occ)
  small_occ = occ[, c("name", "longitude", "latitude", "subspecies")]
  small_occ = unique(small_occ)
  nrow(small_occ)
  
  write.table(
    small_occ,
    "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/OUTPUTS/all_species.occ",
    row.names = F
  )
  
  ## import the merged new localities
  merged = read.csv(
    "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/locality information/all_source_localities_merged_together.txt",
    sep = "\t"
  )
  final_merge = gtools::smartbind(small_occ, merged)
  final_merge$name = gsub("_", " ", final_merge$name)
  final_merge = unique(final_merge)
  nrow(final_merge)
  
  final_merge$genus = final_merge$name
  final_merge$species = final_merge$name
  
  for (row in 1:nrow(final_merge)) {
    if (row %% 1000 == 0) {
      print(paste(row, nrow(final_merge), sep = " / "))
    }
    #print(final_merge$name[row])
    split = unlist(strsplit(final_merge$name[row], " "))
    genus = split[1]
    species = split[2]
    #final_merge$genus[row] = genus
    final_merge$species[row] = species
  }
  
  
  write.table(
    final_merge,
    "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/final_merged_occurrances.occ",
    row.names = F,
    sep = "\t"
  )
  
  
  ## okay data processing is done. now it is time to move on
}

#####

## Step 3: pick by species and move forward

if (FIRSTSTEP <= 3) {
  print("STEP3")
  clean = read.csv(
    "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/final_merged_occurrances.occ",
    sep = "\t",
    numerals = "no.loss",
    stringsAsFactors = F
  )
  sort(unique(clean$name))
  clean$longitude = as.numeric(clean$longitude)
  clean$subspecies[clean$subspecies == ""] = "unknown"
  clean = clean[clean$latitude > 0,]
  clean = clean[clean$longitude < 0,]
  clean = clean[!(is.na(clean$longitude)),]
  
  cleanLoc = clean[, c("name", "longitude", "latitude", "subspecies")]
  
  ## Step 3a: run anomaly detection
  
  detected = detectSpatialOutliers(cleanLoc, 1e-6)
  purged = detected[[1]]
  anom = detected[[2]]
  kept = detected[[3]]
  
  write.table(
    kept,
    "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/final_merged_occurrances_CLEAN.occ",
    row.names = F,
    sep = "\t"
  )
  
  
  
  
  
  #raster::plot(bgLayer, col="grey",colNA="darkgrey",main=paste("Anomalies"))
  plot(cleanLoc$longitude,
       cleanLoc$latitude,
       col = "lightgrey",
       pch = 0)
  points(purged$longitude,
         purged$latitude,
         col = "red",
         pch = "+")
  
  
}

## Step 4: run the full pipeline

if (FIRSTSTEP <= 4) {
  print("STEP4")
  ## to do -- turn this into a fxn or for loop
  
  kept = read.csv(
    "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/locality information/final_merged_occurrances_CLEAN_USA.occ",
    sep = "\t",
    numerals = "no.loss",
    stringsAsFactors = F
  )
  
  ## GENERATE FUNCTIONS
  ## ALL OF THESE FUNCTIONS NEED TO LAYER UNKNOWN POINTS ON LAST
  {
    printPointsPng = function(species, subspecies, bg, loc_good) {
      print("PrintPointsPNG")
      png(paste("Subspecies_assignment_", species, "_USA.png", sep = ""))
      #print(length(subspecies))
      #col = floor(sqrt(length(subspecies)))
      #row = ceiling((sqrt(length(subspecies))))
      mf = n2mfrow(length(subspecies))
      if (mf[[1]] > mf[[2]]) {
        mf = rev(mf)
      }
      par(mfrow = mf)
      
      unkInd = which(levels(loc_good$subspecies) == "unknown")
      notUnkInd = which(levels(loc_good$subspecies) != "unknown")
      loc_good$subspecies = factor(loc_good$subspecies, levels(loc_good$subspecies)[c(unkInd, notUnkInd)])
      
      for (sub in subspecies) {
        print(sub)
        raster::plot(
          bg,
          col = "grey",
          colNA = "darkgrey",
          main = paste("assigned", sub, sep = " "),
          legend = F
        )
        temp = loc_good[loc_good$assigned == sub,]
        if (nrow(temp) > 0) {
          for (assignNum in 1:length(levels(temp$subspecies))) {
            #print(assignNum)
            assignSub = levels(temp$subspecies)[assignNum]
            mycol = palette()[assignNum]
            
            points(temp[temp$subspecies == assignSub, 2:3],
                   col = mycol,
                   pch = assignNum)
            
          }
          
          #points(temp[temp$assigned==sub,2:3],
          #       col=as.factor(temp$subspecies),
          #       pch=as.numeric(as.factor(temp$subspecies)))
          legend(
            "top",
            legend = as.factor(unique(temp$subspecies)),
            pch = unique(as.numeric(as.factor(
              temp$subspecies
            ))),
            bty = "n",
            col = as.factor(unique(temp$subspecies))
          )
        } else {
          print("NO POINTS")
        }
      }
      dev.off()
    }
    
    printPointsPdfGood = function(species, subspecies, bg, loc_good) {
      print("PrintPointsPdfGood")
      pdf(paste(
        "Subspecies_assignment_goodSubspp_",
        species,
        "_USA.pdf",
        sep = ""
      ))
      #print(length(subspecies))
      #col = floor(sqrt(length(subspecies)))
      #row = ceiling((sqrt(length(subspecies))))
      #mf = n2mfrow(length(subspecies))
      #if(mf[[1]] > mf[[2]]){mf = rev(mf)}
      #par(mfrow=mf)
      
      unkInd = which(levels(loc_good$subspecies) == "unknown")
      notUnkInd = which(levels(loc_good$subspecies) != "unknown")
      loc_good$subspecies = factor(loc_good$subspecies, levels(loc_good$subspecies)[c(unkInd, notUnkInd)])
      
      for (sub in subspecies) {
        print(sub)
        
        raster::plot(
          bg,
          col = "grey",
          colNA = "darkgrey",
          main = paste("assigned", sub, sep = " "),
          legend = F
        )
        temp = loc_good[loc_good$assigned == sub,]
        if (nrow(temp) > 0) {
          for (assignNum in 1:length(levels(temp$subspecies))) {
            #print(assignNum)
            assignSub = levels(temp$subspecies)[assignNum]
            mycol = palette()[assignNum]
            
            points(temp[temp$subspecies == assignSub, 2:3],
                   col = mycol,
                   pch = assignNum)
            
          }
          
          #points(temp[,2:3],
          #       col=as.factor(temp$subspecies),
          #       pch=as.numeric(as.factor(temp$subspecies)))
          legend(
            "top",
            legend = as.factor(unique(temp$subspecies)),
            pch = unique(as.numeric(as.factor(
              temp$subspecies
            ))),
            bty = "n",
            col = as.factor(unique(temp$subspecies))
          )
        } else {
          print("NO POINTS")
        }
      }
      dev.off()
    }
    
    printPointsPdfSuspect = function(species, subspecies, bg, loc_suspect) {
      print("PrintPointsSuspect")
      pdf(paste(
        "Subspecies_assignment_suspectSubspp_",
        species,
        "_USA.pdf",
        sep = ""
      ))
      #print(length(subspecies))
      #col = floor(sqrt(length(subspecies)))
      #row = ceiling((sqrt(length(subspecies))))
      #mf = n2mfrow(length(subspecies))
      #if(mf[[1]] > mf[[2]]){mf = rev(mf)}
      #par(mfrow=mf)
      
      unkInd = which(levels(loc_suspect$subspecies) == "unknown")
      notUnkInd = which(levels(loc_suspect$subspecies) != "unknown")
      loc_suspect$subspecies = factor(loc_suspect$subspecies,
                                      levels(loc_suspect$subspecies)[c(unkInd, notUnkInd)])
      
      for (sub in subspecies) {
        print(sub)
        raster::plot(
          bg,
          col = "grey",
          colNA = "darkgrey",
          main = paste("assigned", sub, sep = " "),
          legend = F
        )
        temp = loc_suspect[loc_suspect$assigned == sub,]
        if (nrow(temp) > 0) {
          for (assignNum in 1:length(levels(temp$subspecies))) {
            #print(assignNum)
            assignSub = levels(temp$subspecies)[assignNum]
            mycol = palette()[assignNum]
            
            points(temp[temp$subspecies == assignSub, 2:3],
                   col = mycol,
                   pch = assignNum)
            
          }
          
          
          #points(temp[temp$assigned==sub,2:3],
          #       col=as.factor(temp$subspecies),
          #       pch=as.numeric(as.factor(temp$subspecies)))
          legend(
            "top",
            legend = as.factor(unique(temp$subspecies)),
            pch = unique(as.numeric(as.factor(
              temp$subspecies
            ))),
            bty = "n",
            col = as.factor(unique(temp$subspecies))
          )
        } else {
          print("NO POINTS")
        }
      }
      dev.off()
    }
    
    outputProcessedSpecies = function(processedSpecies,
                                      outputDir,
                                      species,
                                      subspecies,
                                      bg) {
      palette(c(adjustcolor(palette()[1], alpha.f = 0.5), palette()[2:length(palette())]))
      
      loc_suspect = processedSpecies$loc_suspect
      loc_good = processedSpecies$loc_good
      pol = processedSpecies$pol
      
      print("Writing Tables")
      write.table(
        loc_suspect,
        paste(
          paste(
            "SuspectSubspp",
            species,
            paste(subspecies, collapse = " "),
            sep = "_"
          ),
          "_USA.txt",
          sep = ""
        ),
        quote = FALSE,
        sep = "\t",
        row.names = FALSE
      )
      write.table(
        loc_good,
        paste(
          paste(
            "GoodSubspp",
            species,
            paste(subspecies, collapse = " "),
            sep = "_"
          ),
          "_USA.txt",
          sep = ""
        ),
        quote = FALSE,
        sep = "\t",
        row.names = FALSE
      )
      
      print("Saving polygons")
      for (i in 1:length(pol)) {
        obj = pol[[i]]
        name = names(pol)[i]
        #print(name)
        obj2 = sp::SpatialPolygonsDataFrame(obj, data = as.data.frame(rep(1, length(obj))),
                                            match.ID = FALSE)
        #print(obj)
        #print("---")
        rgdal::writeOGR(
          obj = obj2,
          dsn = paste(paste("Polygon", species, name, sep = "_"), "_USA.shp", sep = ""),
          layer = name,
          driver = "ESRI Shapefile"
        )
        #print("end")
      }
      
      print("Printing PNG and PDF files")
      printPointsPng(
        species = species,
        subspecies = subspecies,
        bg = bg,
        loc_good = loc_good
      )
      printPointsPdfGood(
        species = species,
        subspecies = subspecies,
        bg = bg,
        loc_good = loc_good
      )
      printPointsPdfSuspect(
        species = species,
        subspecies = subspecies,
        bg = bg,
        loc_suspect = loc_suspect
      )
      
      palette("default")
      
    }
  }
  
  runOneSpecies = function(spp = "Auriparus flaviceps",
                           dataframe = kept,
                           outPath = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/OUTPUTS/",
                           bg = Env[[1]]) {
    species_dataframe = droplevels(unique(dataframe[dataframe$name == spp,]))
    species_dataframe
    
    subspecies = as.character(sort(unique(species_dataframe$subspecies)))
    subspecies = subspecies[!(subspecies %in% c("", "unknown"))]
    rownames(species_dataframe) = seq(1, nrow(species_dataframe))
    
    
    
    print("Plotting")
    png(paste("Labeled occurences_", spp, "_manual_USA.png", sep = ""))
    par(mfrow = c(2, 1), mar = c(2, 2, 1, 0))
    
    no_unk = species_dataframe[species_dataframe$subspecies != "unknown",]
    raster::plot(
      bg,
      col = "grey",
      colNA = "darkgrey",
      main = spp,
      legend = F
    )
    points(
      as.numeric(no_unk$longitude),
      as.numeric(no_unk$latitude),
      col = as.factor(no_unk$subspecies),
      pch = c(seq(0, length(
        unique(no_unk$subspecies)
      )))[as.factor(no_unk$subspecies)]
    )
    legend(
      "topright",
      legend = as.factor(unique(no_unk$subspecies)),
      pch = as.factor(unique(no_unk$subspecies)),
      bty = "n",
      col = as.factor(unique(no_unk$subspecies))
    )
    
    unk =  species_dataframe[species_dataframe$subspecies == "unknown",]
    raster::plot(
      bg,
      col = "grey",
      colNA = "darkgrey",
      main = "",
      legend = F
    )
    points(
      as.numeric(unk$longitude),
      as.numeric(unk$latitude),
      col = as.factor(unk$subspecies),
      pch = c(seq(0, length(
        unique(unk$subspecies)
      )))[as.factor(unk$subspecies)]
    )
    legend(
      "topright",
      legend = as.factor(unique(unk$subspecies)),
      pch = as.factor(unique(unk$subspecies)),
      bty = "n",
      col = as.factor(unique(unk$subspecies))
    )
    
    dev.off()
    
    ## anomaly detection only on this species
    
    ## don't need to do this twice if calling the full pipeline below
    
    detected_spp = detectSpatialOutliers(species_dataframe, 1e-6)
    purged_spp = detected_spp[[1]]
    anom_spp = detected_spp[[2]]
    kept_spp = detected_spp[[3]]
    
    #raster::plot(bgLayer, col="grey",colNA="darkgrey",main=paste("Anomalies"))
    # plot(species_dataframe$longitude,
    #      species_dataframe$latitude,
    #      col="lightgrey",pch=0)
    # points(purged_spp$longitude,purged_spp$latitude,
    #        col="red",pch="+")
    
    
    processedSpecies = subsppLabelR::databaseToAssignedSubspecies(
      spp = spp,
      subsppList = subspecies,
      pointLimit = 1,
      dbToQuery = c("gbif"),
      quantile = 0.95,
      xmin = -125,
      xmax = -60,
      ymin = 10,
      ymax = 50,
      plotIt = T,
      bgLayer = bg,
      outputDir = paste(outPath, "2.processed/images/", sep = ""),
      datafile = species_dataframe,
      epsilon = 1e-6
    )
    
    outputProcessedSpecies(processedSpecies, outPath, spp, subspecies, bg)
    
  }
  
  for (i in 1:length(unique(kept$name))) {
    #if (!(i %in% c(999999))){
    
    species = unique(kept$name)[i]
    
    focal_taxa = c("Amphispiza bilineata",
                   "Auriparus flaviceps",
                   "Campylorhynchus brunneicapillus",
                   "Cardinalis sinuatus",
                   "Melozone fusca",
                   "Phainopepla nitens",
                   "Polioptila melanura",
                   "Toxostoma crissale",
                   "Toxostoma curvirostre",
                   "Vireo bellii")
    
    # focal_taxa = c(
    #   "Melozone fusca",
    #   "Phainopepla nitens",
    #   "Polioptila melanura",
    #   "Toxostoma crissale",
    #   "Toxostoma curvirostre",
    #   "Vireo bellii"
    # )
    
    if ((species %in% focal_taxa)) {
      print(species)
      ## not working for matching subspp: callipepla california/douglassii
      
      runOneSpecies(spp = species,
                    dataframe = kept,
                    outPath = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/OUTPUTS/")
      print("xxxxxxxxxx")
    }
    
    
  }
  
  
  
  
  
  ## chooing auriparus flaviceps, has most localities
  
  #####
  
  
  
  ## test out k nearest neighbors instead
  {
    good = read.csv(
      "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/OUTPUTS/2.processed/data/GoodSubspp_asiatica_acaciarum flaviceps fraterculus hidalgensis lamprocephalus ornatus sinaloae.txt",
      header = T,
      sep = "\t"
    )
    #good = processedSpecies$loc_good
    susp = read.csv(
      "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/OUTPUTS/2.processed/data/SuspectSubspp_asiatica_acaciarum flaviceps fraterculus hidalgensis lamprocephalus ornatus sinaloae.txt",
      header = T,
      sep = "\t"
    )
    #susp = processedSpecies$loc_suspect
    alls = rbind(good, susp)
    
    knn = class::knn(
      train = good[, c(2:3)],
      test = susp[, c(2:3)],
      k = 8,
      cl = good$assigned
    )
    knn_vc = class::knn.cv(train = good[, c(2:3)],
                           k = 8,
                           cl = good$assigned)
    good$knnvc = knn_vc
    
    training = Auriparus_flaviceps[Auriparus_flaviceps$subspecies != "unknown",]
    testing = Auriparus_flaviceps[Auriparus_flaviceps$subspecies == "unknown",]
    
    knn_train = class::knn.cv(train = training[, c(2:3)],
                              k = 8,
                              cl = training$subspecies)
    knn_test = class::knn(
      train = training[, c(2:3)],
      test = testing[, c(2:3)],
      k = 8,
      cl = training$subspecies
    )
    testing$knn = knn_test
    training$knn = knn_train
    
    Auriparus_knn = rbind(testing, training)
    
    plot(Auriparus_knn$longitude,
         Auriparus_knn$latitude,
         col = as.factor(Auriparus_knn$knn))
    legend(
      "topright",
      legend = as.factor(unique(Auriparus_knn$knn)),
      bty = "n",
      pch = 1,
      col = as.factor(unique(Auriparus_knn$knn))
    )
  }
  #####
  
  ### actually run the niche models?
  
}

## Step 5: spThin

if (FIRSTSTEP <= 5) {
  print("STEP5")
  #filename = "AllSpecies_NotThinned_NoUncertains_WITHINPOLYGON.txt"
  
  focal_taxa = c(
    "Amphispiza_bilineata",
    "Auriparus_flaviceps",
    "Campylorhynchus_brunneicapillus",
    "Cardinalis_sinuatus",
    "Melozone_fusca",
    "Phainopepla_nitens",
    "Polioptila_melanura",
    "Toxostoma_crissale",
    "Toxostoma_curvirostre",
    "Vireo_bellii"
  )
  
  
  for (species in focal_taxa) {
    ## search for the filename contaning focal_taxa
    print(species)
    
    path = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/OUTPUTS/2.processed/data/CLEAN DATA/"
    setwd(path)
    
    spp_df = read.csv(paste(species, "_clean_USA.txt", sep = ""), sep = "\t")
    spp_df = unique(spp_df[, c(1:3)])
    
    ## run spthin on the glommed data
    
    thin(
      loc.data = spp_df,
      lat.col = "latitude",
      long.col = "longitude",
      spec.col = "name",
      thin.par = 10,
      reps = 1,
      locs.thinned.list.return = TRUE,
      write.files = TRUE,
      max.files = 3,
      out.dir = path,
      out.base = species,
      write.log.file = TRUE,
      log.file = paste(species, "_USA_thin.log", sep = "")
    )
    
  }
  
}

## Step 6: after spThin -- run basic ENMs 

if (FIRSTSTEP <= 6) {
  
  ## add any subsequent layers
  
  layers_to_add = c("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/soilgrids_data/SOIL-ACIDITY-ACDWRB_M_ss_250m_ll.tif",
    "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/CEC/CEC_Lakes_and_Rivers_Shapefile/NA_Lakes_and_Rivers/distance_to_river_resamp_clipped.asc",
                    "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/elevation_and_derivatives/roughness_NASA_clipped_resampled.tif",
                    "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/elevation_and_derivatives/elevation_NASA_clipped_resampled.tif",
                    "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/CEC/Land_Cover_2005v2_TIFF/LandCover_2005v2/data/NA_LandCover_2005V2_25haMMU_wgs84.asc",
                    "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/elevation_and_derivatives/NASA_relief_clip_resamp.tif",
                    "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/elevation_and_derivatives/NASA_slope_clip_resample.tif",
                    "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/Landfire_VegType/US_140EVT_20180618/us_140evt_converted_resampled.asc",
                    "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/soilgrids_data/SOIL-TEXTURE-0CM-TEXMHT_M_sl1_250m_ll.tif")
  
  names_layers = c("acidity","rivers","roughness","elevation","landcover","relief","slope","landfire","texture")
  
  # landfire = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/Landfire_VegType/US_140EVT_20180618/us_140evt.asc"
  # landfire_r = raster(landfire)
  # land = raster::crop(landfire_r, ext)
  # writeRaster(land,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/Landfire_VegType/US_140EVT_20180618/us_140evt_CROPPED.asc")
  
  for (i in 1:length(layers_to_add)) {
    print(i)
    lay = raster(layers_to_add[i])
    print(names(lay))
    print("crop")
    lay = raster::crop(lay, ext)
    if (i==1) {Env2 = lay}
    if(res(lay) != res(Env2)) {
    print("resample")
    lay = raster::resample(lay,Env2,method="ngb")
    }
    if (i==1) {ext = raster::extent(lay)}
    print("crop")
    lay = raster::crop(lay, ext)
    print(lay)
    #lay2 = resample(lay,bg,method="ngb")
    #writeRaster(lay2, filename="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/Landfire_VegType/US_140EVT_20180618/us_140evt_converted_resampled.asc", options="INTERLEAVE=BAND", overwrite=TRUE,
    #            format="ascii")
    if (i!= 1) {Env2 = stack(Env2,lay) } else {Env2 = lay}
    

  }
  
  writeRaster(Env2, filename="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/ENMS_multilayer2.tif",
              options="INTERLEAVE=BAND", overwrite=TRUE)
  
  if (res(Env2) != res(Env)) {
    print("resampling worldclim")
    Env = resample(Env,Env2,method="ngb")
  }
  
  Env3 = stack(Env,Env2)
  
  # names(Env3)=c("bio1","bio10","bio11","bio12","bio13","bio14","bio15",
  #         "bio16","bio17","bio18","bio19","bio2","bio3","bio4",
  #         "bio5","bio6","bio7","bio8","bio9",names_layers)
  
  writeRaster(Env3, filename="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/ENMS_multilayer3.tif",
              options="INTERLEAVE=BAND", overwrite=TRUE)
  
  ## first test on one species 
  
  landfire = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/Landfire_VegType/US_140EVT_20180618/us_140evt_converted.tif"
  lf = raster(landfire)
  lay = raster::crop(lf, extent(Env))
  
  writeRaster(lay,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/Landfire_VegType/US_140EVT_20180618/us_140evt_converted_cropped2worldclim.tif")
  
  if (res(lay) != res(Env)) {
  
  lf = resample(lay,Env,method="ngb") 
  
  writeRaster(lf,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/Landfire_VegType/US_140EVT_20180618/us_140evt_converted_resampled2worldclim.asc")
  }
  Env4 = stack(Env,lf)
  writeRaster(Env4,filename="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/ENMS_multilayer4.tif",
              options="INTERLEAVE=BAND", overwrite=TRUE)
  

  Env = stack("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/ENMS_multilayer4.tif")
  
  values(Env$ENMS_multilayer4.29)[values(Env$ENMS_multilayer4.29) == -9999] = NA
  
  writeRaster(Env,filename="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/ENMS_multilayer4.tif",
              options="INTERLEAVE=BAND", overwrite=TRUE)
  
  
}

if (FIRSTSTEP <= 7) {
  
  plotStuff=T
  
  Env = stack("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/ENMS_multilayer4.tif")
  
  
  thinDirectory = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/OUTPUTS/2.processed/data/"
  setwd(thinDirectory)
  
  thinned_csvs = list.files(
    path = thinDirectory,
    pattern = "thin1_USA.csv$",
    full.names = T
  )
  
  ## to do manually:
  #thinned_csvs = c("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/OUTPUTS/2.processed/data/Toxostoma_curvirostre_thin1_nominate.csv",
  #                 "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/OUTPUTS/2.processed/data/Toxostoma_palmeri_thin1.csv")
  #thinfile = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/OUTPUTS/2.processed/data/Toxostoma_crissale_thin1.csv"
  
  outdirectory="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/OUTPUTS/3.enms/USA_ONLY/"
  setwd(outdirectory)
  
  for (thinfile in thinned_csvs) {
    
    print(thinfile)
    thinned = read.csv(thinfile)
    sppname = thinned$name[1]
    loc = cbind(thinned$longitude,thinned$latitude)
    extr = extract(Env, loc) ## get env data at each location
    loc = loc[!is.na(extr[,1]),] ## remove any locations that have no ENV data
    
    ## add in more of this:
    ## https://cran.r-project.org/web/packages/ENMeval/vignettes/ENMeval-vignette.html
    
    ## split data into four equal blocks and different feature classes 
    ## throws errors with hinge 
    res = ENMevaluate(occ=loc, env = Env, method='block', 
                      parallel=T, numCores=4, fc=c("L", "LQ"), 
                      RMvalues=seq(0.5,4,0.5), rasterPreds=T,
                      updateProgress = T)
    ## this uses your occurences and trimmed env data. it will
    ## separate the data into 4 blocks, split sample into 4 equal quads
    ## running in parallel, with 4 cores. 
    ## the feature classes used are Linear, Linear+Quadratic, and Hunge.
    ## the regularization values will range from 1/2 to 4 by 1/2s.
    ## and rasterPreds=T means not running the predict function we ran before
    ## and thus can't have any AICc
    
    ## now find the best models 
    ## minimize omission rate and optimize AUC values
    ## remeber to change ORmin to OR10 depending on sample size
    ## for omission rate, if small sample size and high confidence
    ## do mean.ORmin. if don't know, do mean.OR10
    setsort = res@results[order(res@results[,'avg.test.or10pct']),] ## previously Mean.ORmin
    setsort2 = setsort[order(setsort[,'avg.test.AUC'], decreasing=TRUE),] ## previously Mean.AUC
    top = setsort2[1,]
    print(top)
    write.csv(setsort2,file=paste("AddedPredThin_ResultsTable_",sppname,"_addedLayers_zoomedin.csv",sep=""))
    
    best = which(as.character(res@results[,1]) == as.character(setsort2[1,1]))
    
    if(plotStuff==T) {
      pred.raw = predict(Env, res@models[[best]])
      writeRaster(pred.raw,filename=paste("AddedPredThin_BestModel_",sppname,"_addedLayers_zoomedin.asc",sep=""),
                  format="ascii",overwrite=T)
      png(paste("PredThin_BestModel_",sppname,"_addedLayers_zoomedin.png",sep=""))
      plot(pred.raw, col=viridis::viridis(99),main=sppname)
      points(loc,pch=3,col=rgb(1,0,0,0.1))
      dev.off()
      
      png(paste("PredThin_avg.test.AUC_",sppname,"_addedLayers_zoomedin.png",sep=""))
      eval.plot(res@results, "avg.test.AUC")
      dev.off()
    }
    
    ## with thresh model
    ev.set <- evaluate(thinned[,2:3], res@bg.pts, res@models[[best]], Env)
    th1 = threshold(ev.set) ## omission options
    
    write.csv(th1,file=paste("AddedPredThin_ThreshTable_",sppname,"_addedLayers_zoomedin.csv",sep=""))
    
    if(plotStuff==T) {
      
      p1.nomit = pred.raw>= th1$no_omission ## the highest thresh where no points are omittec
      p1.equal = pred.raw>= th1$equal_sens_spec ## equal sensitivity and specificity according to ROC curve
      p1.spsen = pred.raw>= th1$spec_sens ## equal sensitivity and specificity according to ROC curve
      
      writeRaster(p1.nomit,filename=paste("Thresh_NoOmission_BestModel_",sppname,"_addedLayers_zoomedin.asc",sep=""),
                  format="ascii",overwrite=T)
      writeRaster(p1.equal,filename=paste("Thresh_EqualSensSpec_BestModel_",sppname,"_addedLayers_zoomedin.asc",sep=""),
                  format="ascii",overwrite=T)
      writeRaster(p1.spsen,filename=paste("Thresh_MaxSpecSens_BestModel_",sppname,"_addedLayers_zoomedin.asc",sep=""),
                  format="ascii",overwrite=T)
      png(paste("Thresh_NoOmission_BestModel_",sppname,"_addedLayers_zoomedin.png",sep=""))
      plot(p1.nomit, col=viridis::viridis(99))
      points(loc,pch=3,col=rgb(1,0,0,0.1))
      dev.off()
      png(paste("Thresh_EqualSensSpec_BestModel_",sppname,"_addedLayers_zoomedin.png",sep=""))
      plot(p1.equal, col=viridis::viridis(99))
      points(loc,pch=3,col=rgb(1,0,0,0.1))
      dev.off()
      png(paste("Thresh_MaxSpecSens_BestModel_",sppname,"_addedLayers_zoomedin.png",sep=""))
      plot(p1.spsen, col=viridis::viridis(99))
      points(loc,pch=3,col=rgb(1,0,0,0.1))
      dev.off()
    }
    
    ## and now with mcp background data 
    mcp <- function (xy) { 
      xy <- as.data.frame(coordinates(xy))
      coords.t <- chull(xy[, 1], xy[, 2])
      xy.bord <- xy[coords.t, ]
      xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
      return(SpatialPolygons(list(Polygons(list(Polygon(as.matrix(xy.bord))), 1))))
    }
    
    MCP.locs = mcp(thinned[,2:3])
    
    if(plotStuff==T) {
      png(paste("MinConvexPolygon_",sppname,"_addedLayers_zoomedin.png",sep=""))
      plot(Env[[1]], col=viridis::viridis(99))
      plot(MCP.locs, add=T)
      points(loc,pch=3,col=rgb(1,0,0,0.1))
      dev.off()
    }
    
    ## background points inside polygon only, and also just on landmasses
    envPoly <- rasterToPolygons(Env[[1]], fun=NULL, na.rm=T)
    bg.area.locs <- gIntersection(envPoly, MCP.locs)
    MCP.raster.locs <- rasterize(bg.area.locs, Env[[1]])
    bg.points.locs <- randomPoints(mask = MCP.raster.locs, n = 10000)
    if (plotStuff==T) {
      png(paste("MinConvexPolygon_points_",sppname,"_addedLayers_zoomedin.png",sep=""))
      plot(Env[[1]], col=viridis::viridis(99))
      points(bg.points.locs)
      points(loc,pch=3,col=rgb(1,0,0,0.1))
      dev.off()
    }
    
    ## with only this polygon of bg points
    ## keeps failing -- non conformable arguments
    
    ## debugging
    bgpointsres = ENMevaluate(occ=loc, env = Env, method='block',  
                              bg.coords=bg.points.locs,
                              parallel=T, numCores=4, fc=c("L", "LQ"), 
                              RMvalues=seq(0.5,4,0.5), rasterPreds=T,
                              updateProgress = T)
    
    png(paste("BgPointsPlot_Mean.AUC_",sppname,"_addedLayers_zoomedin.png",sep=""))
    eval.plot(bgpointsres@results, "avg.test.AUC")
    dev.off()
    
    setsort = bgpointsres@results[order(bgpointsres@results[,'avg.test.or10pct']),]
    head(setsort)
    setsort2 = setsort[order(setsort[,'avg.test.AUC'], decreasing=TRUE),]
    head(setsort2)
    top = setsort2[1,]
    write.csv(setsort2,file=paste("BgPoints_ResultsTable_",sppname,"_addedLayers_zoomedin.csv",sep=""))
    
    
    best.bg = which(as.character(bgpointsres@results[,1]) == as.character(setsort2[1,1]))
    if (plotStuff==T) {
      pred.bg = predict(Env, bgpointsres@models[[best.bg]])
      writeRaster(pred.bg,filename=paste("BgPointsPlot_PredRaw_BestModel_",sppname,"_addedLayers_zoomedin.asc",sep=""),
                  format="ascii",overwrite=T)
      png(paste("BgPointsPlot_PredRaw_BestModel_",sppname,"_addedLayers_zoomedin.png",sep=""))
      plot(pred.bg, col=viridis::viridis(99))
      points(loc,pch=3,col=rgb(1,0,0,0.1))
      #points(bg.points.thinlocs,col="red")
      dev.off()
    }
    print(bgpointsres@results[best.bg,])
    
    ev.set <- evaluate(thinned[,2:3], bgpointsres@bg.pts, 
                       bgpointsres@models[[best.bg]], Env)
    th2 = threshold(ev.set) ## omission options
    
    write.csv(th2,file=paste("BgPoints_ThreshTable_",sppname,"_addedLayers_zoomedin.csv",sep=""))
    
    if(plotStuff==T){
      p2.nomit = pred.bg >= th2$no_omission ## the highest thresh where no points are omittec
      p2.equal = pred.bg >= th2$equal_sens_spec ## equal sensitivity and specificity according to ROC curve
      p2.spsen = pred.bg >= th2$spec_sens 
      
      writeRaster(p2.nomit,filename=paste("BgPoints_Thresh_NoOmission_BestModel_",sppname,"_addedLayers_zoomedin.asc",sep=""),
                  format="ascii",overwrite=T)
      writeRaster(p2.equal,filename=paste("BgPoints_Thresh_EqualSensSpec_BestModel_",sppname,"_addedLayers_zoomedin.asc",sep=""),
                  format="ascii",overwrite=T)
      writeRaster(p2.spsen,filename=paste("BgPoints_Thresh_MaxSpecSens_BestModel_",sppname,"_addedLayers_zoomedin.asc",sep=""),
                  format="ascii",overwrite=T)
      png(paste("BgPoints_Thresh_NoOmission_BestModel_",sppname,"_addedLayers_zoomedin.png",sep=""))
      plot(p2.nomit, col=viridis::viridis(99))
      points(loc,pch=3,col=rgb(1,0,0,0.1))
      dev.off()
      png(paste("BgPoints_Thresh_EqualSensSpec_BestModel_",sppname,"_addedLayers_zoomedin.png",sep=""))
      plot(p2.equal, col=viridis::viridis(99))
      points(loc,pch=3,col=rgb(1,0,0,0.1))
      dev.off()
      png(paste("BgPoints_Thresh_MaxSpecSens_BestModel_",sppname,"_addedLayers_zoomedin.png",sep=""))
      plot(p2.spsen, col=viridis::viridis(99))
      points(loc,pch=3,col=rgb(1,0,0,0.1))
      dev.off()
    }
    
    ## add in past projections
    
  }
  
}



