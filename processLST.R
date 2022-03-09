####Process LST rasters and create layers needed for RF downstream analysis####
############Mask the LST raster by the QA raster, apply scale factor, convert to degC, calculate median and sd, write outputs

###REPROCESS LST RASTERS TO FIX PROJECTION ISSUE####
library(parallel)
library(snow)
library(raster)
library(rgdal)
library(stringr)

#set temp dir to slate folder
#install.packages('unixtools',repos='http://www.rforge.net/')
unixtools::set.tempdir("/N/slate/farellam")
tempdir()


processLST <- function(year, tile){
  print(paste("now processing", tile, year, sep=" "))
  #get a list of all LST files
  listST <- list.files(path="/N/project/CoverCrops/CoverCrops/BulkOrderIndianaProvisionalSurfaceTemp/U.S. Landsat 4-8 ARD/", recursive = TRUE, pattern = "*_ST.tif", full.names = TRUE)
  #subset tiles for a single tile/year
  onetile <- listST[grep(tile, listST)]
  nxtyr <- year+1
  fall_oneyr <- c(onetile[grep(paste(year,"11", sep=""), onetile)], onetile[grep(paste(year, "12", sep=""), onetile)], 
                  onetile[grep(paste(nxtyr,"01", sep=""), onetile)], onetile[grep(paste(nxtyr, "02", sep=""), onetile)], onetile[grep(paste(nxtyr, "03", sep=""), onetile)])
  #since files are processed in 2018 and 2019 need to further subset to only include filed acquired during desired date range
  goodfiles <- c(paste(year, "11", sep=""), paste(year, "12", sep=""), paste(nxtyr, "01", sep=""), paste(nxtyr, "02", sep=""), paste(nxtyr, "03", sep=""))
  fall_oneyr <- unique(subset(fall_oneyr, str_sub(fall_oneyr, -32, -27) %in% goodfiles))
  
  #get a list of the QA files
  listQA <- list.files(path="/N/project/CoverCrops/CoverCrops/BulkOrderIndianaProvisionalSurfaceTemp/U.S. Landsat 4-8 ARD/", recursive=TRUE, pattern="*STQA.tif", full.names=TRUE)
  cloud <- listQA[grep(tile, listQA)]
  msk_oneyr <- c(cloud[grep(paste(year,"11", sep=""), cloud)], cloud[grep(paste(year, "12", sep=""), cloud)], 
                 cloud[grep(paste(nxtyr,"01", sep=""), cloud)], cloud[grep(paste(nxtyr, "02", sep=""), cloud)], cloud[grep(paste(nxtyr, "03", sep=""), cloud)])
  msk_oneyr <- unique(subset(msk_oneyr, str_sub(msk_oneyr, -34, -29) %in% goodfiles))
  
  
  
  #make raster stack of the LST rasters for each season x year
  lst.stack <- raster::stack(fall_oneyr)
  #make raster stack of the LST QA rasters for each season x year
  msk.stack <- raster::stack(msk_oneyr)
  
  beginCluster(numcore)
  print("recoding QA raster")
  msk.stack <- msk.stack*0.01 #scale factor: https://prd-wret.s3.us-west-2.amazonaws.com/assets/palladium/production/atoms/files/LSDS-1330-LandsatSurfaceTemperature_ProductGuide-v2.pdf
  msk.stack[msk.stack > 7] <- NA
  gc()
  #mask, reproject, convert to degC, and write raster
  print("masking fall LST")
  fall_mask <- raster::mask(lst.stack, msk.stack)
  
  print("converting fall LST to celsius and applying the scale factor")
  lst_C <- ((fall_mask/10)-273.15)
  
  print("calculating median and sd values")
  med_fall <- clusterR(lst_C, calc, args=list(median, na.rm=TRUE))
  sd_fall <- clusterR(lst_C, calc, args=list(sd, na.rm=TRUE))
  
  print("reprojecting median and sd rasters")
  crs.final <- crs("+proj=utm +zone=16 +ellps=WGS84 +units=m +no_defs")
  med_fall_p <- projectRaster(med_fall, crs=crs.final)
  sd_fall_p <- projectRaster(sd_fall, crs=crs.final)
  
  
  endCluster()
  
  print("writing output rasters")
  writeRaster(med_fall_p, paste("/N/project/CoverCrops/CoverCrops/HLS30_Indiana/utm_proj/LSTmedFall", year, "_", tile, ".tif", sep=""), overwrite=TRUE)
  writeRaster(sd_fall_p, paste("/N/project/CoverCrops/CoverCrops/HLS30_Indiana/utm_proj/LSTsdFall", year, "_", tile, ".tif", sep=""), overwrite=TRUE)
  
}

rasterOptions(tmpdir="/N/slate/farellam/", tmptime = 24, progress="", timer=TRUE, overwrite = T, chunksize=8e10, maxmemory=5e8)
numcore = 16
processLST(tile="023010", year=2014)
processLST(tile="023010", year=2015)
processLST(tile="023010", year=2016)
processLST(tile="023010", year=2017)
processLST(tile="023010", year=2018)
processLST(tile="023010", year=2019)

