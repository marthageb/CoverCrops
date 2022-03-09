#####Create cloud mask from QA layer, apply mask and make raster stack for entire year (Nov-Oct)####
###HLS data downloaded from hls.gsfc.nasa.gov/data/v1.4/L30 for tiles: 16TDL, 16TDK, 16SDJ, 16SDH, 16TEL, 16TEK, 16SEJ, 16SEH, 16TFL, 16TFK, 16SFJ, 16SFH for dates: November 2014 - October 2020
###Before running this code .hdf4 files need to be converted to .tif files with 'process_hdf4.R' 

library(stringr)
library(raster)

#set temp dir to slate folder
#install.packages('unixtools',repos='http://www.rforge.net/')
unixtools::set.tempdir("/N/slate/farellam")
tempdir()

#function to stack each band by the whole year (Nov - Oct), mask by the QA layer, and apply the scale factor
stack_year <- function(tile, year){
  all_beg <- Sys.time()
  yr_wd <- paste("/N/project/CoverCrops/CoverCrops/HLS_L30_v1.4/Tiffs/", tile, "_", year, "_Tiffs", sep="")
  setwd(yr_wd)
  yr_files <- list.files(yr_wd, pattern="tif$", full.names = TRUE)
  
  #subset files to only include Nov - Dec
  doy <- as.numeric(str_sub(yr_files, -12, -10))
  novdec <- doy[doy > 305]
  total <- length(novdec)/11
  novdec <- novdec[1:total]
  novdec <- paste("00", novdec, sep="")
  novdec <- str_sub(novdec, -3, -1)
  novdec <- rep(novdec, each=11) #repeat each doy 11 times (because there's 11 bands)
  
  beg_name <- unique(str_sub(yr_files, 1, -13))
  beg_name <- rep(beg_name, total)
  end_name <- unique(str_sub(yr_files, -9,-1))
  end_name <- rep(end_name, length(beg_name))
  
  cur_yr_files <- paste(beg_name, novdec, end_name, sep="")
  
  #get the files for the next year
  nxtyr <- year + 1
  yr_wd <- paste("/N/project/CoverCrops/CoverCrops/HLS_L30_v1.4/Tiffs/", tile, "_", nxtyr, "_Tiffs", sep="")
  setwd(yr_wd)
  yr_files <- list.files(yr_wd, pattern="tif$", full.names = TRUE)
  #subset files to only include Jan - Oct
  doy <- as.numeric(str_sub(yr_files, -12, -10))
  janoct <- doy[doy < 304]
  total <- length(janoct)/11
  janoct <- janoct[1:total]
  janoct <- paste("00", janoct, sep="")
  janoct <- str_sub(janoct, -3, -1)
  janoct <- rep(janoct, each=11) #repeat each doy 11 times (because there's 11 bands)
  
  beg_name <- unique(str_sub(yr_files, 1, -13))
  beg_name <- rep(beg_name, total)
  end_name <- unique(str_sub(yr_files, -9,-1))
  end_name <- rep(end_name, length(beg_name))
  
  nxt_yr_files <- paste(beg_name, janoct, end_name, sep="")
  
  #get all files names for each bad for the entire year 
  b2_files <- c(cur_yr_files[grep("Band2", cur_yr_files)],  nxt_yr_files[grep("Band2", nxt_yr_files)])
  b3_files <- c(cur_yr_files[grep("Band3", cur_yr_files)],  nxt_yr_files[grep("Band3", nxt_yr_files)])
  b4_files <- c(cur_yr_files[grep("Band4", cur_yr_files)],  nxt_yr_files[grep("Band4", nxt_yr_files)])
  b5_files <- c(cur_yr_files[grep("Band5", cur_yr_files)],  nxt_yr_files[grep("Band5", nxt_yr_files)])
  b6_files <- c(cur_yr_files[grep("Band6", cur_yr_files)],  nxt_yr_files[grep("Band6", nxt_yr_files)])
  b7_files <- c(cur_yr_files[grep("Band7", cur_yr_files)],  nxt_yr_files[grep("Band7", nxt_yr_files)])
  b10_files <- c(cur_yr_files[grep("Band10", cur_yr_files)],  nxt_yr_files[grep("Band10", nxt_yr_files)])
  b11_files <- c(cur_yr_files[grep("Band11", cur_yr_files)],  nxt_yr_files[grep("Band11", nxt_yr_files)])
  qa_files <- c(cur_yr_files[grep("QA", cur_yr_files)],  nxt_yr_files[grep("QA", nxt_yr_files)])
  
  #Create mask....
  beginCluster(clustnum)
  start <- Sys.time()
  
  qa <- raster::stack(qa_files)
  
  #first determine if there's clouds or not (see 'QA bit math.xls' for calculation of bit values)
  clouds_fxn <- function(x){(floor(x/2)) - (floor((x/2)/2) * 2)}
  clouds <- clusterR(qa, clouds_fxn)
  print('made cloud mask')
  
  gc()
  #then determine if there's cloud shadow
  #shadow_fxn <- function(y){(floor(y/8)) - (floor((y/8)/2) * 2)}
  #shadow <- clusterR(qa, shadow_fxn)
  #print('made cloud shadow mask')
  
  #gc()
  #then determine if there's snow/ice
  #snow_fxn <- function(z){(floor(z/16)) - (floor((z/16)/2) * 2)}
  #snow <- clusterR(qa, snow_fxn)
  #print('made snow mask')
  
  #add the clouds, cloud shadow, and snow/ice layers together
  #qa <- clusterR(clouds, fun=function(x, y, z) {return(x + y)},
                # args = list(y = shadow))
 
  
  #any cells with values > 0 (aka they have clouds, cloud shadows, and/or snow/ice) assign NA
  #qa <- clusterR(qa, fun=function(x) {x[x>0] = NA; return(x)})
   qa <- clusterR(clouds, fun=function(x) {x[x>0] = NA; return(x)})
   endCluster()
  outname <- paste("qa", tile, year, sep="_")
  
  #write the QA stack
  writeRaster(qa, paste("/N/project/CoverCrops/CoverCrops/HLS_L30_v1.4/stacks/cloud_mask/", outname, ".tif", sep=""))
  end <- Sys.time()
  print(paste(end-start, "time difference", sep=" "))
  print("QA stack written, now processing B3")
  
  gc()
  b_stack <- raster::stack(b2_files)
  names(b_stack) <- paste(str_sub(b2_files, -16, -13), str_sub(b2_files, -12, -10), sep="_")
  #apply the scale factor-- see Table 9 on the user guide: https://hls.gsfc.nasa.gov/wp-content/uploads/2019/01/HLS.v1.4.UserGuide_draft_ver3.1.pdf
  scale_fxn <- function(y){
    y * 0.0001
  }
  #apply the scale factor
  #beginCluster(clustnum)
  #b_stack <- clusterR(b_stack, scale_fxn)
  #endCluster()
  #apply the QA mask
  b_mask <-raster::mask(b_stack, qa)
  outname <- paste("b2", tile, year, sep="_")
  writeRaster(b_mask, paste("/N/project/CoverCrops/CoverCrops/HLS_L30_v1.4/stacks/cloud_mask/", outname, ".tif", sep=""))
  print("B2 stack written, now processing B3")
  
  gc()
  b_stack <- raster::stack(b3_files)
  names(b_stack) <- paste(str_sub(b3_files, -16, -13), str_sub(b3_files, -12, -10), sep="_")
  #apply the scale factor
  #beginCluster(clustnum)
  #b_stack <- clusterR(b_stack, calc, args=list(fun=scale_fxn))
  #endCluster()
  #apply the QA mask
  b_mask <-raster::mask(b_stack, qa)
  outname <- paste("b3", tile, year, sep="_")
  writeRaster(b_mask, paste("/N/project/CoverCrops/CoverCrops/HLS_L30_v1.4/stacks/cloud_mask/", outname, ".tif", sep=""))
  print("B3 stack written, now processing B4")
  
  gc()
  b_stack <- raster::stack(b4_files)
  names(b_stack) <- paste(str_sub(b4_files, -16, -13), str_sub(b4_files, -12, -10), sep="_")
  #apply the scale factor
  #beginCluster(clustnum)
  #b_stack <- clusterR(b_stack, scale_fxn)
  #endCluster()
  #apply the QA mask
  b_mask <-raster::mask(b_stack, qa)
  outname <- paste("b4", tile, year, sep="_")
  writeRaster(b_mask, paste("/N/project/CoverCrops/CoverCrops/HLS_L30_v1.4/stacks/cloud_mask/", outname, ".tif", sep=""))
  print("B4 stack written, now processing B5")
  
  gc()
  b_stack <- raster::stack(b5_files)
  names(b_stack) <- paste(str_sub(b5_files, -16, -13), str_sub(b5_files, -12, -10), sep="_")
  #apply the scale factor
  #beginCluster(clustnum)
  #b_stack <- clusterR(b_stack, scale_fxn)
  #endCluster()
  #apply the QA mask
  b_mask <-raster::mask(b_stack, qa)
  outname <- paste("b5", tile, year, sep="_")
  writeRaster(b_mask, paste("/N/project/CoverCrops/CoverCrops/HLS_L30_v1.4/stacks/cloud_mask/", outname, ".tif", sep=""))
  print("B5 stack written, now processing B6")
  
  b_stack <- raster::stack(b6_files)
  names(b_stack) <- paste(str_sub(b6_files, -16, -13), str_sub(b6_files, -12, -10), sep="_")
  #apply the scale factor
  #beginCluster(clustnum)
  #b_stack <- clusterR(b_stack, scale_fxn)
  #endCluster()
  #apply the QA mask
  b_mask <-raster::mask(b_stack, qa)
  outname <- paste("b6", tile, year, sep="_")
  writeRaster(b_mask, paste("/N/project/CoverCrops/CoverCrops/HLS_L30_v1.4/stacks/cloud_mask/", outname, ".tif", sep=""))
  print("B6 stack written, now processing B7")
  
  gc()
  b_stack <- raster::stack(b7_files)
  names(b_stack) <- paste(str_sub(b7_files, -16, -13), str_sub(b7_files, -12, -10), sep="_")
  #apply the scale factor
  #beginCluster(clustnum)
  #b_stack <- clusterR(b_stack, scale_fxn)
  #endCluster()
  #apply the QA mask
  b_mask <-raster::mask(b_stack, qa)
  outname <- paste("b7", tile, year, sep="_")
  writeRaster(b_mask, paste("/N/project/CoverCrops/CoverCrops/HLS_L30_v1.4/stacks/cloud_mask/", outname, ".tif", sep=""))
  print("B7 stack written, now processing B10")
  
  gc()
  b_stack <- raster::stack(b10_files)
  names(b_stack) <- paste(str_sub(b10_files, -16, -13), str_sub(b10_files, -12, -10), sep="_")
  #apply the scale factor-- ***bands 10 and 11 have a different scale factor than the other bands***
  scale_fxn <- function(x){
    x * 0.01
  }
  #beginCluster(clustnum)
  #b_stack <- clusterR(b_stack, scale_fxn)
  #endCluster()
  #apply the QA mask
  b_mask <-raster::mask(b_stack, qa)
  outname <- paste("b10", tile, year, sep="_")
  writeRaster(b_mask, paste("/N/project/CoverCrops/CoverCrops/HLS_L30_v1.4/stacks/cloud_mask/", outname, ".tif", sep=""))
  print("B10 stack written, now processing B11")
  
  gc()
  b_stack <- raster::stack(b11_files)
  names(b_stack) <- paste(str_sub(b11_files, -16, -13), str_sub(b11_files, -12, -10), sep="_")
  #apply the scale factor
  #beginCluster(clustnum)
  #b_stack <- clusterR(b_stack, scale_fxn)
  #endCluster()
  #apply the QA mask
  b_mask <-raster::mask(b_stack, qa)
  outname <- paste("b11", tile, year, sep="_")
  writeRaster(b_mask, paste("/N/project/CoverCrops/CoverCrops/HLS_L30_v1.4/stacks/cloud_mask/", outname, ".tif", sep=""))
  print("B11 stack written!")
  
  all_end <- Sys.time()
  print(paste('total time = ', all_end - all_beg, sep=""))
}


#increasing chunk size can speed things up when dealing w large rasters
rasterOptions(tmpdir="/N/slate/farellam/", tmptime = 24, progress="", timer=FALSE, overwrite = T, chunksize=8e10, maxmemory=5e8)

clustnum = 12


stack_year(tile = "16TDL", year = 2017)
stack_year(tile = "16TDL", year = 2018)
