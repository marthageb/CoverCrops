#### Create stacks of the HLS data needed for RF analysis####
######## 1. load the stacks for each band get/apply band names
######## 2. calculate NDVI for the entire year (Nov- Oct). convert any NDVI values >1 or <0 to NA to mask out water/low veg areas.
######## 3. subset fullyear band stacks to off season (Nov- Mar)
######## 4. calculate VIs for the off season (NDTI, STI, Rn, Gn, Bn, ExG, ExGR)
######## 5. calculate summary stats (min, max, median) for the desired bands and VIs
######## 6. calculate NDVI amplitude and ratio
######## 7. stack all RF bands and write output. 
################ RF bands: b3_med, b3_min, b3_max, b5_med, b5_min, b5_max, b6_med, b6_min, b6_max, b10_med, b10_min, b10_max, 
################ b11_med, b11_min, b11_max, ndvi_med, ndvi_min, ndvi_max, ndti_med, ndti_min, ndti_max, sti_med, sti_min, sti_max, 
################ ExG_med, ExG_min, ExG_max, ExGR_med, ExGR_min, ExGR_max, full_ndvi_max, ndvi_amp, ndvi_ratio

library(raster)
library(doParallel)
library(stringr)

#set temp dir to slate folder
#install.packages('unixtools',repos='http://www.rforge.net/')
unixtools::set.tempdir("/N/slate/farellam")
tempdir()

####HLS RFdata###
#get dates from file names
names_fxn <- function(tile, year){
  yr_wd <- paste("/N/project/CoverCrops/CoverCrops/HLS_L30_v1.4/Tiffs/", tile, "_", year, "_Tiffs", sep="")
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
  
  #get all files names for each band for the entire year 
  file_names <- c(cur_yr_files[grep("Band3", cur_yr_files)],  nxt_yr_files[grep("Band3", nxt_yr_files)])
  return(file_names)
}




make_stacks <- function(tile, year) {
  gc()
  print(paste("now processing", tile, year, sep=" "))
  
  setwd("/N/project/CoverCrops/CoverCrops/HLS_L30_v1.4/stacks/cloud_mask/")
  all <- list.files("/N/project/CoverCrops/CoverCrops/HLS_L30_v1.4/stacks/cloud_mask/")
  tile_files <- all[grep(tile, all)]
  tile_files <- tile_files[grep(year, tile_files)]
  
  b2 <- stack(tile_files[grep("b2", tile_files)])
  b3 <- stack(tile_files[grep("b3", tile_files)])
  b4 <- stack(tile_files[grep("b4", tile_files)])
  b5 <- stack(tile_files[grep("b5", tile_files)])
  b6 <- stack(tile_files[grep("b6", tile_files)])
  b7 <- stack(tile_files[grep("b7", tile_files)])
  b10 <- stack(tile_files[grep("b10", tile_files)])
  b11 <- stack(tile_files[grep("b11", tile_files)])
  
  names(b2) <- paste(str_sub(names_fxn(tile=tile, year=year), -16, -13), str_sub(names_fxn(tile=tile, year=year), -12, -10), sep="_")
  names(b3) <- paste(str_sub(names_fxn(tile=tile, year=year), -16, -13), str_sub(names_fxn(tile=tile, year=year), -12, -10), sep="_")
  names(b4) <- paste(str_sub(names_fxn(tile=tile, year=year), -16, -13), str_sub(names_fxn(tile=tile, year=year), -12, -10), sep="_")
  names(b5) <- paste(str_sub(names_fxn(tile=tile, year=year), -16, -13), str_sub(names_fxn(tile=tile, year=year), -12, -10), sep="_")
  names(b6) <- paste(str_sub(names_fxn(tile=tile, year=year), -16, -13), str_sub(names_fxn(tile=tile, year=year), -12, -10), sep="_")
  names(b7) <- paste(str_sub(names_fxn(tile=tile, year=year), -16, -13), str_sub(names_fxn(tile=tile, year=year), -12, -10), sep="_")
  names(b10) <- paste(str_sub(names_fxn(tile=tile, year=year), -16, -13), str_sub(names_fxn(tile=tile, year=year), -12, -10), sep="_")
  names(b11) <- paste(str_sub(names_fxn(tile=tile, year=year), -16, -13), str_sub(names_fxn(tile=tile, year=year), -12, -10), sep="_")
  
  #get raster band names for the off season (Nov - March)
  band_names <- names(b3)
  curr_year <- band_names[grep(year, band_names)]
  days <- as.numeric(str_sub(curr_year, -3, -1))
  nov <- days[days < 334]
  nov <- paste("00", nov, sep="")
  nov <- str_sub(nov, -3, -1)
  nov_bands <- paste("X", year, "_", nov, sep="")
  
  dec <- days[days > 334]
  dec <- paste("00", dec, sep="")
  dec <- str_sub(dec, -3, -1)
  dec_bands <- paste("X", year, "_", dec, sep="")
  
  nxt_year <- band_names[grep(year+1, band_names)]
  days <- as.numeric(str_sub(nxt_year, -3, -1))
  jan <- days[days < 31]
  jan <- paste("00", jan, sep="")
  jan <- str_sub(jan, -3, -1)
  jan_bands <- paste("X", year+1, "_", jan, sep="")
  
  feb <- days[days < 59]
  feb <- feb[feb > 31]
  feb <- paste("00", feb, sep="")
  feb <- str_sub(feb, -3, -1)
  feb_bands <- paste("X", year+1, "_", feb, sep="")
  
  mar <- days[days < 90]
  mar <- mar[mar > 59]
  mar <- paste("00", mar, sep="")
  mar <- str_sub(mar, -3, -1)
  mar_bands <- paste("X", year+1, "_", mar, sep="")
  
  wint_bands <- c(dec_bands, jan_bands)
  
  days <- days[days< 90]
  days <- paste("00", days, sep="")
  days <- str_sub(days, -3, -1)
  nxt_year <- paste("X", year+1, "_", days, sep="")
  
  off_seas <- c(curr_year, nxt_year)
  
  beginCluster(clustnum)
  
  #calculate ndvi for the whole year
  ndvi_full <- clusterR(b5, fun=function(x, y) {return((x-y)/(x+y))},
                        args = list(y=b4))
  
  
  #set NDVI < 0 to NA to remove water pixels and > 1 to NA
  
  ndvi_full <- clusterR(ndvi_full, fun=function(x){x[x>1] = NA; return(x)})
  ndvi_full <- clusterR(ndvi_full, fun=function(x){x[x<0] = NA; return(x)})
  
  print("finished with NDVI for the full year")
  gc()
  names(ndvi_full) <- paste(str_sub(names_fxn(tile=tile, year=year), -16, -13), str_sub(names_fxn(tile=tile, year=year), -12, -10), sep="_")
  
  #subset stacks to just the off season
  fall_ndvi <- raster::subset(ndvi_full, off_seas)
  fall_b2 <- raster::subset(b2, off_seas)
  fall_b3 <- raster::subset(b3, off_seas)
  fall_b4 <- raster::subset(b4, off_seas)
  fall_b5 <- raster::subset(b5, off_seas)
  fall_b6 <- raster::subset(b6, off_seas)
  fall_b7 <- raster::subset(b7, off_seas)
  fall_b10 <- raster::subset(b10, off_seas)
  fall_b11 <- raster::subset(b11, off_seas)
  
  
  #calculate ndti and sti
  fall_ndti <- clusterR(fall_b6, fun=function(x,y) {return((x-y)/(x+y))},
                        args = list(y = fall_b7))
  fall_sti <- clusterR(fall_b6, fun=function(x,y) {return(x/y)},
                       args = list(y = fall_b7))
  names(fall_ndti) <- names(fall_ndvi)
  names(fall_sti) <- names(fall_ndvi)
  
  
  print("ndti and sti for the off season calculated")
  gc()
  
  #calculate normalized R, G, B
  fall_Rn <- clusterR(fall_b4, fun=function(x, y, z) {return(x/(x+y+z))},
                      args = list(y=fall_b3, z=fall_b2))
  fall_Gn <- clusterR(fall_b3, fun=function(x, y, z) {return(x/(x+y+z))},
                  args = list(y=fall_b4, z=fall_b2))
  fall_Bn <- clusterR(fall_b2, fun=function(x, y, z) {return(x/(x+y+z))},
                      args = list(y=fall_b3, z=fall_b4))
  
  names(fall_Rn) <- names(fall_ndvi)
  names(fall_Gn) <- names(fall_ndvi)
  names(fall_Bn) <- names(fall_ndvi)
  print("normalized R, G, B for the off season calculated")
  gc()
  
  #calculate Excess Green and Excess Green - Red
  fall_ExG <- clusterR(fall_Gn, fun=function(x, y, z) {return((2*x) - y - z)},
                       args = list(y=fall_Rn, z=fall_Bn))
  fall_ExGR <- clusterR(fall_ExG, fun=function(x, y, z) {return(x - (1.4*y) - z)},
                        args = list(y=fall_Rn, z=fall_Gn))
  
  names(fall_ExG) <- names(fall_ndvi)
  names(fall_ExGR) <- names(fall_ndvi)
  print("excess Green and excess Green-Red for the off season calculated")
  gc()
  endCluster()
  
  
  
  #calcualte median, min, and max vals
  med_fun <- function(x){median(x, na.rm=TRUE)}
  
  beginCluster(clustnum)
  b3_med <- clusterR(fall_b3, calc, args=list(fun=med_fun))
  b3_max <- clusterR(fall_b3, calc, args=list(max, na.rm=TRUE))
  b3_min <- clusterR(fall_b3, calc, args=list(min, na.rm=TRUE))
  print("done with B3")
  
  gc()
  b5_med <- clusterR(fall_b5, calc, args=list(fun=med_fun))
  b5_max <- clusterR(fall_b5, calc, args=list(max, na.rm=TRUE))
  b5_min <- clusterR(fall_b5, calc, args=list(min, na.rm=TRUE))
  print("done with B5")
  
  gc()
  b6_med <- clusterR(fall_b6, calc, args=list(fun=med_fun))
  b6_max <- clusterR(fall_b6, calc, args=list(max, na.rm=TRUE))
  b6_min <- clusterR(fall_b6, calc, args=list(min, na.rm=TRUE))
  print("done with B6")
  
  gc()
  b10_med <- clusterR(fall_b10, calc, args=list(fun=med_fun))
  b10_max <- clusterR(fall_b10, calc, args=list(max, na.rm=TRUE))
  b10_min <- clusterR(fall_b10, calc, args=list(min, na.rm=TRUE))
  print("done with B10")
  
  gc()
  b11_med <- clusterR(fall_b11, calc, args=list(fun=med_fun))
  b11_max <- clusterR(fall_b11, calc, args=list(max, na.rm=TRUE))
  b11_min <- clusterR(fall_b11, calc, args=list(min, na.rm=TRUE))
  print("done with B11")
  
  gc()
  ndvi_med <- clusterR(fall_ndvi, calc, args=list(fun=med_fun))
  ndvi_max <- clusterR(fall_ndvi, calc, args=list(max, na.rm=TRUE))
  ndvi_min <- clusterR(fall_ndvi, calc, args=list(min, na.rm=TRUE))
  print("done with NDVI")
  
  gc()
  ndti_med <- clusterR(fall_ndti, calc, args=list(fun=med_fun))
  ndti_max <- clusterR(fall_ndti, calc, args=list(max, na.rm=TRUE))
  ndti_min <- clusterR(fall_ndti, calc, args=list(min, na.rm=TRUE))
  print("done with NDTI")
  
  gc()
  sti_med <- clusterR(fall_sti, calc, args=list(fun=med_fun))
  sti_max <- clusterR(fall_sti, calc, args=list(max, na.rm=TRUE))
  sti_min <- clusterR(fall_sti, calc, args=list(min, na.rm=TRUE))
  print("done with STI")
  
  gc()
  ExG_med <- clusterR(fall_ExG, calc, args=list(fun=med_fun))
  ExG_max <- clusterR(fall_ExG, calc, args=list(max, na.rm=TRUE))
  ExG_min <- clusterR(fall_ExG, calc, args=list(min, na.rm=TRUE))
  print("done with ExG")
  
  gc()
  ExGR_med <- clusterR(fall_ExGR, calc, args=list(fun=med_fun))
  ExGR_max <- clusterR(fall_ExGR, calc, args=list(max, na.rm=TRUE))
  ExGR_min <- clusterR(fall_ExGR, calc, args=list(min, na.rm=TRUE))
  print("done with ExG-R")
  
  #ndvi_med <- clusterR(fall_ndvi, calc, args=list(fun=med_fun))
  #ndti_med <- clusterR(fall_ndti, calc, args=list(fun=med_fun))
  #sti_med <- clusterR(fall_sti, calc, args=list(fun=med_fun))
  #max_ndvi <- clusterR(fall_ndvi, calc, args=list(max, na.rm=TRUE))
  #min_ndvi <- clusterR(fall_ndvi, calc, args=list(min, na.rm=TRUE))
  full_ndvi_max <- clusterR(ndvi_full, calc, args=list(max, na.rm=TRUE))
  
  endCluster()
  gc()
  #calculate ndvi amplitude and ratio
  cl <- makeCluster(clustnum)
  registerDoParallel(cl)
  ndvi_amp <- full_ndvi_max - ndvi_max
  ndvi_ratio <- ndvi_max / full_ndvi_max
  stopCluster(cl)
  
  #combine all of the RF layers together into a single stack
  stack_r <- stack(b3_med, b3_min, b3_max,
                   b5_med, b5_min, b5_max,
                   b6_med, b6_min, b6_max,
                   b10_med, b10_min, b10_max,
                   b11_med, b11_min, b11_max,
                   ndvi_med, ndvi_min, ndvi_max,
                   ndti_med, ndti_min, ndti_max,
                   sti_med, sti_min, sti_max,
                   ExG_med, ExG_min, ExG_max,
                   ExGR_med, ExGR_min, ExGR_max,
                   full_ndvi_max, ndvi_amp, ndvi_ratio)
  
  names(stack_r) <- c("b3_med", "b3_min", "b3_max", 
                      "b5_med", "b5_min", "b5_max",
                      "b6_med", "b6_min", "b6_max",
                      "b10_med", "b10_min", "b10_max",
                      "b11_med", "b11_min", "b11_max",
                      "ndvi_med", "ndvi_min", "ndvi_max",
                      "ndti_med", "ndti_min", "ndti_max",
                      "sti_med", "sti_min", "sti_max",
                      "exG_med", "exG_min", "exG_max",
                      "exGR_med", "exGR_min", "exGR_max",
                      "full_ndvi_max", "ndvi_amp", "ndvi_ratio")
  
  writeRaster(stack_r,paste("/N/project/CoverCrops/CoverCrops/data/spatial/RFstacks/HLS_stacks/cloudmsk/RFstack_months_", tile, "_", year, ".grd", sep=""), format="raster", overwrite = TRUE)
  print(paste("done with", tile, year, sep=" "))
}



setwd("/N/project/CoverCrops/CoverCrops/HLS_L30_v1.4/stacks/cloud_mask/")


rasterOptions(tmpdir="/N/slate/farellam/", tmptime = 1.25, progress="", timer=TRUE, overwrite = T, chunksize=8e10, maxmemory=5e8)

clustnum = 12

make_stacks(tile="16SDG", year = 2017)
make_stacks(tile="16SDG", year = 2018)
make_stacks(tile="16SDG", year = 2019)




