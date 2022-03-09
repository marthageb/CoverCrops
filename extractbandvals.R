###ALTERNATE CODE TO 'HLS_RFstacks.R' Instead of doing VIs and summary statistic calculations on raster stacks, band reflectance values extracted for transect locations and THEN VIs and summary statistics calculated on extracted .csv values.
#Steps:
  # 1.	raster stacks (output from ‘process_HLS.R’ loaded and band values (B2, B3, B4, B5, B6, B7, B10, B11, QA) extracted for transect locations
  # 2.	NDVI calculated for the full year
  # 3.	Band values subset for the CC/off-season – Nov  March
  # 4.	VIs calculated for CC/off-season (NDTI, STI, Rn, Gn, Bn, ExG, ExGR)
  # 5.	min, med, and max values calculated for #s 4 & 2
  # 6.	GDD calculated for min, med, and max NDVI values (# of days since Nov. 1 to each of these dates)
  # 7.	Output for each tile/year written: data/RFdata/transect_extract/TILE_YEAR.csv
                      

library(rgdal)
library(raster)
library(exactextractr)
library(stringr)

#load transect locations .shp file
sites <- readOGR("/N/project/CoverCrops/CoverCrops/data/spatial/25mbuffs.shp")
sites.df <- as.data.frame(sites)
sites.df$ID <- rownames(sites.df)
sites.df <- sites.df[,-c(5,6)]
#ctransform sites crs to crs of other rater layers
mycrs <- CRS('+proj=utm +zone=16 +ellps=WGS84 +units=m +no_defs')
sites <- spTransform(sites, mycrs)

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

##Function to calculate min and max GDD
minmax_fxn <- function(y){
  data <- minmax[y,]
  data$minyr <- substr(data$mindate, 7, 10)
  data$mindoy <- str_sub(data$mindate, -3, -1)
  data$minorig <- as.Date(as.numeric(data$mindoy), origin=as.Date(paste(data$minyr, "-01-01", sep="")))
  data$min_GDD <- length(seq(from=data$date1, to=data$minorig, by='day'))-1
  
  data$maxyr <- substr(data$maxdate, 7, 10)
  data$maxdoy <- str_sub(data$maxdate, -3, -1)
  data$maxorig <- as.Date(as.numeric(data$maxdoy), origin=as.Date(paste(data$maxyr, "-01-01", sep="")))
  data$max_GDD <- length(seq(from=data$date1, to=data$maxorig, by='day'))-1
  
  gdd.results <- cbind(data$sampleID, data$min_GDD, data$max_GDD)
  return(gdd.results)
}

###function to extract row(s) that are used to calculate the median value
which.median = function(x) {
  if (length(na.omit(x)) %% 2 != 0) {
    which(x == median(x, na.rm=TRUE))
  }  else if (length(na.omit(x)) %% 2 == 0) {
    a = sort(x)[c(length(na.omit(x))/2, length(na.omit(x))/2+1)]
    c(which(x == a[1]), which(x == a[2]))
  }
}


##function to calculate median GDD 
gddmed_fxn <- function(i){
  samplID <- colnames(ndvi.t)[i]
  data <- as.data.frame(ndvi.t[,i])
  data$date <- rownames(data)
  colnames(data) <- c("ndvi", "date")
  data$year <- substr(data$date, 7, 10)
  data$doy <- str_sub(data$date, -3, -1)
  data$origdate <- as.Date(as.numeric(data$doy), origin=as.Date(paste(data$year, "-01-01", sep="")))
  medndvi <- data[which.median(data$ndvi),]
  date1 <- as.Date(paste(min(data$year), "-11-01", sep=""))
  if (nrow(medndvi) == 1){
    gdd <- length(seq(from=date1, to=medndvi$origdate, by='day'))-1
  } else{
    dates <- medndvi$origdate
    gdd1 <- length(seq(from=date1, to=dates[1], by='day'))-1
    gdd2 <- length(seq(from=date1, to=dates[2], by='day'))-1
    gdd <- (gdd1 + gdd2) / 2
  }
  gdd.results <- cbind(samplID, gdd)
  return(gdd.results)
}

##function to load raster stacks, assign names, extract values for transect locs#
stack_extract <- function(tile, year) {
  gc()
  print(paste("now processing", tile, year, sep=" "))
  
  setwd("/N/project/CoverCrops/CoverCrops/HLS_L30_v1.4/stacks/cloud_mask/")
  all <- list.files("/N/project/CoverCrops/CoverCrops/HLS_L30_v1.4/stacks/cloud_mask/")
  tile_files <- all[grep(tile, all)]
  tile_files <- tile_files[grep(year, tile_files)]
  
  b2 <- raster::stack(tile_files[grep("b2", tile_files)])
  b3 <- raster::stack(tile_files[grep("b3", tile_files)])
  b4 <- raster::stack(tile_files[grep("b4", tile_files)])
  b5 <- raster::stack(tile_files[grep("b5", tile_files)])
  b6 <- raster::stack(tile_files[grep("b6", tile_files)])
  b7 <- raster::stack(tile_files[grep("b7", tile_files)])
  b10 <- raster::stack(tile_files[grep("b10", tile_files)])
  b11 <- raster::stack(tile_files[grep("b11", tile_files)])
  
  names(b2) <- paste(str_sub(names_fxn(tile=tile, year=year), -16, -13), str_sub(names_fxn(tile=tile, year=year), -12, -10), sep="_")
  names(b3) <- paste(str_sub(names_fxn(tile=tile, year=year), -16, -13), str_sub(names_fxn(tile=tile, year=year), -12, -10), sep="_")
  names(b4) <- paste(str_sub(names_fxn(tile=tile, year=year), -16, -13), str_sub(names_fxn(tile=tile, year=year), -12, -10), sep="_")
  names(b5) <- paste(str_sub(names_fxn(tile=tile, year=year), -16, -13), str_sub(names_fxn(tile=tile, year=year), -12, -10), sep="_")
  names(b6) <- paste(str_sub(names_fxn(tile=tile, year=year), -16, -13), str_sub(names_fxn(tile=tile, year=year), -12, -10), sep="_")
  names(b7) <- paste(str_sub(names_fxn(tile=tile, year=year), -16, -13), str_sub(names_fxn(tile=tile, year=year), -12, -10), sep="_")
  names(b10) <- paste(str_sub(names_fxn(tile=tile, year=year), -16, -13), str_sub(names_fxn(tile=tile, year=year), -12, -10), sep="_")
  names(b11) <- paste(str_sub(names_fxn(tile=tile, year=year), -16, -13), str_sub(names_fxn(tile=tile, year=year), -12, -10), sep="_")
  
  #extract Band Vals for transect locs
  print('Extracting Band Vals')
  b2vals <- exact_extract(b2, sites, 'mean')
  b3vals <- exact_extract(b3, sites, 'mean')
  b4vals <- exact_extract(b4, sites, 'mean')
  b5vals <- exact_extract(b5, sites, 'mean')
  b6vals <- exact_extract(b6, sites, 'mean')
  b7vals <- exact_extract(b7, sites, 'mean')
  b10vals <- exact_extract(b10, sites, 'mean')
  b11vals <- exact_extract(b11, sites, 'mean')
  
  gc()
  #Calcualte fullyr NDVI
  fullndvi <- (b5vals-b4vals)/(b5vals+b4vals)
  
  #set NDVI < 0 to NA to remove water pixels and > 1 to NA 
  fullndvi[fullndvi<0] <- NA
  fullndvi[fullndvi>1] <- NA
  
  
  #subset the remaining bands to the off season (Nov - March)
  doy <- str_sub(colnames(b4vals), -3, -1)
  yr <- str_sub(colnames(b4vals), -8, -5)
  dates <- cbind.data.frame(yr, doy)
  curr <- subset(dates, yr==year)
  nxt <- subset(dates, yr==year+1)
  nxt <- subset(nxt, as.numeric(doy) < 90)
  yrs <- unlist(list(curr$yr, nxt$yr))
  days <- unlist(list(curr$doy, nxt$doy))
  wanted <- paste("mean.X", yrs, "_", days, sep="")
  col.num <- which(colnames(b4vals) %in% wanted)
  b2off <- b2vals[,col.num]
  b3off <- b3vals[,col.num]
  b4off <- b4vals[,col.num]
  b5off <- b5vals[,col.num]
  b6off <- b6vals[,col.num]
  b7off <- b7vals[,col.num]
  b10off <- b10vals[,col.num]
  b11off <- b11vals[,col.num]
  ndvioff <- fullndvi[,col.num]
  
  #calculate other VIs
  print('calculating VIs')
  ndti <- (b6off - b7off)/(b6off + b7off)
  sti <- b6off/b7off
  Rn <- b4off/(b2off+b3off+b4off)
  Gn <- b3off/(b2off+b3off+b4off)
  Bn <- b2off/(b2off+b3off+b4off)
  ExG <- (2*Gn) - Rn - Bn
  ExGR <- ExG - (1.4*Rn) - Gn
  
  #create a matrix to store the final vals
  final <- data.frame(nrow=nrow(b4off),0)
  
  #write function to assign NA in max calculation if all values in row are NA
  rowMax <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)
  rowMin <- function(x) ifelse( !all(is.na(x)), min(x, na.rm=T), NA)
  rowMed <- function(x){median(x, na.rm=TRUE)}
  
  #calcualte min, max, and med values
  print('calculating min, max, median vals')
  b3min <- apply(b3off, 1, rowMin)
  b3max <- apply(b3off, 1, rowMax)
  b3med <- apply(b3off, 1, rowMed)
  b5min <- apply(b5off, 1, rowMin)
  b5max <- apply(b5off, 1, rowMax)
  b5med <- apply(b5off, 1, rowMed)
  b6min <- apply(b6off, 1, rowMin)
  b6max <- apply(b6off, 1, rowMax)
  b6med <- apply(b6off, 1, rowMed)
  b10min <- apply(b10off, 1, rowMin)
  b10max <- apply(b10off, 1, rowMax)
  b10med <- apply(b10off, 1, rowMed)
  b11min <- apply(b11off, 1, rowMin)
  b11max <- apply(b11off, 1, rowMax)
  b11med <- apply(b11off, 1, rowMed)
  ndvimin <- apply(ndvioff, 1, rowMin)
  ndvimax <- apply(ndvioff, 1, rowMax)
  ndvimed <- apply(ndvioff, 1, rowMed)
  ndtimin <- apply(ndti, 1, rowMin)
  ndtimax <- apply(ndti, 1, rowMax)
  ndtimed <- apply(ndti, 1, rowMed)
  stimin <- apply(sti, 1, rowMin)
  stimax <- apply(sti, 1, rowMax)
  stimed <- apply(sti, 1, rowMed)
  ExGmin <- apply(ExG, 1, rowMin)
  ExGmax <- apply(ExG, 1, rowMax)
  ExGmed <- apply(ExG, 1, rowMed)
  ExGRmin <- apply(ExGR, 1, rowMin)
  ExGRmax <- apply(ExGR, 1, rowMax)
  ExGRmed <- apply(ExGR, 1, rowMed)
  full_ndvimax <- apply(fullndvi, 1, rowMax)
  ndviamp <- full_ndvimax - ndvimax
  ndviratio <- full_ndvimax / ndvimax
  
  #bind all rows together
  final <- cbind.data.frame(sites.df, b3min, b3max, b3med, b5min, b5max, b5med, b6min, b6max, b6med, b10min, b10max, b10med, b11min, b11max, b11med,
                            ndvimin, ndvimax, ndvimed, ndtimin, ndtimax, ndtimed, stimin, stimax, stimed, ExGmin, ExGmax, ExGmed, ExGRmin, ExGRmax, ExGRmed,
                            full_ndvimax, ndviamp, ndviratio)
  
  rownames(ndvioff) <- rownames(final)
  #remove NA observations
  final <- final[complete.cases(final),]
  final$sampleID <- rownames(final)
  
  #subset ndvi to only the rows with observations
  ndvigood <- ndvioff[rownames(final),]
  
  #min and max GDD
  print('getting min and max GDD')
  
  mindate <- apply(ndvigood, 1, function(x) colnames(ndvigood[which.min(x)]))
  maxdate <- apply(ndvigood, 1, function(x) colnames(ndvigood[which.max(x)]))
  
  minmax <- cbind.data.frame(mindate, maxdate)
  minmax$sampleID <- rownames(minmax)
  minmax$date1 <- as.Date(paste(year, "-11-01", sep=""))
  
  rows <- seq(1:nrow(minmax))
  minmaxgdd <- do.call(rbind.data.frame, lapply(rows, minmax_fxn))
  names(minmaxgdd) <- c('sampleID', "min_GDD", "max_GDD")
  
  
  #median GDD
  print('getting median GDD')
  
  ndvi.t <- t(ndvigood)
  cols <- seq(1:ncol(ndvi.t))
  medgdd <- do.call(rbind.data.frame, lapply(cols, gddmed_fxn))
  names(medgdd) <- c('sampleID', "med_GDD")
  
  final <- merge(final, minmaxgdd, by='sampleID')
  final <- merge(final, medgdd, by="sampleID")
  
  final$year <- year
  write.csv(paste("/N/project/CoverCrops/CoverCrops/data/RFdata/transect_extract/", tile, "_", year, ".csv", sep=""))
  return(final)
}


#apply the stack extract function
stack_extract(tile="16TEL", year = 2014)
stack_extract(tile="16TEL", year = 2015)
stack_extract(tile="16TEL", year = 2016)
stack_extract(tile="16TEL", year = 2017)
stack_extract(tile="16TEL", year = 2018)
stack_extract(tile="16TEL", year = 2019)


