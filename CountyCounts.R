####Get pixel counts of CC, Conv, and NoTill from the RFmodels at the county level for agriculture areas only.
    ## 1.	Process the NLCD raster—crop to the extent of the RFmodels and resample to allow for the two data products to align
    ## 2.	Load the 2- and 3-class RF predictions (output from applyRF.R), mask to agriculture areas only (output from #1), mask to the county level.
    ## 3.	Get pixel counts of the RF 2- and 3-class classifications 
    ## 4. calculate acreage for pixel counts, combine 2- and 3-class CC acreage for each county/year into a single dataframe, combine with ISDA estimates, and write output.
                                                                                                   

library(raster)
library(rgdal)
library(stringr)
library(dplyr)
library(plyr)

#load the nlcd landcover raster
nlcd <- raster("/N/project/CoverCrops/CoverCrops/data/spatial/NLCD/NLCD2016_HLScrop.tif")
#Croplands == 81,82
nlcd[nlcd<80] <- NA


countyco_fxn <- function(county, year) {
  #load the RF prediction
  rf2class <- raster(paste("/N/project/CoverCrops/CoverCrops/data/spatial/RFresults/2class/", county, year, ".tif", sep=""))
  rf3class <- raster(paste("/N/project/CoverCrops/CoverCrops/data/spatial/RFresults/3class/", county, year, ".tif", sep=""))
  
  #load the county shapefile
  copoly <- readOGR(paste("/N/project/CoverCrops/CoverCrops/data/spatial/INcounties/indv_cos/", county, ".shp", sep=""))
  copoly_t <- spTransform(copoly, crs(rf2class)) #transform to HLS CRS
  
  #crop the nlcd raster to the extent of the rf prediction
  nlcd_crop <- crop(nlcd, rf2class)
  #Then we need to resample the nlcd raster so extents can match up
  cropmsk <- resample(nlcd_crop, rf2class, method='ngb')
  
  #now, mask the RF prediction to ag areas only
  rfag2class <- mask(rf2class, cropmsk)
  rfag3class <- mask(rf3class, cropmsk)
  
  #then, mask the prediction to the county
  rfco2class <- mask(rfag2class, copoly_t)
  rfco3class <- mask(rfag3class, copoly_t)
  
  plot(rfco2class) #look at the data
  plot(rfco3class) #look at the data
  writeRaster(rfco2class, paste("/N/project/CoverCrops/CoverCrops/data/spatial/RFresults/RF2class_", county, year, ".tif",sep=""))
  writeRaster(rfco3class, paste("/N/project/CoverCrops/CoverCrops/data/spatial/RFresults/RF3class_", county, year, ".tif", sep=""))
              
  
  #Then summarize the number of pixels classified as each class for ag areas only in whichever county.
  freq2 <- as.data.frame(freq(rfco2class))
  freq3 <- as.data.frame(freq(rfco3class))
  
  
  RF2class <- readRDS("/N/project/CoverCrops/CoverCrops/RFmodels/RF2class_Conv_CC.rds")
  levels2 <- rbind(as.data.frame(RF2class$levels), "NA")
  
  RF3class <- readRDS("/N/project/CoverCrops/CoverCrops/RFmodels/RF3class_HLS_cloudmsk_25mbuff.rds")
  levels3 <- rbind(as.data.frame(RF3class$levels), "NA")
  
  
  final2 <- cbind.data.frame(levels2, freq2)
  final2 <- final2[,-2]
  names(final2) <- c("class", "count")
  print(final2)
  write.csv(final2, paste("/N/project/CoverCrops/CoverCrops/RFCoCounts/2class_", county, year, ".csv", sep=""))
  
  final3 <- cbind.data.frame(levels3, freq3)
  final3 <- final3[,-2]
  names(final3) <- c("class", "count")
  print(final3)
  write.csv(final3, paste("/N/project/CoverCrops/CoverCrops/RFCoCounts/3class_", county, year, ".csv", sep=""))
}

countyco_fxn(county='decatur', year= 2019)
countyco_fxn(county='gibson', year= 2014)
countyco_fxn(county='gibson', year= 2017)
countyco_fxn(county='gibson', year= 2018)
countyco_fxn(county='gibson', year= 2019)
countyco_fxn(county='hamilton', year= 2014)
countyco_fxn(county='stjoe', year= 2019)


#######################################################################################################################
#################Combine county count for each year and model into a single data frame#######################
#######################################################################################################################

##function to calculate CC acreage from pixel #
onedf <- function(file){
  data <- read.csv(file)
  year <- str_sub(file, -8, -5)
  prediction <- substr(file, 1, 6)
  #get the county
  county <- str_split(file, pattern="_")
  county <- county[[1]][[2]]
  county <- str_split(county, pattern="20")
  county <- county[[1]][[1]]
  
  #calculate the acreage of the CC prediction
  cc <- subset(data, class=="CC")
  area <- cc$count * 30^2
  #from the area, calculate the acreage
  acres <- area / 4046.8564224
  
  countyacres <- cbind.data.frame(year, prediction, acres)
  names(countyacres)[names(countyacres) == 'acres'] <- county
  
  return(countyacres)
}

setwd("/N/project/CoverCrops/CoverCrops/RFCoCounts/")
#list files created in the 'countyco_fxn'
allcos <- list.files("/N/project/CoverCrops/CoverCrops/RFCoCounts/")

#subset for individual counties and apply the function
benton_co <- allcos[grep('benton', allcos)]
benton <- lapply(benton_co, onedf) %>% bind_rows()

decatur_co <- allcos[grep('decatur', allcos)]
decatur <- lapply(decatur_co, onedf) %>% bind_rows()

gibson_co <- allcos[grep('gibson', allcos)]
gibson <- lapply(gibson_co, onedf) %>% bind_rows()

hamilton_co <- allcos[grep('hamilton', allcos)]
hamilton <- lapply(hamilton_co, onedf) %>% bind_rows()

posey_co <- allcos[grep('posey', allcos)]
posey <- lapply(posey_co, onedf) %>% bind_rows()

stjoe_co <- allcos[grep('stjoe', allcos)]
stjoe <- lapply(stjoe_co, onedf) %>% bind_rows()

warren_co <- allcos[grep('warren', allcos)]
warren <- lapply(warren_co, onedf) %>% bind_rows()

washington_co <- allcos[grep('washington', allcos)]
washington <- lapply(washington_co, onedf) %>% bind_rows()

white_co <- allcos[grep('white', allcos)]
white <- lapply(white_co, onedf) %>% bind_rows()

whitley_co <- allcos[grep('whitley', allcos)]
whitley <- lapply(whitley_co, onedf) %>% bind_rows()

#then, combine all of the individual counties into a single df
all <- join_all(list(benton, decatur,gibson, hamilton, posey, stjoe, warren, washington, white, whitley))

#load Indiana State Dept of Ag CC estimates
isda <- read.csv("/N/project/CoverCrops/CoverCrops/data/ISDA_CCbyCo.csv")
isda$year <- str_sub(isda$siteyear, -4, -1)
isda$county <- gsub('[[:digit:]]+', '', isda$siteyear)

isda <- isda[,-1]
isda_wide <- reshape(isda, idvar="year", timevar="county", direction="wide")

colnames(isda_wide) <- gsub('ISDA.', '', colnames(isda_wide))

isda_wide$prediction <- "ISDA"

final <- rbind.data.frame(all, isda_wide)


write.csv(final, "/N/project/CoverCrops/CoverCrops/data/CCpred_yearxco.csv")
