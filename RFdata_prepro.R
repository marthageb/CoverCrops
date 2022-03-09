###Get the data formatted for RF analysis. Extract raster vals for the LST and HLS stacks for sample locs. Then, combine extracted vals with windshield survey data and write output.



library(dplyr)
library(stringr)
library(rgdal)
library(raster)
library(sf)

#set temp dir to slate folder
#install.packages('unixtools',repos='http://www.rforge.net/')
unixtools::set.tempdir("/N/slate/farellam")
tempdir()
#####################################################################################################
##########################extract LST values#########################################################
#####################################################################################################
sites <- readOGR("/N/project/CoverCrops/CoverCrops/data/spatial/25mbuffs.shp")
sites.df <- as.data.frame(sites)
sites.df$ID <- rownames(sites.df)
sites.df <- sites.df[,-c(5,6)]


#make a list of the LST rasters that overlap windshield survey observations-- see pg. 142 of black notebook
overlap <- c("2014_021008", "2014_021009", "2014_022007", "2014_022008", "2014_022009", "2014_022010", "2014_023008", "2014_023009", 
             "2015_021008", "2015_021009", "2015_021010", "2015_022007", "2015_022008", "2015_022009", "2015_022010", "2015_023008", "2015_023009",
             "2016_021008", "2016_021009", "2016_021010", "2016_022007", "2016_022008", "2016_022009", "2016_022010", "2016_023008", "2016_023009",
             "2017_022007", "2017_022008", "2017_022009", "2017_022010", "2017_023008", "2017_023009",
             "2018_022007", "2018_022008", "2018_022009", "2018_022010", "2018_023008", "2018_023009",
             "2019_022007", "2019_022008", "2019_022009", "2019_022010", "2019_023008")

files <- list.files("/N/project/CoverCrops/CoverCrops/HLS30_Indiana/", full.names = TRUE)

#select only the median or sd rasters
med_files <- files[grep("LSTmed", files)]
sd_files <- files[grep("LSTsd", files)]

#keep only the desired tiles
med_files <- subset(med_files, str_sub(med_files, -15, -5) %in% overlap)
sd_files <- subset(sd_files, str_sub(sd_files, -15, -5) %in% overlap)

#create rasters
med_r <- lapply(med_files, raster)
sd_r <- lapply(sd_files, raster)

#transform sites to the same projection as the LST rasters
sites <- spTransform(sites, crs(med_r[[1]]))

#function to extract median LST values 
extract_med <- function(x) {
  LSTval <- raster::extract(x, sites, fun=mean, df=TRUE)
  LSTval <- LSTval[complete.cases(LSTval),]
  LSTdf <- merge(sites.df, LSTval, by='ID')
  id <- strsplit(names(LSTdf[6]), split="_")
  LST_ID <- paste(id[[1]][2], id[[1]][3], sep="_")
  #LST_ID <- paste(id[[1]][2], str_sub(id[[1]][1], -4,-1), sep="_")
  if (nrow(LSTdf) > 0) {
    LSTdf$LST_ID <- LST_ID
    names(LSTdf)[6] <- "LST_med"
    print(unique(LSTdf$LST_ID))
    return(LSTdf)
  }
}
med_results <- lapply(med_r, extract_med) %>% bind_rows()
#med_final <- med_results

    #need to reshape data so that different monthly means for each fall are in separate columns-- ONLY NEED TO DO THIS IS UNSING MONTHLY OFF-SEASON MEDIAN VALS (VS. TOTAL OFF-SEASON MEDIAN VALS)
    med_results$month <- sapply(strsplit(med_results$LST_ID, split="_"), '[', 1)
    med_results$year <- as.numeric(str_sub(med_results$month, -4, -1))
    med_results$month <- substr(med_results$month, 1, 3)
    #update year values for feb and march to be the correct fall year
    med_results$year <- ifelse(med_results$month == "feb", med_results$year - 1, med_results$year)
    med_results$year <- ifelse(med_results$month == "mar", med_results$year - 1, med_results$year)
    med_results$tile <- sapply(strsplit(med_results$LST_ID, split="_"), '[', 2)
    med_results$LST_ID <- paste("fall", med_results$year, "_", med_results$tile, sep="")
    med_results <- med_results[,-c(9,10)]
    med_results <- reshape(med_results, idvar=c("ID", "site", "Field_No", "Latitude", "Longitude", "LST_ID"), timevar="month", direction='wide')
    
    
med_results$uid <- paste(med_results$site, med_results$Field_No, med_results$LST_ID, sep="")
med_results <- med_results[complete.cases(med_results),]


#function to extract sd LST values 
extract_sd <- function(x) {
  LSTval <- raster::extract(x, sites, fun=mean, df=TRUE)
  LSTval <- LSTval[complete.cases(LSTval),]
  LSTdf <- merge(sites.df, LSTval, by='ID')
  LSTdf$LST_ID <- substr(names(LSTdf[6]), 6, 20)
  names(LSTdf)[6] <- "LST_sd"
  print(unique(LSTdf$LST_ID))
  return(LSTdf)
}

sd_results <- lapply(sd_r, extract_sd) %>% bind_rows()
sd_results$uid <- paste(sd_results$site, sd_results$Field_No, sd_results$LST_ID, sep="")
sd_results <- sd_results[complete.cases(sd_results),]
names(sd_results)
names(med_results)

results <- merge(med_results, sd_results, by=c("uid", "ID", "site", "Field_No", "Latitude", "Longitude", "LST_ID"))
  
  
  
write.csv(results, "/N/project/CoverCrops/CoverCrops/data/RFdata/LST.csv")

#####################################################################################################
##############################extract other Landsat RF data #########################################
#####################################################################################################
sites <- readOGR("/N/project/CoverCrops/CoverCrops/data/spatial/25mbuffs.shp")
sites.df <- as.data.frame(sites)
sites.df$ID <- rownames(sites.df)
sites.df <- sites.df[,-c(5,6)]

#RFstacks <- list.files("/N/project/CoverCrops/CoverCrops/data/spatial/RFstacks/new_stacks", full.names = TRUE, pattern=".grd") #for the LANDSAT 8 C2L2 data
stacks <- list.files("/N/project/CoverCrops/CoverCrops/data/spatial/RFstacks/HLS_stacks/cloudmsk/", full.names = TRUE, pattern=".grd") # for the HLS data
#select the one value per month stacks
#stacks <- stacks[grep("_months_", stacks)]
#select a single tile
onetile <- "16TFL"
RFstacks <- (stacks[grep(onetile, stacks)])


#transform sites to the projection of the landsat rasters
crsrast <- raster(RFstacks[1])
sites <- spTransform(sites, crs(crsrast))  

clustnum <- 8

beginCluster(clustnum)

#function to extract values from raster data
##takes about 5 mins per file (48 stacks x 5 = 4hrs)
extractvals <- function(x) {
  year <- str_sub(x, -8, -5)
  tile <- str_sub(x, -14, -10)
  print(paste("now processing", tile, year, sep=" "))
  rast <- raster::stack(x)
  vals <- raster::extract(rast, sites, fun=mean, df=TRUE)
  vals2 <- vals[complete.cases(vals),]
  vals_df <- merge(sites.df, vals2, by='ID')
  if (nrow(vals_df) > 0) {
    vals_df$year <- str_sub(x, -8, -5)
    vals_df$tile <- str_sub(x, -15, -10)
    return(vals_df)
  }
}

results <- lapply(RFstacks, extractvals) %>% bind_rows()

endCluster()

nrow(results)
#write.csv(results, "/N/project/CoverCrops/CoverCrops/data/RFdata/TIRSC2L2_021032.csv")
write.csv(results, paste("/N/project/CoverCrops/CoverCrops/data/RFdata/HLS/cloudmsk_", onetile, ".csv", sep=""))
unique(results$site)



#combine all of the tiles into a single dataset
setwd("/N/project/CoverCrops/CoverCrops/data/RFdata/HLS/")
file_list <- list.files()
file_list <- file_list[grep("cloudmsk_", file_list)]
#file_list <- file_list[grep("mo_vals", file_list)]
for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- read.csv(file, header=TRUE)
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <-read.csv(file, header=TRUE)
    dataset<-rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
}






write.csv(dataset, "/N/project/CoverCrops/CoverCrops/data/RFdata/HLS/cloudmsk_alltiles.csv")

#####################################################################################################
#########################combine extracted RF predictor vars with survey data########################
#####################################################################################################
#load survey data (output from 'format_windshield.R')
survey <- read.csv("/N/project/CoverCrops/CoverCrops/data/surveydata.csv")
#change 'StJoseph' to 'StJoe'
survey$site[survey$site == "StJoseph"] <- 'StJoe'
#make unique ID column
survey$uid <- paste(survey$date, survey$site, survey$Field_No, sep="")
#get rid of rows with NA CoverCrop survey values
survey <- survey[!is.na(survey$Cover_Crop),]
#check to see id uid values are truly unique -- want value to be 0
anyDuplicated(survey$uid)


#load LST data 
lst1 <- read.csv("/N/project/CoverCrops/CoverCrops/data/RFdata/LST.csv")
lst2 <- read.csv("/N/project/CoverCrops/CoverCrops/data/RFdata/fall_singleval/LST.csv")
lst <- merge(lst1, lst2[,c(2,9,10)], by="uid")

#lst <- subset(lst, LST_ID== "fall2015_021010")
lst$date <- substr(lst$LST_ID, 1, 8)
#capitilize 'Fall'
lst$date <- str_to_title(lst$date)
#make unique ID column
lst$uid <- paste(lst$date, lst$site, lst$Field_No, sep="")

#get average values when there's more than one unique ID
RFlst <- lst %>%
  group_by(uid) %>%
  summarize(sdLST = mean(LST_sd.x), med_fullLST = mean(LST_med),
            med_novLST = mean(LST_med.nov), med_wintLST = mean(LST_med.win),
            med_febLST = mean(LST_med.feb), med_marLST = mean(LST_med.mar))

#when doing single val in off season
RFlst <- lst %>% group_by(uid) %>% summarize(sdLST = mean(LST_sd), medLST = mean(LST_med)) #when not doing monthly vals

#load other Landsat predictor vars
pred <- read.csv("/N/project/CoverCrops/CoverCrops/data/RFdata/HLS/cloudmsk_alltiles.csv")

pred$date <- paste("Fall", pred$year, sep="")
#make uniqueID column
pred$uid <- paste(pred$date, pred$site, pred$Field_No, sep="")

colnames(pred) <- sub("mean.", "", colnames(pred)) 

#get average values when there's more than one unique ID
RFpred <- pred %>%
  group_by(uid) %>%
  summarize(medB3 = mean(b3_med), minB3 = mean(b3_min), maxB3 = mean(b3_max),
            medB5 = mean(b5_med), minB5 = mean(b5_min), maxB5 = mean(b5_max),
            medB6 = mean(b6_med), minB6 = mean(b6_min), maxB6 = mean(b6_max),
            medB10 = mean(b10_med), minB10 = mean(b10_min), maxB10 = mean(b10_max),
            medB11 = mean(b11_med), minB11 = mean(b11_min), maxB11 = mean(b11_max),
            medNDVI = mean(ndvi_med), minNDVI = mean(ndvi_min), maxNDVI = mean(ndvi_max),
            medNDTI = mean(ndti_med), minNDTI = mean(ndti_min), maxNDTI = mean(ndti_max),
            medSTI = mean(sti_med), minSTI = mean(sti_min), maxSTI = mean(sti_max),
            medExG = mean(exG_med), minExG = mean(exG_min), maxExG = mean(exG_max),
            medExGR = mean(exGR_med), minExGR = mean(exGR_min), maxExGR = mean(exGR_max),
            full_maxNDVI = mean(full_ndvi_max), ampNDVI = mean(ndvi_amp), ratioNDVI = mean(ndvi_ratio))
            

#for monthly vals
RFpred <- pred %>% group_by(uid) %>% summarize(b3_nov=mean(b3_nov), b3_wint=mean(b3_wint), b3_feb=mean(b3_feb), b3_mar= mean(b3_mar),
                                               b5_nov=mean(b5_nov), b5_wint=mean(b5_wint), b5_feb=mean(b5_feb), b5_mar= mean(b5_mar),
                                               b6_nov=mean(b6_nov), b6_wint=mean(b6_wint), b6_feb=mean(b6_feb), b6_mar= mean(b6_mar),
                                               b10_nov=mean(b10_nov), b10_wint=mean(b10_wint), b10_feb=mean(b10_feb), b10_mar= mean(b10_mar),
                                               b11_nov=mean(b11_nov), b11_wint=mean(b11_wint), b11_feb=mean(b11_feb), b11_mar= mean(b11_mar),
                                               ndvi_nov=mean(ndvi_nov), ndvi_wint=mean(ndvi_wint), ndvi_feb=mean(ndvi_feb), ndvi_mar= mean(ndvi_mar),
                                               ndti_nov=mean(ndti_nov), ndti_wint=mean(ndti_wint), ndti_feb=mean(ndti_feb), ndti_mar= mean(ndti_mar),
                                               sti_nov=mean(sti_nov), sti_wint=mean(sti_wint), sti_feb=mean(sti_feb), sti_mar= mean(sti_mar))
#combine all of the datasets together
RFdata <- merge(RFlst, RFhls, by="uid")
RFdata <- merge(survey, RFdata ,by="uid")
names(RFdata)

unique(RFdata$site)
#data2 <- subset(RFdata ,site %in% c("Benton", "Gibson", "Posey", "Warren", "White"))
#data2 <- subset(data2, date=="Fall2015")

#write the dataset
write.csv(RFdata, "/N/project/CoverCrops/CoverCrops/data/RFdata/cloudmsk_RFdata_25mbuff.csv")



