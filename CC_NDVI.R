####NDVI analysis for the cover crop season for the different cover crop types. 
    ## 1. ndviextract: function to extract B4 and B5 vals for the transect locs., calculate NDVI, subset NDVI values to the cover crop season (Nov- March)
    ## 2. add uid with site and date, aggregate ndvi values by uid, combine ndvi values with CC type
    ## 3. plotting: for each CC type, gap fill missing NDVI values, calcualte CC type average NDVI and sd vals, plot data



library(rgdal)
library(raster)
install.packages("exactextractr", repos = "http://cran.us.r-project.org")
library(exactextractr)
library(stringr)
library(tibble)
library(extrafont)
library(zoo)
library(dplyr)
library(ggplot2)

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

##function to load raster stacks, assign names, extract values for transect locs, and calcualte NDVI#
ndviextract <- function(tile, year) {
  gc()
  print(paste("now processing", tile, year, sep=" "))
  
  setwd("/N/project/CoverCrops/CoverCrops/HLS_L30_v1.4/stacks/cloud_mask/")
  all <- list.files("/N/project/CoverCrops/CoverCrops/HLS_L30_v1.4/stacks/cloud_mask/")
  tile_files <- all[grep(tile, all)]
  tile_files <- tile_files[grep(year, tile_files)]
  
  b4 <- raster::stack(tile_files[grep("b4", tile_files)])
  b5 <- raster::stack(tile_files[grep("b5", tile_files)])
 
  names(b4) <- paste(str_sub(names_fxn(tile=tile, year=year), -16, -13), str_sub(names_fxn(tile=tile, year=year), -12, -10), sep="_")
  names(b5) <- paste(str_sub(names_fxn(tile=tile, year=year), -16, -13), str_sub(names_fxn(tile=tile, year=year), -12, -10), sep="_")
  
  
  #extract Band Vals for transect locs
  print('Extracting Band Vals')
  b4vals <- exact_extract(b4, sites, 'mean')
  b5vals <- exact_extract(b5, sites, 'mean')
  
  gc()
  #Calcualte fullyr NDVI
  fullndvi <- (b5vals-b4vals)/(b5vals+b4vals)
  
  #set NDVI < 0 to NA to remove water pixels and > 1 to NA 
  fullndvi[fullndvi<0] <- NA
  fullndvi[fullndvi>1] <- NA
  
  #subset the ndvi to the off season (Nov - March)
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
  ndvioff <- fullndvi[,col.num]
  
  #assign site.df row names to the ndvioff df
  rownames(ndvioff) <- rownames(sites.df)
  
  #remove observations with all NA values
  ndviobs <- ndvioff[rowSums(is.na(ndvioff)) != ncol(ndvioff), ]
  
  #subset sites to only the rows with ndvi observations
  sitesobs <- sites.df[rownames(ndviobs),]
  
  #combine extracted NDVI values with site info
  final <- cbind.data.frame(sitesobs, ndviobs)
  
  #save the dataframe
  write.csv(final, paste("/N/project/CoverCrops/CoverCrops/data/NDVI/", tile, "_", year, ".csv", sep=""))
}

#apply the stack extract function
ndviextract(tile="16TDL", year = 2014)
ndviextract(tile="16TDL", year = 2015)
ndviextract(tile="16TDL", year = 2016)
ndviextract(tile="16TDK", year = 2014)
ndviextract(tile="16TDK", year = 2015)
ndviextract(tile="16TDK", year = 2016)

ndviextract(tile="16SDH", year = 2015)
ndviextract(tile="16SDH", year = 2016)

ndviextract(tile="16SDG", year = 2015)
ndviextract(tile="16SDG", year = 2016)



###############################################################################
#####################GET DATA READY FOR PLOTTING##############################
##load extracted ndvi data, 
setwd("/N/project/CoverCrops/CoverCrops/data/NDVI/") 


files <- list.files("/N/project/CoverCrops/CoverCrops/data/NDVI/")
files <- files[grep("16", files)]
#oneyr <- files[grep(year, files)]

#make an empty list to hold outputs
final <- list()
alldates <- c()

#read each file, create UID column, also get dates from column names
for (i in files){
  print(i)
  x <- read.csv(i)
  dates <- colnames(x)
  dates <- dates[grep('mean.X', dates)]
  years <- substr(dates, 7,10)
  minyear <- min(as.numeric(years))
  uid <- paste("Fall", minyear, x$site, x$Field_No, sep="")
  x <- cbind.data.frame(uid, x)
  final[[i]] <- x
  alldates <- c(alldates, dates)
}

#get the dates from all the files
dates <- unique(alldates)

#add a column of NAs when date doesn't exist
adddate <- function(file){
  add_column(file, !!!dates[setdiff(names(dates), names(file))])
  
}
#apply the 'adddate' function and bind all results together into a single dataset
cordates <- lapply(final, adddate) %>% bind_rows()
#get rid of unwanted columns
cordates <- cordates[,c(1,4:5, 7:196)]

#get average values when there's more than one unique ID
by.uids <- cordates %>% group_by(uid)
ndvi <- by.uids %>%  summarize_each(funs(mean(., na.rm=TRUE)))


#load the RF data
rfdata <- read.csv("/N/project/CoverCrops/CoverCrops/data/RFdata/cloudmsk_RFdata_25mbuff.csv")
rfdata[sapply(rfdata, is.infinite)] <- NA
rfdata <- rfdata[rowSums(is.na(rfdata[,17:51]))!=35,]
#select only wanted columns
rfdata <- rfdata[c(2, 12, 14)]

#make RF classes
rfdata$RF2class <- "CC"
rfdata$RF2class[rfdata$Cover_Crop == "N"] <- "NoCC" 
rfdata$RF2class <- as.factor(rfdata$RF2class)

#make empty column to hold 3class classes
rfdata$RF3class <- NA
#determine if fields had No Till 
rfdata$RF3class[rfdata$Tillage == 'N' | rfdata$Tillage == 'n'] <- 'NoTill'
#determine if fields had Conventional Tillage
rfdata$RF3class[rfdata$Tillage == 'C' | rfdata$Tillage == 'c'] <- 'ZConv'
#if a field had CC in the 2 class model it should also have CC in the 3 class model
rfdata$RF3class[rfdata$RF2class == 'CC'] <- 'CC'
rfdata$RF3class <- as.factor(rfdata$RF3class)
levels(as.factor(rfdata$RF3class))

#look at CCtype
cc <- subset(rfdata, RF3class=="CC")
summary(as.factor(cc$Cover_Crop))

#define CC type
cc$cc.agg <- cc$Cover_Crop
cc$cc.agg[cc$Cover_Crop %in% c("0", "0B", "5N", "AB", "AL", "ann", "AP", "ba", "BC", "BG", "BLL", "BLS",
                               "bo", "BO", "BP", "BW", "cb", "CL", "CP", "CR", "GP", "H", "l", "L", "LO",
                               "O", "OB", "ob", "OB", "OBS", "OR", "Q", "R", "U", "X", "Y", "P")] <- 'Mix/Other'
cc$cc.agg[cc$Cover_Crop %in% c("a", "A")] <- "RyeGrass"
cc$cc.agg[cc$Cover_Crop == "B"] <- "Brassica"
cc$cc.agg[cc$Cover_Crop %in% c("c", "C")] <- "CerealRye"
cc$cc.agg[cc$Cover_Crop == "G"] <- "WinterGrain"
#cc$cc.agg[cc$Cover_Crop == "P"] <- "Mix"
cc$cc.agg[cc$Cover_Crop %in% c("w", "W")] <- "Wheat"
summary(as.factor(cc$cc.agg))

#get no till and conv
other <- subset(rfdata, RF3class %in% c("NoTill", "ZConv"))
other$cc.agg <- other$RF3class

rfall <- rbind.data.frame(cc, other)


plotdata <- merge(rfall, ndvi, by='uid')

plotdata$seas <- substr(plotdata$uid, 1, 8)

###########################PLOTTING############################
#make a dataframe with CC type and plotting colors
alltypes <- unique(plotdata$cc.agg)
cols <- c("#69585F", "#CCDDD3", "#EF959D", "#60656F", "#FFF9A5", "#B48B7D", "#AFC2D5", "#DFEFCA")
cols.df <- cbind.data.frame(alltypes, cols)

plotndvi <- function(cctype, season){
  #subset a single CC type
  onetyp <- subset(plotdata, cc.agg == cctype)
  
  #for a single year
  onetyp <- subset(onetyp, seas== season)
  
  alldates <- colnames(onetyp)
  alldates <- as.data.frame(alldates[grep("mean.X", alldates)])
  names(alldates) <- 'colnames'
  alldates$year <- as.numeric(substr(alldates$colnames, 7, 10))
  alldates$doy <- as.numeric(substr(alldates$colnames, 12, 14))
  thisyr <- as.numeric(substr(season, 5, 8))
  current <- subset(alldates, year == thisyr)
  current <- subset(current, doy > 95)
  nextyr <- subset(alldates, year == thisyr +1)
  nextyr <- subset(nextyr, doy < 95)
  allyear <- rbind.data.frame(current, nextyr)
  
  ids <- onetyp[,c(1:8)]
  other <- onetyp[,c(9:199)]
  other <-other[,allyear$colnames]
  onetyp <- cbind.data.frame(ids, other)

  
  #get rid of observations with less than 4 observations from Nov - March
  onetyp_dates <- onetyp[,9:ncol(onetyp)]
  typgood <- onetyp[ncol(onetyp_dates) - rowSums(is.na(onetyp[,9:ncol(onetyp)])) > 4,]
  
      ####gap filling missing NDVI values####
      #create an empty dataframe to hold output
      gapfilled <- NULL
       
      for (i in 1:nrow(typgood)){
        #select a single observation
        oneobs <- typgood[i,]
        uid <- oneobs$uid
        oneobs <- oneobs[1,9:(ncol(oneobs))]
        oneobs <- as.data.frame(t(oneobs))
        names(oneobs) <- 'V1'
        #get the date
        oneobs$date <- substr(rownames(oneobs), 7, 14)
        oneobs$year <- substr(oneobs$date, 1, 4)
        oneobs$doy <- substr(oneobs$date, 6, 8)
        oneobs$origdate <- as.Date(as.numeric(oneobs$doy), origin=as.Date(paste(oneobs$year, "-01-01", sep="")))
        
        #select a single year -- probably want to do this earlier and eliminate obs that don't have a certain # of observations for that year
        oneyr <- oneobs[oneobs$origdate >= "2014-11-01" & oneobs$origdate <= "2020-03-31",]
        
        
        zooVals <- zoo(oneyr$V1, oneyr$origdate)
        zooapx <- as.data.frame(na.approx(zooVals))
        names(zooapx) <- 'zooapx'
        zooapx$date <- as.Date(rownames(zooapx))
        
        #add a row of NAs when date doesn't exist
        dates <- oneyr$origdate
        missing <- as.data.frame(as.Date(setdiff(dates, zooapx$date)))
        names(missing) <- 'date'
        if (nrow(missing) > 0){
          missing$zooapx <- "NA"
          zoofinal <- rbind.data.frame(missing, zooapx)
        } else {
          zoofinal <- zooapx
        }
        
        
        #create the final df
        #zoofinal <- rbind.data.frame(missing, zooapx)
        zoofinal <- merge(zoofinal, oneyr[,c(2,5)], by.x = "date", by.y="origdate" )
        zoonames <- paste("mean.X", zoofinal$date.y, sep="")
        zoofinal <- t(zoofinal[,2])
        colnames(zoofinal) <- zoonames
        nogaps <- cbind.data.frame(uid, zoofinal)
        
        gapfilled <- rbind(gapfilled, nogaps)
        
      }
  uid <- gapfilled$uid
  gapfilled <- lapply(gapfilled[,2:ncol(gapfilled)], as.numeric)
  gapfilled <- cbind.data.frame(uid, gapfilled)
  
  typgood <- merge(typgood[,c(1:8)], gapfilled, by="uid")    

  #calculate avg and sd
  end <- as.numeric(ncol(typgood))
  avgndvi  <- as.data.frame(colMeans(typgood[,c(9:end)], na.rm=TRUE))
  sdndvi <- as.data.frame(apply(typgood[,c(9:end)],2,sd, na.rm=TRUE))
  avgndvi <- cbind.data.frame(avgndvi, sdndvi)
  colnames(avgndvi) <- c('ndvi', 'sd')
  ##get dates
  avgndvi$date <- rownames(avgndvi)
  avgndvi$date <- substr(avgndvi$date, 7, 14)
  avgndvi$doy <- substr(avgndvi$date, 6, 8)
  avgndvi$year <- substr(avgndvi$date, 1, 4)
  #avgndvi$year <- 2014
  #avgndvi$year[as.numeric(avgndvi$doy) < 95] <- 2015
  
  avgndvi$origdate <- as.Date(as.numeric(avgndvi$doy), origin=as.Date(paste(avgndvi$year, "-01-01", sep="")))

  #remove rows with incomplete obs
  avgndvi <- avgndvi[complete.cases(avgndvi[ , 1:2]),]
  
  avgndvi$lower <- avgndvi$ndvi - avgndvi$sd
  avgndvi$upper <- avgndvi$ndvi + avgndvi$sd
  
  avgndvi$type <- cctype
  
  #make the plot
  #plot labels
  if(cctype=="ZConv"){
    lbl <- "ConventionalTillage"
  } else {
    lbl <- cctype
  }
  #lbl1 <- lbl
  lbl1 <- paste(lbl, " ", season, sep="")
  #nxtyear <- 2015
  nxtyear <- as.numeric(str_sub(season, -4, -1))+1
  xlabloc <- paste(nxtyear, "-03-01", sep="")
  #thisyr <- 2014
  thisyr <- as.numeric(str_sub(season, -4, -1))
 
  
  if(cctype=="Mix/Other"){
    outname <- "Mix"
  } else {
    outname <- cctype
  }
  
  #filename <- paste(outname, "_Allyrs.png", sep="")
  filename <- paste(outname, "_",  season, "_gapfilled", ".png", sep="")
  print(filename)  
  
  write.csv(avgndvi, paste("/N/project/CoverCrops/CoverCrops/data/NDVI/plotdata/", outname, "_",  season, "_gapfilled.csv", sep=""))
  
  
  #plot colors
  col <- cols.df[cols.df$alltypes == cctype, ]
  col <- col[1,2]
  
  #ribbon values
  eb <- aes(ymax = upper, ymin = lower)
  
  plot1 <- ggplot(data = avgndvi, aes( origdate, ndvi )) + 
    geom_line(color=col) +
    geom_ribbon(eb, fill=col, alpha=0.5)+
    ylim(0,0.9)+
    xlim(as.Date(paste(thisyr, "-11-01", sep="")), as.Date(paste(nxtyear, "-03-31", sep="")))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))  + #for no background or grid
    labs(y = "NDVI", x="Date") +
    annotate("text", label = lbl1, parse=FALSE, x = as.Date(xlabloc), y = 0.9, size = 3.5)
    
    
  plot1
  ggsave(filename, device='png', width=10, height=7, plot=plot1, dpi = 300, units="cm")
}

setwd("/N/project/CoverCrops/CoverCrops/data/NDVI/plots/")
alltypes <- unique(plotdata$cc.agg)
lapply(alltypes, plotndvi, season='Fall2015')


######plotting wheat, ryegrass, no till, and conventional till on the same graph
files <- list.files("/N/project/CoverCrops/CoverCrops/data/NDVI/plotdata/", full.names = TRUE)
plotdata <- do.call(rbind, lapply(files, read.csv))

summary(as.factor(plotdata$type))
plotdata <- subset(plotdata, type %in% c("Wheat", "RyeGrass", "NoTill", "ZConv"))


eb <- aes(ymax = upper, ymin = lower)


plot1 <- ggplot(data = plotdata, aes(x = as.Date(origdate), y = ndvi, color=type, fill=type)) +
  geom_line() +
  scale_colour_manual(values=c("#69585F","#DFEFCA", "#FFF9A5", "#60656F")) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill=type), alpha=0.5)+
  scale_fill_manual(values=c("#69585F","#DFEFCA", "#FFF9A5", "#60656F"))+
  ylim(0,0.9)+
  xlim(as.Date(paste("2015", "-11-01", sep="")), as.Date(paste("2016", "-03-31", sep="")))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))  + #for no background or grid
  labs(y = "NDVI", x="Date")

plot1