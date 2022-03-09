###Move observations closer to road. Locations chosen by REU students did not fall in desired fields-- too far from road points. 
###Calculate new coordinates, create .shp file for all survey location, and make 25m buffer .shp file around new points


library(rgdal)
library(stringr)
library(sp)
library(tidyverse)
library(sf)
library(dplyr)

###########Step 1: get the coordinates and Field IDs for the center of the 120m buffers###########
#load the current coordinates I'm using (the center of the 120m buffers)
#data <- read.csv("/N/project/CoverCrops/CoverCrops/data/RFdata/HLS/cloudmsk_alltiles.csv")
data <- read.csv("/N/project/CoverCrops/CoverCrops/data/surveydata_RandL.csv")
#select a single site
site <- subset(data, site=="Whitley")
#select field, long, and lat columns 
field <- site[,c(6:8)]
#create column with L or R fields 
field$RL <- str_sub(field$Field_No, -1, -1)

#subset L and R fields, then create 'ident' column
Rfield <- subset(field, RL=="R")
Lfield <- subset(field, RL=="L")
Rfield$ident <- as.numeric(strsplit(Rfield$Field_No, "R"))
Lfield$ident <- as.numeric(strsplit(Lfield$Field_No, "L"))

#remove duplicated observations
Rfield <- unique(Rfield)
Lfield <- unique(Lfield)

#remove unwanted coulmns
Rfield <- Rfield[,c(1,2,5)]
Lfield <- Lfield[,c(1,2,5)]
colnames(Rfield) <- c("buffY", "buffX", "ident")
colnames(Lfield) <- c("buffY", "buffX", "ident")

Rfield$ident <- paste(Rfield$ident, "R", sep="")
Lfield$ident <- paste(Lfield$ident, "L", sep="")


###########Step 2: get the coordinates and Field IDs for the transect coords (on the road)###########
#navigate to the folder mallory shared via google drive with the transect and REU data
setwd("/N/project/CoverCrops/CoverCrops/data/IN cover crop transect data/")


#load the road transect points
#Benton County: readOGR("Buffer and transect data joins/Benton/Benton_Pts.shp")
#Decatur County: readOGR("Buffer and transect data joins/Decatur/1_Decatur_points.shp")
#Gibson County: readOGR("Buffer and transect data joins/Gibson/Gibson_Pts.shp")
#Hamilton: readOGR("Buffer and transect data joins/Hamilton/1_Hamilton_points.shp")
#Posey: readOGR("Buffer and transect data joins/Posey/Posey_Pts.shp")
#St.Joseph: readOGR("Buffer and transect data joins/StJoseph/1_StJoseph_points.shp")
#Warren: readOGR("Buffer and transect data joins/Warren/Warren_Pts.shp")
#Washington: readOGR("Buffer and transect data joins/Washington/1_Washington_points.shp")
#White: readOGR("Buffer and transect data joins/White/White_Pts.shp")
#Whitley: readOGR("Buffer and transect data joins/Whitley/1_Whitley_points.shp")
road <- readOGR("Buffer and transect data joins/Whitley/1_Whitley_points.shp")
road.df <- as.data.frame(road)
#select ident, lat, and long columns
road.df <- road.df[,c(2, 11, 10)]
#Change names to match R and L field names
colnames(road.df) 
colnames(road.df) <- c("ident", "roadY", "roadX")
##SKIP TO LINE 77 if L and R are already included on 'ident'

Lroad <- road.df
Lroad$ident <- paste(Lroad$ident, "L", sep="")

Rroad <- road.df
Rroad$ident <- paste(Rroad$ident, "R", sep="")


###########Step 3: combine the data frames and calculate coords 30m from road points###########
Rdf <- merge(Rroad, Rfield, by="ident")
Ldf <- merge(Lroad, Lfield, by="ident")
all.df <- rbind.data.frame(Rdf, Ldf)


###if R and L are already on 'ident'##
fields <- rbind.data.frame(Rfield, Lfield)
all.df <- merge(road.df, fields, by="ident")


#to find the point that is 30m away from the road: https://stackoverflow.com/questions/1934210/finding-a-point-on-a-line
#distance between the start and end points (should be ~ 0.00181)
d = sqrt((all.df$buffX - all.df$roadX)^2 + (all.df$buffY - all.df$roadY)^2)

#segment ratio-- we want the point that is apx 30m (0.0003 decimal degrees) away from point 1
r = 0.0003/d

#find the point that divides the segment into the ratio (1-r):r
all.df$x3 <- r * all.df$buffX + (1 - r) * all.df$roadX
all.df$y3 <- r * all.df$buffY + (1 - r) * all.df$roadY

final <- all.df[,c(1,6,7)]
colnames(final) <- c('ident', 'newX', 'newY')

#make sure the number of fields match up
length(unique(site$Field_No))
nrow(final)

#write the new coords
write.csv(final, "/N/project/CoverCrops/CoverCrops/data/newcoords_Whitley.csv")

####################################################################################################################
########################################CREATE A SHAPEFILE FOR ALL SURVEY LOCS######################################
####################################################################################################################
#list the files with the new coords for all counties
files <- list.files("/N/project/CoverCrops/CoverCrops/data/", pattern="newcoords")

setwd("/N/project/CoverCrops/CoverCrops/data/")
allcords <- function(x){
  county <- unlist(str_split(x, "_"))[2]
  county <- unlist(str_split(county, ".csv"))[1]
  coords <- read.csv(x)
  coords <- coords[,-1]
  coords$county <- county
  return(coords)
}

all_cos <- lapply(files, allcords) %>% bind_rows()

#get crs 
crs.df <- readOGR("/N/project/CoverCrops/CoverCrops/data/IN cover crop transect data/Buffer and transect data joins/Benton/Benton_Pts.shp")

mycrs <- crs(crs.df)
final_shp <- SpatialPointsDataFrame(all_cos[,c(2,3)], all_cos, proj4string = mycrs)

names(final_shp) <- c("Field_No", "Longitude", "Latitude", "site")
writeOGR(final_shp, "/N/project/CoverCrops/CoverCrops/data/spatial/sites.shp", layer="sites", driver = "ESRI Shapefile", overwrite_layer = TRUE)


####################################################################################################################
##################################CREATE a 25M BUFFER AROUND SURVEY LOC POINTS######################################
####################################################################################################################
#load survey loc .shp file
locs <- readOGR("/N/project/CoverCrops/CoverCrops/data/spatial/sites.shp")
#create UID
locs$uid <- paste(locs$site, locs$Field_No, sep="")

#get in a projected coordinate system that has meters as the units
locs_t <- spTransform(locs, CRS('+init=epsg:4326'))
locs_p <- spTransform(locs_t, CRS("+proj=utm +zone=16 +datum=WGS84"))
#convert to sf object
locs_sf <- st_as_sf(locs_p)

#create buffers
buffs <- st_buffer(locs_sf, dist=25)

#convert back to a spatial object
buffs.polys <- as(buffs, "Spatial")

#write the buffer polygons
writeOGR(buffs.polys, "/N/project/CoverCrops/CoverCrops/data/spatial/25mbuffs.shp", layer='data', driver="ESRI Shapefile")

