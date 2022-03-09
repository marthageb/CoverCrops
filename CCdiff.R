####Make .shp files with CC acreage differences from ISDA observed acreage in the 2- and 3-class models. 


library(rgdal)
library(raster)
library(extrafont)
library(ggplot2)
library(lattice)
require(gridExtra)

setwd("C:/Users/farellam/OneDrive - Indiana University/Documents/CoverCrops/")
incos <- readOGR("data/spatial/INcounties.shp")
plot(incos)

#subset just the counties included in my analysis
mycos <- subset(incos, NAME %in% c("Benton", "Decatur", "Gibson", "Hamilton", "Posey", "St. Joseph", "Warren", "Washington", "White", "Whitley"))

#load the CC predictions
pred <- read.csv('/N/project/CoverCrops/CoverCrops/data/CCpred_yearxco.csv')
pred <- pred[,-1]


onepred <- subset(pred, prediction=="2class")
rownames(onepred) <- onepred$year
onepred <- onepred[,-c(1:2)]
names(onepred)<-str_to_title(names(onepred))
names(onepred)[names(onepred) == 'Stjoe'] <- "St. Joseph"


#onepred_t <- as.data.frame(t(onepred))
#names(oneyr_t)[1] <- "twoclass"
#names(oneyr_t)[2] <- "threeclass"
#onepred_t$NAME <- rownames(onepred_t)


isda <- subset(pred, prediction=="ISDA")
rownames(isda) <- isda$year
isda <- isda[,-c(1:2)]
names(isda)<-str_to_title(names(isda))
names(isda)[names(isda) == 'Stjoe'] <- "St. Joseph"
#isda_t <- as.data.frame(t(isda))
#names(oneyr_t)[1] <- "twoclass"
#names(oneyr_t)[2] <- "threeclass"
#isda_t$NAME <- rownames(isda_t)

final <- onepred-isda
avg <- colMeans(final)
min <- apply(final,2,min)
max <- apply(final,2,max)

final2 <- cbind.data.frame(avg, min, max)
final2$NAME <- rownames(final2)

#make a data frame with NLCD cropland acres for each county
NAME <- c("Benton", "Decatur", "Gibson", "Hamilton", "Posey", "St. Joseph", "Warren", "Washington", "White", "Whitley")
acres <- c(235291.5, 159391.1, 222650.4, 131887.9, 184701.6, 155122.0, 173500.4, 133807.0, 268745.3, 154038.7)
acres.df <- cbind.data.frame(NAME, acres)

final2 <- merge(final2, acres.df, by="NAME")
final2$per <- final2$avg/final2$acres * 100


cccos <- merge(mycos, final2, by="NAME")
cc.df <- as.data.frame(cccos)
#writeOGR(cccos, "C:/Users/farellam/OneDrive - Indiana University/Documents/CoverCrops/data/spatial/cccos2class.shp", driver="ESRI Shapefile", layer="cccos")
