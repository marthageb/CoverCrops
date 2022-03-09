####Apply the 2- and 3-class RF models at the county-level. 
    ##For each county, combine HLS stack with LSTsd and LSTmed (need to transform CRS and resample LST rasters). Mosaic HLS and LST tiles when more than one tile cover transect locations.
    ##Load the RF2class and RF3class models and then apply to RFstacks. Save 2- and 3-class predictions



library(raster)
options("rgdal_show_exportToProj4_warnings"="none")
library(rgdal)
library(doParallel)

#define the HLS and LST tiles
HLStile <- c("16SEH", "16SEJ")
LSTtile <- '022010'
year <- "2015" #define the year
county <- 'washington' #define the county

####GET LST raster to 'fit' the other raster layers####
#load the HLS raster stack
hls1 <- stack(paste("/N/project/CoverCrops/CoverCrops/data/spatial/RFstacks/HLS_stacks/cloudmsk/RFstack_months_", HLStile[1], "_", year, ".grd", sep=""))
hls2 <- stack(paste("/N/project/CoverCrops/CoverCrops/data/spatial/RFstacks/HLS_stacks/cloudmsk/RFstack_months_", HLStile[2], "_", year, ".grd", sep=""))
hls.all <- list(hls1, hls2)
names(hls.all) <- NULL
hls.all$fun <- mean
hls <- do.call(mosaic,hls.all)


#load the LST raster, for some reason 'projectRaster' in R didn't reproject the raster to the correct location so instead, I used ArcPro to reproject the LST rasters. ArcPro project saved: C:/Users/farellam/Documents/ArcGIS/Projects/reProjLST
#ArcPro reprojected rasters saved: slate/farellam/CoverCrops/data/spatial/seasLSTYEAR_ArcReProj.tif
lstmed <- raster(paste("/N/project/CoverCrops/CoverCrops/HLS30_Indiana/all_fall/LSTmedfall", year, "_", LSTtile, ".tif", sep=""))
lstsd <- raster(paste("/N/project/CoverCrops/CoverCrops/HLS30_Indiana/LSTsdfall", year, "_", LSTtile, ".tif", sep=""))
lst <- stack(lstsd, lstmed)
#This is to fix the crs issue I created in 'processLST.R'. Need to load the original LST raster, redefine the projection to the wonky projection. Then reassign the crs of the median and SD LST raster this wonky CRS. And can THEN transform median and SD LST to the HLS projection.
oldwd <- "/N/project/CoverCrops/CoverCrops/BulkOrderIndianaProvisionalSurfaceTemp/U.S. Landsat 4-8 ARD/"
ogfiles <- list.files(oldwd, recursive = TRUE, pattern = "*_ST.tif", full.names = FALSE)
onetile <- ogfiles[grep(LSTtile, ogfiles)]
oneyear <- onetile[grep(year, onetile)]
one <- oneyear[[1]]
lstog <- raster(paste(oldwd, one, sep=""))
badcrs <- "+proj=longlat +datum=WGS84 +no_defs +datum=WGS84"
mybad <- projectRaster(lstog, crs=badcrs)
#now, redefine the lst CRS and transform to HLS CRS
crs(lst) <- crs(mybad)
lst_p <- projectRaster(lst, crs=crs(hls), method='bilinear')
print("lst and hls raster stacks have the same projection")

#load the county polygon
copoly <- readOGR(paste("/N/project/CoverCrops/CoverCrops/data/spatial/INcounties/indv_cos/", county, ".shp", sep=""))
copoly_t <- spTransform(copoly, crs(hls)) #transform to HLS CRS

#clip the HLS and LST data to the county area
cohls <- raster::crop(hls, copoly_t)
colst <- raster::crop(lst_p, copoly_t)


#resample the lst raster so that it is on the same grid as the raster stack -- takes about 1 min for a tile on 12 cores
beginCluster(12)
lst.resamp <- resample(colst, cohls)

#stack the two raster stacks together
rfstack <- stack(lst.resamp, cohls)

#write the raster stack
writeRaster(rfstack, paste("/N/project/CoverCrops/CoverCrops/data/spatial/RFdata/RFstack_" ,year, county, ".grd", sep=""), format="raster")

print(paste("RFstack written for ", county, " ", year, sep=""))

names(rfstack) <- c("sdLST", "medLST", "medB3", "minB3", "maxB3", "medB5", "minB5", "maxB5", "medB6", "minB6", "maxB6", "medB10", "minB10", "maxB10", "medB11", "minB11", "maxB11",
                    "medNDVI", "minNDVI", "maxNDVI", "medNDTI", "minNDTI", "maxNDTI", "medSTI", "minSTI", "maxSTI", "medExG", "minExG", "maxExG", "medExGR", "minExGR", "maxExGR",
                    "full_maxNDVI", "ampNDVI", "ratioNDVI")

#########################################################################################
##########################APPLY THE 2 and 3 CLASS MODEL###################################
#########################################################################################
#load the RFmodel
RF2class <- readRDS("/N/project/CoverCrops/CoverCrops/RFmodels/RF2class_Conv_CC.rds")
RF3class <- readRDS("/N/project/CoverCrops/CoverCrops/RFmodels/RF3class_HLS_cloudmsk_25mbuff.rds")

#apply the RFmodel  
#trying two diff parallel processing fxns. OPTION 1:
#beginCluster(5)
#twoclass <- clusterR(RFstack, fun = predict, args = list(model = RF2class))
#endCluster()
start <- Sys.time()
#OPTION 2:
cl <- makePSOCKcluster(12)
registerDoParallel(cl)
twoclass <- predict(rfstack, RF2class)
threeclass <- predict(rfstack, RF3class)
stopCluster(cl)
end <- Sys.time()
print(paste("took ", end-start, " to run the 2 and 3 class models on ", county, " county", sep=""))


plot(threeclass)


#writeRaster(twoclass, paste("/N/slate/farellam/CoverCrops/data/spatial/RFresults/3class/", seas, "_", pathrow,".tif", sep=""), format="GTiff", overwrite=TRUE)
outdir <- "/N/project/CoverCrops/CoverCrops/data/spatial/RFresults/"
raster::writeRaster(twoclass, paste(outdir, "2class/", county, year, ".tif", sep=""), format="GTiff", overwrite=TRUE)  
raster::writeRaster(threeclass, paste(outdir, "3class/", county, year, ".tif", sep=""), format="GTiff", overwrite=TRUE)  



freq2 <- as.data.frame(freq(twoclass))
freq3 <- as.data.frame(freq(threeclass))

levels2 <- rbind(as.data.frame(RF2class$levels), "NA")
levels3 <- rbind(as.data.frame(RF3class$levels), "NA")

final2 <- cbind.data.frame(levels2, freq2)
final2 <- final2[,-2]
names(final2) <- c("class", "count")
final2

final3 <- cbind.data.frame(levels3, freq3)
final3 <- final3[,-2]
names(final3) <- c("class", "count")
final3
