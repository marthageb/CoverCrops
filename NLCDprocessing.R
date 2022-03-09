####pre-process the NLCD dataset for downstream analysis: reproject to HLS crs and crop to extent covered by HLS tiles.


library(raster)
library(rgdal)

#Landcover_Rast <- raster("/Users/mallory/Documents/Temp_Project/landcvi020l_nt00016/landcover_proj.tif")
rasterOptions(tmptime = 24,progress="text",timer=TRUE,overwrite = T,chunksize=2e+08,maxmemory=1e+8)
Landcover_Rast <- raster("/N/project/CoverCrops/CoverCrops/data/spatial/NLCD/NLCD_2016/nlcd_2016_land_cover_l48_20210604.img")
plot(Landcover_Rast)
dataType(Landcover_Rast)="INT4S"
extent(Landcover_Rast)
extent(e2)
#e2 <- extent(500000, 2300000, 177285, 2900000)
e2 <- extent(600000, 975000, 1600000, 2200000)

cropped <- crop(Landcover_Rast, e2)
plot(cropped)

#Trying to project 
#load hls raster (the raster you want this CRS to match)
hls1 <- stack("/N/project/CoverCrops/CoverCrops/data/spatial/RFstacks/HLS_stacks/cloudmsk/RFstack_months_16TDL_2014.grd")
Landcover <- projectRaster(cropped, crs=crs(hls1))
plot(Landcover)

#crop landcover raster to extent of HLS tiles
ext <- extent(399960, 709800, 4090200, 4700040)
LC_crop <- crop(Landcover, ext)
writeRaster(LC_crop, "/N/project/CoverCrops/CoverCrops/data/spatial/NLCD/NLCD2016_HLScrop.tif")
