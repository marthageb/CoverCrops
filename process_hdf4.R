###convert HLS .hdf4 files into .tif files for downstream processing (single .tif per band x date) 
###**if processing on RED need to start R with ‘launchR2.sh’** to ensure correct hdf4 drivers are installed

library(rgdal)
library(raster)
library(gdalUtils)
#install.packages('unixtools',repos='http://www.rforge.net/')
unixtools::set.tempdir("/N/slate/farellam")

################################################################################################################################################
##############################convert .hdf files to .tif files##################################################################################
################################################################################################################################################
##takes about 15 minutes for a single title.

convert_hdf4 <- function(tile, year){
  workdir <- "/N/project/CoverCrops/CoverCrops/HLS_L30_v1.4"
  hdfs1 = list.files(paste(workdir, tile, year, sep="/"), pattern="hdf$")
  setwd(paste(workdir, tile, year, sep="/"))
  
  
  outdir <- paste(workdir, "/", tile, "_", year, "_Tiffs/", sep="")
  dir.create(outdir)
  gtiffs1 = gsub("hdf","tif",hdfs1) #
  fromSRS = "+proj=utm +zone=16 +ellps=WGS84 + datum=WGS84 + units=m + no_defs" # original HDF SRS
  toSRS = "+proj=longlat +datum=WGS84 +no_defs" # desired GeoTIFF SRS
  
  beg <- Sys.time()
  for(i in 1:length(hdfs1)){
    gc()
    outname <- paste(outdir, "Band1_", gtiffs1[i], sep="")
    gdal_translate(hdfs1[i], outname, sd_index=1)
    #gdalwarp(outname, outname, s_srs=fromSRS, t_srs=toSRS, srcnodata=-1000, dstnodata=NA, overwrite = T) # project geotiffs
    print(paste(substr(hdfs1[i], 9, 22), ': done with B1', sep=""))
    
    outname <- paste(outdir, "Band2_", gtiffs1[i], sep="")
    gdal_translate(hdfs1[i], outname, sd_index=2)
    #gdalwarp(outname, outname, s_srs=fromSRS, t_srs=toSRS, srcnodata=-1000, dstnodata=NA, overwrite = T) # project geotiffs
    print(paste(substr(hdfs1[i], 9, 22), ': done with B2', sep=""))
    
    outname <- paste(outdir, "Band3_",gtiffs1[i], sep="")
    gdal_translate(hdfs1[i], outname, sd_index=3)
    #gdalwarp(outname, outname, s_srs=fromSRS, t_srs=toSRS, srcnodata=-1000, dstnodata=NA, overwrite = T) # project geotiffs
    print(paste(substr(hdfs1[i], 9, 22), ': done with B3', sep=""))
    
    outname <- paste(outdir, "Band4_",gtiffs1[i], sep="")
    gdal_translate(hdfs1[i], outname, sd_index=4)
    #gdalwarp(outname, outname, s_srs=fromSRS, t_srs=toSRS, srcnodata=-1000, dstnodata=NA, overwrite = T) # project geotiffs
    print(paste(substr(hdfs1[i], 9, 22), ': done with B4', sep=""))
    
    outname <- paste(outdir, "Band5_",gtiffs1[i], sep="")
    gdal_translate(hdfs1[i], outname, sd_index=5)
    #gdalwarp(outname, outname, s_srs=fromSRS, t_srs=toSRS, srcnodata=-1000, dstnodata=NA, overwrite = T) # project geotiffs
    print(paste(substr(hdfs1[i], 9, 22), ': done with B5', sep=""))
    
    outname <- paste(outdir, "Band6_",gtiffs1[i], sep="")
    gdal_translate(hdfs1[i], outname, sd_index=6)
    #gdalwarp(outname, outname, s_srs=fromSRS, t_srs=toSRS, srcnodata=-1000, dstnodata=NA, overwrite = T) # project geotiffs
    print(paste(substr(hdfs1[i], 9, 22), ': done with B6', sep=""))
    
    outname <- paste(outdir, "Band7_",gtiffs1[i], sep="")
    gdal_translate(hdfs1[i], outname, sd_index=7)
    #gdalwarp(outname, outname, s_srs=fromSRS, t_srs=toSRS, srcnodata=-1000, dstnodata=NA, overwrite = T) # project geotiffs
    print(paste(substr(hdfs1[i], 9, 22), ': done with B7', sep=""))
    
    outname <- paste(outdir, "Band9_",gtiffs1[i], sep="")
    gdal_translate(hdfs1[i], outname, sd_index=8)
    #gdalwarp(outname, outname, s_srs=fromSRS, t_srs=toSRS, srcnodata=-1000, dstnodata=NA, overwrite = T) # project geotiffs
    print(paste(substr(hdfs1[i], 9, 22), ': done with B9', sep=""))
    
    outname <- paste(outdir, "Band10_", gtiffs1[i], sep="")
    gdal_translate(hdfs1[i], outname, sd_index=9)
    #gdalwarp(outname, outname, s_srs=fromSRS, t_srs=toSRS, srcnodata=-1000, dstnodata=NA, overwrite = T) # project geotiffs
    print(paste(substr(hdfs1[i], 9, 22), ': done with B10', sep=""))
    
    outname <- paste(outdir, "Band11_",gtiffs1[i], sep="")
    gdal_translate(hdfs1[i], outname, sd_index=10)
    #gdalwarp(outname, outname, s_srs=fromSRS, t_srs=toSRS, srcnodata=-1000, dstnodata=NA, overwrite = T) # project geotiffs
    print(paste(substr(hdfs1[i], 9, 22), ': done with B11', sep=""))
    
    outname <- paste(outdir, "QA_",gtiffs1[i], sep="")
    gdal_translate(hdfs1[i], outname, sd_index=11)
    #gdalwarp(outname, outname, s_srs=fromSRS, t_srs=toSRS, srcnodata=-1000, dstnodata=NA, overwrite = T) # project geotiffs
    print(paste(substr(hdfs1[i], 9, 22), ': done QA', sep=""))
  }
  end <- Sys.time()
  diff <- end - beg
  print(paste('time to process ', length(hdfs1), " files = ", diff, sep=""))
}


#yr = 2020
convert_hdf4(tile = '16SDH', year = 2019)
