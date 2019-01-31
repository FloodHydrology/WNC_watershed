####################################################################################
# Name: WNC Threshold-Based Watershed Delineation
# Coder: C. Nathan Jones
# Date: 20 Jan 2019
# Purpose: Create function to delineate watershed based on area thresholds
####################################################################################

####################################################################################
# Step 1: Setup Worskspace ---------------------------------------------------------
####################################################################################
# Clear Memory  
rm(list=ls(all=TRUE))

# Load packages 
library(sf)         # for spatial analysis
library(raster)     # for spatial analysis
library(tidyverse)  # for data wrangling

#define relevant directories
wbt_dir<-     "C:\\WBT/whitebox_tools"
scratch_dir<- "C:\\ScratchWorkspace/"
data_dir<-    "//storage.research.sesync.org/njones-data/Research Projects/WNC_watershed/spatial_data/RawData/"
output_dir<-   "//storage.research.sesync.org/njones-data/Research Projects/WNC_watershed/spatial_data/DerivedData/"

#Download data (these are from 1_basic_watershed_delineation.R)
dem<-raster(paste0(data_dir,"WS2_DEM1.tif"))
watershed<-raster(paste0(output_dir, "watershed.tif"))

#Crop DEM to watershed
dem<-dem*watershed
dem<-crop(dem, extent(c(1353165,1357325,961094,967555.7)))

#Plot
plot(dem)

####################################################################################
# Step 2: DEM Preprocessing ---------------------------------------------------------
####################################################################################
#Export DEM and stream layer to local working directory
writeRaster(dem, 
            paste0(scratch_dir,"dem.tif"), 
            overwrite=T)

#Gaussian Filter
system(paste(paste(wbt_dir), 
             "-r=GaussianFilter", 
             paste0("--wd=",scratch_dir),
             "-i='dem.tif'", 
             "-o='dem_filter.tif'",
             "--sigma=3"))

#Fill "single cell" depressions
system(paste(paste(wbt_dir),
             "-r=FillSingleCellPits",
             paste0("--wd=",scratch_dir),
             "--dem='dem_filter.tif'",
             "-o='dem_breach_minor.tif'"))

#Breach larger depressions
system(paste(paste(wbt_dir), 
             "-r=BreachDepressions", 
             paste0("--wd=",scratch_dir),
             "--dem='dem_breach_minor.tif'", 
             "-o='dem_breach_major.tif'"))

#Conduct flow accumulation and flow direction analysis
system(paste(paste(wbt_dir), 
             "-r=D8FlowAccumulation", 
             "--out_type='cells'",
             paste0("--wd=",scratch_dir),
             "--dem='dem_breach_major.tif'", 
             "-o='fac.tif'"))
system(paste(paste(wbt_dir), 
             "-r=D8Pointer", 
             paste0("--wd=",scratch_dir),
             "--dem='dem_breach_major.tif'", 
             "-o='fdr.tif'",
             "--out_type=sca"))

#Retreive fac and fdr rasters for subsequent steps
fac<-raster(paste0(scratch_dir,"fac.tif"))
  fac@crs<-dem@crs
fdr<-raster(paste0(scratch_dir,"fdr.tif"))
  fdr@crs<-dem@crs

####################################################################################
# Step 3: Creat delineation function------------------------------------------------
####################################################################################
#Create function to 
threshold_fun<-function(fac,         #flow accumulation raster
                        fdr,         #flow direction raster
                        threshold,   #threshold for watersheds [in map units]
                        scratch_dir, #scratch workspace
                        output_dir,  #where to store the output rasters!
                        output_name){
  
  #Estimate number of cells for threshold
  threshold<-threshold/res(fac)[1]/res(fac)[2]
  
  #use threshold to define the flow_net
  flowgrid<-fac
  flowgrid[flowgrid<threshold]<-NA
  flowgrid@crs<-dem@crs
  
  #Create function to find local min
  min_fun<-function(x, ...){
    #define center cell value
    center<-x[ceiling(length(x)/2)]
    
    #define min value in window
    min<-min(x, na.rm=T)
    
    #If the center cell is the local minumum, return "1"
    output<-ifelse(min==center, 1,0)
    
    #define output cell value
    output
  }
  
  #use function to find local minima in raster with 3x3 moving window
  outlet_grd<-focal(x=flowgrid, w=matrix(1,3,3), fun=min_fun, na.rm=T)

  #convert outlet_grd to points
  outlet<-data.frame(rasterToPoints(outlet_grd)) %>%
    filter(layer==1)
  outlet<-st_as_sf(outlet, coords=c("x","y"), crs=paste(fac@crs))

  #Write fdr and points to scratch dir
  writeRaster(fdr,
              paste0(scratch_dir, "fdr_fun.tif"),
              overwrite=T)
  st_write(outlet, paste0(scratch_dir,"outlet_fun.shp"), delete_layer=T)
  
  #Delineate the watersheds
  system(paste(paste(wbt_dir),
               "-r=Watershed", 
               paste0("--wd=",scratch_dir),
               "--d8_pntr='fdr_fun.tif'", 
               "--pour_pts='outlet_fun.shp'",
               "-o='watershed.tif"))
  
  #load watershed shape into memory
  sheds<-raster(paste0(scratch_dir,"watershed.tif"))
  writeRaster(sheds,paste0(output_dir,output_name))

  #export watershed grd (you have to reread b/c raster does not store data in memory)
  sheds<-raster(paste0(output_dir,output_name))
  sheds
}

#Estimate watersheds that are 1%,5%,and 10% the size of the original watershed
#thresholds
one<-cellStats(watershed, sum)*0.01*3.25^2
five<-cellStats(watershed, sum)*0.05*3.25^2
ten<-cellStats(watershed, sum)*0.10*3.25^2
#delineation
sheds_1pcnt <-threshold_fun(fac, fdr, one,  scratch_dir, output_dir, "test_01pct.tif")
sheds_5pcnt <-threshold_fun(fac, fdr, five, scratch_dir, output_dir, "test_05pct.tif")
sheds_10pcnt<-threshold_fun(fac, fdr, ten,  scratch_dir, output_dir, "test_10pct.tif")