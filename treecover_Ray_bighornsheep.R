#Code written by Cynthia L. Norton, Univserity of Arizona, 2024 
#Detect tree cover more than 1m

setwd("//snow/projects/RaBET_CLN/RaBET_upscale/Utah")
options(scipen = 100, digits = 4)

###DOWNLOADING 3DEP LAZ FILES
# Path to the text file with URLs
urls_file <- "//snow/projects/RaBET_CLN/RaBET_upscale/Utah/MOABsite/data/3DEP.txt"

# Directory to save files
save_dir <- "//snow/projects/RaBET_CLN/RaBET_upscale/Utah/MOABsite/data/LiDAR"

# Read URLs from the text file
urls <- readLines(urls_file)
options(timeout=500)

# Function to download files
download_laz <- function(url, dest_dir) {
  file_name <- basename(url)
  dest_file <- file.path(dest_dir, file_name)
  download.file(url, dest_file, mode = "wb")
}

# Download all files bulk
for (url in urls) {
  download_laz(url, save_dir)
}

# Optional: Check the files
list.files(save_dir)

warnings()


folder = list.files('//snow/projects/RaBET_CLN/RaBET_upscale/Utah/ONAQsite/data/LiDAR/',pattern="*.laz$", full.names=TRUE) #change file path


#baselayer points to rasters
#Canopy height model
process_lidar_file <- function(file) {
  # Read the LAS file
  las <- lidR::readLAS("//snow/projects/RaBET_CLN/RaBET_upscale/Utah/ONAQsite/data/LiDAR/USGS_LPC_UT_Central_QL2_2018_12TUK6644_LAS_2019.laz")
  
  las_classify <- classify_ground(las, algorithm = csf(sloop_smooth = TRUE, class_threshold = 1, cloth_resolution = 1, rigidness = 2))
  
  # Filter ground and tree points
  las_ground <- filter_poi(las_classify, Classification == 2)
  las_tree <- filter_poi(las_classify, Classification != 2)
  
  # Generate DTM
  las_DTM <- grid_terrain(las_ground, res = 1, algorithm = knnidw(k = 10, p = 2))
  
  # Normalize height
  las_normalized_height <- normalize_height(las_classify, algorithm = las_DTM)
  
  # Filter AGL points
  las_AGL_clean <- filter_poi(las_normalized_height, Z >= 0.1, Z <= 40)
  
  # Create CHM
  las_CHM <- rasterize_canopy(las_normalized_height, res = 1)
  
  # Locate tree tops
  las_ttops <- locate_trees(las_CHM, lmf(5))
  
  # Segment trees
  las_segmented <- segment_trees(las_AGL_clean, dalponte2016(las_CHM, las_ttops))
  
  terra::crs(las_CHM) <-  simplified_crs

  # Export metrics as TIFF files
  file_base <- tools::file_path_sans_ext(basename(file))

  terra::writeRaster(las_CHM, file.path(output_folder, paste0("CHM_date",file_base, '.tif')), overwrite = TRUE) ##change date
  #slope
  #hillshade
  #aspect
}

lapply(folder, process_lidar_file) #change folder to list of files



##MOSAIC 


##subset filter rasters based on criterias wanted
#filter CHM

##CHM and species###
CHM[CHM <= 1] = NA #removing canopy less than or equal to 1m
CHM_project<- terra::project(CHM,hillshade) #projecting so in the same CRS
CHM_hillshade <- mask(hillshade, CHM) #masking hillshade based on CHM
writerast(CHM_hillshade, file = "CHM_masked_hillshade_112124.tif") #export



#