# Test of Flight

library(rhdf5)
library(raster)
library(lubridate)
library(doParallel)
library(tictoc)
library(lfe)
library(dplyr)
library(caTools)
library(tictoc)

tic()
filepath <- "20220327t204552_SanBernardinoCA_L1_B102_V01.hdf5"
hdr_file <- "20220327t204552_SanBernardinoCA_L1_B102_V01.geo.hdr"

setwd("/Users/Ryan/Documents/Super Resolution/") 

# Gives structure of the hdf5 file
h5ls(filepath)
# Pixel_geolocation: 4 x 512 x 3307
# First layer: Latitude
# Second layer: Longitude
# Third layer: Height (measurement unknown)
# Fourth layer: # of steps taken during ray-casting?

# Read hdf5 file
pixel <- h5read(filepath, 'pixel_geolocation')
# pixel is a large array of 3-dimensions (4 x 512 x 3307)

# Attempt to load hdr file
file.exists(hdr_file)
# read.ENVI for some reason adds ".hdr" except when I need it to work haha
# test <- read.ENVI(hdr_file)

# Save as a data frame
pixel_df_old <- as.data.frame(pixel)
# 4 x 1693184

# Transposes the data frame
pixel_df_new <- as.data.frame(t(pixel_df_old))
head(pixel_df_new)

# Adds column names
pixel_df_new <- setNames(pixel_df_new, c('Latitude','Longitude','Height','Steps'))
head(pixel_df_new)

write.csv(pixel_df_new, "SB_test_location.csv", row.names = FALSE)
toc()
# Produced a 97 MB file
# Too large to open in file preview (Max: 5 MB)
# Processing time took 21.93 sec according to tictoc. Most of the computational time was due to transposing data and writing it out to csv file.