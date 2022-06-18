# This file loads the LST dataframe!

library(rhdf5)
library(raster)
library(lubridate)
library(doParallel)
library(tictoc)
library(lfe)
library(dplyr)
library(caTools)
library(reshape)

tic()
LST_filepath <- "20220327t204552_SanBernardinoCA_L2_B102_V01.hdf5"

setwd("/Users/Ryan/Documents/Super Resolution/") 

# Gives structure of the hdf5 file
h5ls(LST_filepath)
# L2_LST: 512 x 3307

# Read hdf5 file
LST_pixel <- h5read(LST_filepath, 'L2_LST')
# pixel is a large array of 2-dimensions (512 x 3307)

# Save as a data frame
LST_pixel_df_old <- as.data.frame(LST_pixel)
# 512 x 3307


# Melts the data frame into one column array
LST_pixel_df_new <- LST_pixel_df_old %>% melt() %>% dplyr::rename(LST = value) %>% select(-1)
head(LST_pixel_df_new)
toc()