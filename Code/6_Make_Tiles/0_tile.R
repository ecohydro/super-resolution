#####################
# Tile images of an arbitrary size into equal size tiles. 
# Equal margins are taken on each side if necessary. 
# Anna Boser Aug 2 2022
#####################

# load packages -----------------------------------------------------------

library(tidyverse)
library(raster)
library(here)
library(dplyr)
library(parallel)
library(stringr)
library(tictoc)

# parameters and dirs --------------------------------------------------------------

#input directories
input_dir <- here("Data", "Final")
variables <- c("LST", "RGB")

# output file
output_dir = here("Data", "Tiles")

# list all the files
filenames <- list.files(here(input_dir, variables[1])) # same in rgb and lst

# if this script has already been run and some runs are already processed: 
if (file.exists(here(output_dir, "done_files.txt"))){
  done_files <- read.csv(here(output_dir, "done_files.txt"), header = FALSE)[,1]
  filenames <- filenames[!(filenames %in% done_files)]
}

# tile_size
tile_size = 512

# maximum proportion of na values in a file to still be written
na_threshold = .75

# function to generate a tile from x and y starts --------------------------------------------------------------
make_tiles <- function(x_start, y_start, file){

  # throw an error if RGB is the first variable because this one won't let you know which images to remove 
  # based on not having enough lst/emmissivity/surface radiance data
  if (variables[1] == "RGB"){
    stop("RGB cannot be the first variable in the variables list")
  }

  # name of tile
  run <- substr(file, 1, nchar(file)-4) # filename minus extension
  tilename <- paste0(run, "_tile", floor(x_start/tile_size), ",", floor(y_start/tile_size), ".tif")

  variable <- variables[1]

  image <- raster::brick(here(input_dir, variable, file))
  tile <- raster::crop(image, extent(image, x_start, x_start + tile_size, y_start, y_start + tile_size))
  # check the proportion of NA values in the LST file
  prop_na <- length(values(tile)[is.na(values(tile))])/length(values(tile))

  # write the raster only if there is less than na values than the na_threshold
  if (prop_na < na_threshold){

    # write the raster for the first variable
    raster::writeRaster(tile, here(output_dir, variable, tilename), overwrite = TRUE)

    # create and write rasters for all other variables
    for (variable in variables[2:length(variables)]){
        image <- raster::brick(here(input_dir, variable, file))
        tile <- raster::crop(image, extent(image, x_start, x_start + tile_size, y_start, y_start + tile_size))
        if (variable == 'RGB'){
          raster::writeRaster(tile, here(output_dir, variable, tilename), overwrite = TRUE, datatype='INT1U') # save as 8-bit
        } else {
          raster::writeRaster(tile, here(output_dir, variable, tilename), overwrite = TRUE)
        }
    }
    print(paste("Done with tile", tilename))
  } else {
    print(paste("Insufficient cover; passing tile", tilename))
  }
}

# function to tile one image --------------------------------------------------------------

tile_image <- function(file){

  tic()

  image <- raster::brick(here(input_dir, variables[1], file))

  #-----------------------------#
  # get the extents of the tiles
    xmin <- extent(image)[1]
    ymin <- extent(image)[3]

    nrow <- dim(image)[1]
    ncol <- dim(image)[2]

    # get the number of tiles in the x direction that fits
    x_tiles <- floor(nrow/tile_size)

    # get the number of tiles in the y direction
    y_tiles <- floor(ncol/tile_size)

    # skip this run if it isn't big enough
    if (x_tiles == 0 | y_tiles == 0){
      print(paste("skipping run because it is too small for a single tile:", file))
    } else {
      x_margin <- floor((nrow - (tile_size*x_tiles))/2)
      y_margin <- floor((nrow - (tile_size*x_tiles))/2)

      x_starts <- x_margin + seq(0, x_tiles*tile_size-1, tile_size)
      y_starts <- y_margin + seq(0, y_tiles*tile_size-1, tile_size)

      # generate output folder if it doesn't exist yet
      if (!dir.exists(output_dir)){
        for (variable in variables){
          dir.create(here(output_dir, variable), recursive = TRUE)
        }
      }

      for (y_start in y_starts){
        lapply(x_starts, make_tiles, y_start, file)
      }
    }

 write(file, here(output_dir, "done_files.txt"), append = TRUE)
 print(paste0("done with run", file))
 toc()
}

# call function to tile one image on all runs --------------------------------------------------------------
# lapply(filenames, tile_image)

# in parallel --------------------------------------------------------------
no_cores <- parallel::detectCores()
cl <- makeCluster(no_cores, type="FORK")
parLapply(cl, filenames, tile_image)
stopCluster(cl)