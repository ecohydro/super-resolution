#####################
# Take the target and basemap files and aggregate and resample them to the input and output size. 
# Files saved under Tiles/Input and Tiles/Output
# Anna Boser Aug 7 2022
#####################

library(here)
library(raster)
library(parallel)

# parameters and dirs -------------------------------------------------------------- 

#input directories
input_dir <- here("Data", "Tiles")
target = "LST"
basemap = "RGB"

# tile_size of the inputs
tile_size = 672 # divisible by 14 (ECOSTRESS to 5m) and 24 (Landsat)

# output dir
output_dir = here("Data", "For_CNN")

# list all the files
filenames <- c()
for (variable in c(target, basemap)){
    filenames <- c(filenames, unlist(list.files(here(input_dir, paste0(variable, tile_size)))))
}
print(length(filenames)/2)
filenames = filenames[duplicated(filenames)] # ensure they are in all variable folders -- for example this might not be the case if we dropped some RGB images because of missing data (see 0_tile.R)
print(length(filenames))

# input and output size
input_size = 70
output_size = 10
pixel_size = 5

# if this script has already been run and some runs are already processed: 
if (file.exists(here(output_dir, paste0("done_files", tile_size, "_", input_size, "_", output_size, ".txt")))){
  done_files <- read.csv(paste0("done_files", tile_size, "_", input_size, "_", output_size, ".txt"), header = FALSE)[,1]
  filenames <- filenames[!(filenames %in% done_files)]
}

input_scale = as.integer(input_size/pixel_size)
output_scale = as.integer(output_size/pixel_size)

# function to create the input and output from a file --------------------------------------------------------------

in_out_put <- function(file){

    # read in the target file
    target_im <- raster::raster(here(input_dir, paste0(target, tile_size), file))

    # read in the basemap
    basemap_im <- raster::brick(here(input_dir, paste0(basemap, tile_size), file))

    # generate the model output -------------------
    target_outsize = raster::aggregate(target_im, output_scale)

    # save the model output
    output_dir_outputs = here(output_dir, "outputs", paste0(target, "_", tile_size, "_", output_size))
    if (!dir.exists(output_dir_outputs)){
        dir.create(output_dir_outputs, recursive = TRUE)
    }
    raster::writeRaster(target_outsize, here(output_dir_outputs, file), overwrite = TRUE)

    # generate the model inputs -------------------

    # aggregated and resized target
    target_insize = raster::aggregate(target_im, input_scale)
    target_insize = raster::resample(target_insize, target_outsize, method = "ngb")

    # save the model target input
    output_dir_target_inputs = here(output_dir, "inputs", paste0(target, "_", tile_size, "_", input_size, "_", output_size))
    if (!dir.exists(output_dir_target_inputs)){
        dir.create(output_dir_target_inputs, recursive = TRUE)
    }
    raster::writeRaster(target_insize, here(output_dir_target_inputs, file), overwrite = TRUE)

    # aggregate basemap
    basemap_outsize = raster::aggregate(basemap_im, output_scale)

    # save the model basemap input
    output_dir_basemap_inputs = here(output_dir, "inputs", paste0(basemap, "_", tile_size, "_", output_size))
    if (!dir.exists(output_dir_basemap_inputs)){
        dir.create(output_dir_basemap_inputs, recursive = TRUE)
    }
    raster::writeRaster(basemap_outsize, here(output_dir_basemap_inputs, file), overwrite = TRUE, datatype='INT1U') # 8-bit

    # add the mean and standard deviation of each band to a running file -------------------

    # target info
    target_mean_sd = here(output_dir, paste0("target_mean_sd_files", tile_size, "_", input_size, "_", output_size, ".txt"))
    if (!file.exists(target_mean_sd)){
        write(paste("file", "mean", "sd"), target_mean_sd, append = TRUE)
    }
    write(paste(file, mean(target_im[], na.rm = TRUE), sd(target_im[], na.rm = TRUE)), target_mean_sd, append = TRUE)

    # basemap info
    basemap_mean_sd = here(output_dir, paste0("basemap_mean_sd_files", tile_size, "_", input_size, "_", output_size, ".txt"))
    if (!file.exists(basemap_mean_sd)){
        means = paste(paste0("mean", 1:dim(basemap_im)[3]), collapse = " ")
        sds = paste(paste0("sd", 1:dim(basemap_im)[3]), collapse = " ")
        write(paste("file", means, sds), basemap_mean_sd, append = TRUE)
    }
    means = paste(lapply(1:dim(basemap_im)[3], function(n){mean(basemap_im[[n]][], na.rm = TRUE)}), collapse = " ")
    sds = paste(lapply(1:dim(basemap_im)[3], function(n){sd(basemap_im[[n]][], na.rm = TRUE)}), collapse = " ")
    write(paste(file, means, sds), basemap_mean_sd, append = TRUE)
    

    # add this file to the list of processed files -------------------
    write(file, here(output_dir, paste0("done_files", tile_size, "_", input_size, "_", output_size, ".txt")), append = TRUE)
    print(paste0("done with run", file))
}

# call function to tile one image on all runs --------------------------------------------------------------
for (file in filenames){
    in_out_put(file)
}
# lapply(filenames, in_out_put)

# in parallel -- this is super fast but also will write so fast that appending to files gets messed up--------------------------------------------------------------
# no_cores <- parallel::detectCores()
# cl <- makeCluster(no_cores, type="FORK")
# parLapply(cl, filenames, in_out_put)
# stopCluster(cl)