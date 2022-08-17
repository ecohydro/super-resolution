### For each run, stitch relevant PLANET files and convert to rasters for model

library(here)
library(sf)
library(raster)
library(tictoc)

# Crop the PLANET data to the extent of the run
run_keys <- read.csv2(here("Data","Intermediate","run_keys_3.csv"))[2:4]
completed_runs <- read.table(here("Data","Final", "complete_run.txt"))[,1]
fat_runs <- read.table(here("Data","Final","fat_runs.txt"))[,1]

counter <- 0 

for (j in 1:nrow(run_keys)){
    run_id <- run_keys[j,1]
    # If the run is already completed, move on to the next run
    if (run_id %in% completed_runs){
        print(paste(run_id, "is already completed, next!"))
        next
    }
    # If the run is too big, move on to the next run
    if (run_id %in% fat_runs){
        print(paste(run_id, "is already too big, next!"))
        next
    }
    print(paste(run_id,"new run!"))
    shapefile_dir <- here("Data","Intermediate","HyTES_sf", run_id, paste0(run_id,".shp"))
    # If the shapefile doesn't exist, move on to next observation
    if (file.exists(shapefile_dir) == FALSE){
        print(paste(run_id, "has no shapefile, skipped!"))
        next
    }
    rundate <- run_keys[j,2]
    file_names <- run_keys[j,3]
    # If there is more than one associated PLANET file, make a proper vector of the multiple files
    if (substr(file_names, 1, 1) == 'c'){
        file_names <- unlist(strsplit(substr(run_keys[j,3], 3, nchar(run_keys[j,3])-1), ", "))
    }
    first = TRUE
    # Creates PLANET data set by initializing brick with the first observation. Then, merges more PLANET data to original brick.
    tic()
    if (length(file_names) > 4){
        write(run_id, file = here("Data","Final", "fat_runs.txt"), append = TRUE)
        print(paste(run_id, "is too big, it's been recorded. next!"))
        next
    }
    for (i in 1:length(file_names)){
        if (first == FALSE){
            PLANET_old <- brick(here("Data","Raw","PLANET", rundate, file_names[i]))
            PLANET <- merge(PLANET, PLANET_old)
        }
        if (first == TRUE){
            PLANET <- brick(here("Data","Raw","PLANET", rundate, file_names[i]))
            first = FALSE
        }
    }
    toc()
    # Reads in shapefile, transforms HyTES data to same crs as PLANET and crops PLANET to the extent of HyTES data
    tic()
    HyTES_sf <- st_read(shapefile_dir)
    toc()
    tic()
    HyTES_new <- st_transform(HyTES_sf, crs(PLANET))
    toc()
    tic()
    extent_intersection <- raster::intersect(extent(PLANET), extent(HyTES_new))
    if (is.null(extent_intersection)){
        print("No overlap between Planet and HyTES; skipping")
        next
    }
    PLANET_crop <- crop(PLANET, extent_intersection)
    toc()
    # Resample the Emissivity LST and surface radiance to the ~5m grid
    tic()
    HyTES_raster_ver <- rasterize(HyTES_new, PLANET_crop, fun = mean)
    toc()
    # Save a raster brick of 3 layers for variables defined above
    LST_layer <- HyTES_raster_ver[['LST']]
    Emis_layer <- HyTES_raster_ver[['Emssvty']]
    Surf_Rad_layer <- HyTES_raster_ver[['Srfc_Rd']]
    # Take save a brick of it with RGB from PLANET (3 layers)
    PLANET_raster <- PLANET_crop
    writeRaster(LST_layer, here("Data", "Final", "LST", paste0(run_id,".tif")))
    writeRaster(PLANET_raster, here("Data", "Final", "RGB", paste0(run_id,".tif")))
    writeRaster(Emis_layer, here("Data", "Final", "Emis", paste0(run_id, ".tif")))
    writeRaster(Surf_Rad_layer, here("Data", "Final", "Srfc_Rad", paste0(run_id, ".tif")))
    write(run_id, file = here("Data","Final", "complete_run.txt"), append = TRUE)
    print(paste(run_id, "done!"))
    counter <- counter + 1
    print(paste("We made", counter, "raster files so far!"))
}
