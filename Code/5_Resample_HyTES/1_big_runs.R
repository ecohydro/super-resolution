### For runs with excessive amounts of associated PLANET files, make an individual raster for each PLANET file

library(here)
library(sf)
library(raster)
library(tictoc)
library(parallel)
library(stringr)

# Read in all relevant data
run_keys <- read.csv2(here("Data","Intermediate","run_keys.csv"))[2:4]
completed_runs <- read.table(here("Data","Final", "complete_run.txt"))[,1]
fat_runs <- read.table(here("Data","Final","fat_runs.txt"))[,1] #"20210823KirunaSE3" 20170608GarfieldCO2_822-2538
fat_runs <- fat_runs[!(fat_runs %in% completed_runs)]

# get the planet crs for a random planet raster
planet_crs <- crs(raster(here("Data","Raw","PLANET", "2021-08", "1146-1554.tiff")))

# a function that processes the hytes data over a planet image
process_by_planet <- function(planet_filename, HyTES_p_crs){
    # read in planet
    PLANET <- brick(here("Data","Raw","PLANET", rundate, planet_filename))

    # Unique planet ID
    PLANET_ID <- str_split(planet_filename, '.t')[[1]][1]

    # crop both planet and hytes to the overlapping extent
    extent_intersection <- raster::intersect(extent(PLANET), extent(HyTES_p_crs))
    if (is.null(extent_intersection)){
        print("No overlap between Planet and HyTES; skipping")
        return()
    }

    HyTES_crop <- st_crop(HyTES_p_crs, extent_intersection)
    PLANET_crop <- crop(PLANET, extent_intersection)

    if (nrow(HyTES_crop) == 0){
        print("No HyTES observations; skipping")
        return()
    } 

    #Resample Hytes to the planet raster
    HyTES_raster_ver <- rasterize(HyTES_crop, PLANET, fun = mean)

    # For some reason sometimes the HyTES isn't perfectly cropped and adds more past where there is RGB info
    HyTES_raster_ver <- crop(HyTES_raster_ver, extent(PLANET_crop))

    # generate and save layers
    LST_layer <- HyTES_raster_ver[['LST']]
    Emis_layer <- HyTES_raster_ver[['Emssvty']]
    Surf_Rad_layer <- HyTES_raster_ver[['Srfc_Rd']]

    writeRaster(LST_layer, here("Data", "Final", "LST", paste0(run_id,"_",PLANET_ID,".tif")), overwrite=TRUE)
    writeRaster(PLANET_crop, here("Data", "Final", "RGB", paste0(run_id,"_",PLANET_ID,".tif")), overwrite=TRUE)
    writeRaster(Emis_layer, here("Data", "Final", "Emis", paste0(run_id,"_",PLANET_ID,".tif")), overwrite=TRUE)
    writeRaster(Surf_Rad_layer, here("Data", "Final", "Srfc_Rad", paste0(run_id,"_",PLANET_ID,".tif")), overwrite=TRUE)

    print(paste0(run_id, ": Finished ", PLANET_ID, " out of ", length(file_names)))

    return()
    
}

for(i in 1:length(fat_runs)){
    run_id <- fat_runs[i]
    run_info <- run_keys[run_keys$run_id == run_id, ]
    if (run_id %in% completed_runs){
        print(paste(run_id, "is already completed, next!"))
        next
    }
    print(paste(run_id,"new run!"))
    shapefile_dir <- here("Data","Intermediate","HyTES_sf", run_id, paste0(run_id,".shp"))
    # If the shapefile doesn't exist, move on to next observation
    if (file.exists(shapefile_dir) == FALSE){
        print(paste(run_id, "skipped!"))
        next
    }
    rundate <- run_info[2]
    file_names <- run_info[3]
    if (substr(file_names, 1, 1) == 'c'){
        file_names <- unlist(strsplit(substr(run_info[3], 3, nchar(run_info[3])-1), ", "))
    }

    HyTES_sf <- st_read(shapefile_dir)
    HyTES_p_crs <- st_transform(HyTES_sf, planet_crs)

    # process everything in parallel. May need to reduce the number of cores to prevent the system from running out of memory
    no_cores <- 8 #parallel::detectCores()
    cl <- makeCluster(no_cores, type="FORK")
    parLapply(cl, file_names, process_by_planet, HyTES_p_crs)
    stopCluster(cl)

    # if you don't want to process in parallel, go for sequential:
    # for (file in file_names){
    #     process_by_planet(file, HyTES_p_crs)
    # }

    write(run_id, file = here("Data","Final", "complete_run.txt"), append = TRUE)
}