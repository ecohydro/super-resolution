### For runs with excessive amounts of associated PLANET files, make an individual raster for each PLANET file

library(here)
library(sf)
library(raster)
library(tictoc)

# Read in all relevant data
run_keys <- read.csv2(here("Data","Intermediate","run_keys_2.csv"))[2:4]
completed_runs <- read.table(here("Data","Final", "complete_run.txt"))[,1]
fat_runs <- "20210823KirunaSE2" #read.table(here("Data","Final","fat_runs.txt"))[,1]


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
    first <- TRUE
    for (j in 1:length(file_names)){
        if (first == FALSE){
            PLANET <- brick(here("Data","Raw","PLANET", rundate, file_names[j]))
            HyTES_crop <- st_crop(HyTES_new, extent(PLANET))
            PLANET_crop <- crop(PLANET, extent(HyTES_new))
            HyTES_raster_ver <- rasterize(HyTES_crop, PLANET, fun = mean)
            LST_layer <- HyTES_raster_ver[[4]]
            Emis_layer <- HyTES_raster_ver[[5]]
            Surf_Rad_layer <- HyTES_raster_ver[[6]]
            writeRaster(LST_layer, here("Data", "Final", "LST", paste0(run_id,"_",j,".tif")))
            writeRaster(PLANET_crop, here("Data", "Final", "RGB", paste0(run_id,"_",j,".tif")))
            writeRaster(Emis_layer, here("Data", "Final", "Emis", paste0(run_id,"_",j,".tif")))
            writeRaster(Surf_Rad_layer, here("Data", "Final", "Srfc_Rad", paste0(run_id,"_",j,".tif")))
            print(paste0(run_id, ": Finished ", j, " out of ", length(file_names)))
        }
        else{
            PLANET <- brick(here("Data","Raw","PLANET", rundate, file_names[j]))
            HyTES_sf <- st_read(shapefile_dir)
            HyTES_new <- st_transform(HyTES_sf, crs(PLANET))
            HyTES_crop <- st_crop(HyTES_new, extent(PLANET))
            PLANET_crop <- crop(PLANET, extent(HyTES_new))
            HyTES_raster_ver <- rasterize(HyTES_crop, PLANET, fun = mean)
            LST_layer <- HyTES_raster_ver[[4]]
            Emis_layer <- HyTES_raster_ver[[5]]
            Surf_Rad_layer <- HyTES_raster_ver[[6]]
            writeRaster(LST_layer, here("Data", "Final", "LST", paste0(run_id,"_",j,".tif")))
            writeRaster(PLANET_crop, here("Data", "Final", "RGB", paste0(run_id,"_",j,".tif")))
            writeRaster(Emis_layer, here("Data", "Final", "Emis", paste0(run_id,"_",j,".tif")))
            writeRaster(Surf_Rad_layer, here("Data", "Final", "Srfc_Rad", paste0(run_id,"_",j,".tif")))
            print(paste0(run_id, ": Finished ", j, " out of ", length(file_names)))
            first <- FALSE
        }
    }
    write(run_id, file = here("Data","Final", "complete_run.txt"), append = TRUE)
}
