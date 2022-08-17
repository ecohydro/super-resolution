### This file converts .hdf5 data frames into shape files and saves them in the same folder.

library(here)
library(dplyr)
library(hdf5r)
library(raster)
library(reshape)
library(sf)
library(tictoc)

# Reads in part of segment_metadata.csv file to get unqiue run_IDs
# Metadatas completed : _1, _2, 
run_name <- unique(read.csv(here("Data","Segment_Metadata","segment_metadata_3.csv"))[3])

# Planck's function
Planck_fun <- function(emis, temp, wave){
    surface_radiance <- emis * ((3.74 * 10^(-16)) / ((wave*10^(-6))^5)) * (1/(exp(.0143879/((wave*10^(-6)) * temp))-1)) * 10^(-6)
    return(surface_radiance)
}

for (i in 1:nrow(run_name)){
    tic()
    print(paste("New Run:", run_name[i,]))
    # If the file is too old, move on to the next file
    if (substr(run_name[i,], 4, 4) < '7' & substr(run_name[i,], 4, 4) > '3'){
        print(paste(run_name[i,], "is too old, next!"))
    }
    sf_name <- paste0(run_name[i,],".shp")
    # Checks to see if shapefile has already been created
    if (file.exists(here("Data", "Intermediate", "HyTES_sf", run_name[i,], sf_name))){
        print(paste(run_name[i,],"already done, next!"))
        next
    }
    # Initialize empty data frame for each run
    run_df <- data.frame(matrix(ncol = 7, nrow = 0))
    # Extract all file names of each run folder
    file_names <- list.files(here("Data","Intermediate","HyTES", run_name[i,]))
    for (j in 1:length(file_names)){
        # Extracts directory name for each L1/L2 file, or .csv if file is from 2016
        file_dir <- here("Data", "Intermediate", "HyTES", run_name[i,], file_names[j])
        # Special case for runs from 2016 or older
        #if (substr(run_name[i,], 4, 4) == '6'){
        #     is_L2 <- grepl('.L2', file_names[j])
        #     is_alt <- grepl('geoAlt', file_names[j])
        #     is_lat <- grepl('geoLat', file_names[j])
        #     is_lon <- grepl('geoLong', file_names[j])
        #     if (is_alt){
        #         alt <- read.csv(file_dir)
        #         alt_vec <- alt %>% melt() %>% dplyr::rename(Height = value) %>% subset(select= 2)
        #     }
        #     if (is_lat){
        #         lat <- read.csv(file_dir)
        #         lat_vec <- lat %>% melt() %>% dplyr::rename(Latitude = value) %>% subset(select= 2)
        #     }
        #     if (is_lon){
        #         lon <- read.csv(file_dir)
        #         lon_vec <- lon %>% melt() %>% dplyr::rename(Longitude = value) %>% subset(select= 2)
        #         geo_loc_df <- data.frame(lat_vec, lon_vec, alt_vec)
        #     }
        #     if (is_L2){
        #         # Opens L2 HDF5 file
        #         L2_data <- H5File$new(file_dir, mode = 'r')
        #         if ('l2_land_surface_temperature' %in% names(L2_data)){
        #             # Saves LST data as a dataframe
        #             LST_df <- data.frame(L2_data[['l2_land_surface_temperature']][,])
        #             # Saves Emissivity data as a dataframe (around 11.5 micrometers)
        #             emissivity_df <- data.frame(L2_data[['l2_emissivity']][202,,])
        #             # Saves exact wavelength value
        #             lambda <- L2_data[['l2_emissivity_wavelengths']][202]
        #             # Checks if L1 and L2 data frames have same dimensions
        #             if (ncol(LST_df) != L1_data[['pixel_geolocation']]$dims[3]){
        #                 # Removes last n columns of L2 data frame to match dim of L1 data frame
        #                 for (i in 1:(ncol(LST_df)-L1_data[['pixel_geolocation']]$dims[3])){
        #                     LST_df_new <- LST_df_new %>% dplyr::select(-last_col())
        #                     emissivity_df_new <- emissivity_df_new %>% dplyr::select(-last_col())
        #                 }
        #             }
        #         }
        #         else{
        #             # Saves LST data as a dataframe
        #             LST_df <- data.frame(L2_data[['L2_LST']][,])
        #             # Saves Emissivity data as a dataframe (around 11.5 micrometers)
        #             emissivity_df <- data.frame(L2_data[['L2_Emissivity']][202,,])
        #             # Saves exact wavelength value
        #             lambda <- L2_data[['L2_Emissivity_Wavelengths']][202]
        #         }
        #         # Setting first and last 5 rows (bands) to NA for future trimming
        #         LST_df[1:5,] <- NA
        #         LST_df[(nrow(LST_df)-4):nrow(LST_df),] <- NA
        #         # Melts both LST and Emissivity data to one column vectors
        #         LST_vec <- LST_df %>% melt() %>% dplyr::rename(LST = value) %>% subset(select= 2)
        #         emissivity_vec <- emissivity_df %>% melt() %>% dplyr::rename(Emissivity = value) %>% subset(select= 2)
        #         # Apply Planck function to get surface radiance
        #         surf_rad_vec <- mapply(Planck_fun, emissivity_vec, LST_vec, lambda)
        #         colnames(surf_rad_vec) <- c('Surface_Radiance')
        #         # Trim edges and combines L1 and L2 data to a data frame if they are of the same dimensions
        #         if(nrow(geo_loc_df) == nrow(LST_vec)){
        #             loc_LST_emis <- geo_loc_df %>% cbind(LST_vec, emissivity_vec, surf_rad_vec) %>% na.omit()
        #             # Adds segment data to run data frame
        #             run_df <- run_df %>% rbind(loc_LST_emis)
        #         }
        #         # Else it removes the observation from the run data frame, but stores file name to a text file
        #         else{
        #             write(run_name[i,], file = here("Data","Intermediate", "wrong_dim.txt"), append = TRUE) 
        #         }
        #     }
        #}
        #else {
            # Boolean indicator for L1 and L2
            is_L1 <- grepl('_L1', file_names[j])
            is_L2 <- grepl('_L2', file_names[j])
            if (is_L1) {
                # Opens L1 HDF5 file
                L1_data <- H5File$new(file_dir, mode = 'r')
                # Saves pixel geolocation as a data frame
                geo_loc_df <- data.frame(Latitude = as.vector(L1_data[['pixel_geolocation']][1,,]),
                                         Longitude = as.vector(L1_data[['pixel_geolocation']][2,,]),
                                         Height = as.vector(L1_data[['pixel_geolocation']][3,,]),
                                         Steps = as.vector(L1_data[['pixel_geolocation']][4,,]))                       
            }
            else if (is_L2) {
                # Opens L2 HDF5 file
                L2_data <- H5File$new(file_dir, mode = 'r')
                if ('l2_land_surface_temperature' %in% names(L2_data)){
                    # Saves LST data as a dataframe
                    LST_df <- data.frame(L2_data[['l2_land_surface_temperature']][,])
                    # Saves Emissivity data as a dataframe (around 11.5 micrometers)
                    emissivity_df <- data.frame(L2_data[['l2_emissivity']][202,,])
                    # Saves exact wavelength value
                    lambda <- L2_data[['l2_emissivity_wavelengths']][202]
                    # Checks if L1 and L2 data frames have same dimensions
                    if (ncol(LST_df) != L1_data[['pixel_geolocation']]$dims[3]){
                        # Removes last n columns of L2 data frame to match dim of L1 data frame
                        for (i in 1:(ncol(LST_df)-L1_data[['pixel_geolocation']]$dims[3])){
                            LST_df_new <- LST_df_new %>% dplyr::select(-last_col())
                            emissivity_df_new <- emissivity_df_new %>% dplyr::select(-last_col())
                        }
                    }
                }
                else{
                    # Saves LST data as a dataframe
                    LST_df <- data.frame(L2_data[['L2_LST']][,])
                    # Saves Emissivity data as a dataframe (around 11.5 micrometers)
                    emissivity_df <- data.frame(L2_data[['L2_Emissivity']][202,,])
                    # Saves exact wavelength value
                    lambda <- L2_data[['L2_Emissivity_Wavelengths']][202]
                }
                # Setting first and last 5 rows (bands) to NA for future trimming
                LST_df[1:5,] <- NA
                LST_df[(nrow(LST_df)-4):nrow(LST_df),] <- NA
                # Melts both LST and Emissivity data to one column vectors
                LST_vec <- LST_df %>% melt() %>% dplyr::rename(LST = value) %>% subset(select= 2)
                emissivity_vec <- emissivity_df %>% melt() %>% dplyr::rename(Emissivity = value) %>% subset(select= 2)
                # Apply Planck function to get surface radiance
                surf_rad_vec <- mapply(Planck_fun, emissivity_vec, LST_vec, lambda)
                colnames(surf_rad_vec) <- c('Surface_Radiance')
                # Trim edges and combines L1 and L2 data to a data frame if they are of the same dimensions
                if(nrow(geo_loc_df) == nrow(LST_vec)){
                    loc_LST_emis <- geo_loc_df %>% cbind(LST_vec, emissivity_vec, surf_rad_vec) %>% na.omit()
                    # Adds segment data to run data frame
                    run_df <- run_df %>% rbind(loc_LST_emis)
                }
                # Else it removes the observation from the run data frame, but stores file name to a text file
                else{
                    write(run_name[i,], file = here("Data","Intermediate", "wrong_dim.txt"), append = TRUE) 
                }
            }
        #}
    }
    print(paste("L1 and L2 data done:", run_name[i,]))
    # Convert data frame to shapefile and filters out extreme LST values
    run_sf <- st_as_sf(x = run_df, coords = c("Longitude", "Latitude"), crs = "WGS84")
    filtered_sf <- run_sf %>% filter(LST > 250 & LST < 375)
    # Write shapefile to a new Intermediate folder called HyTES_sf 
    dir.create(here("Data", "Intermediate", "HyTES_sf", run_name[i,]))
    st_write(filtered_sf, here("Data", "Intermediate", "HyTES_sf", run_name[i,], sf_name))
    print(paste("Run Done:", run_name[i,]))
    toc()
}

