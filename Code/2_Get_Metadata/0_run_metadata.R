### This file extracts all relevant metadata information for each run writes information to Run_Metadata folder.

library(here)
library(dplyr)
library(stringr)
library(tidyr)
library(lubridate)
library(hdf5r)

# Creata master segment metadata dataframe
for (sheet in list.files(here("Data","Segment_Metadata"))){
    # Read in .csv file
    seg_df <- read.csv(here("Data","Segment_Metadata", sheet))
    # Initialize master segment data set
    if (sheet == 'segment_metadata_1.csv'){
        final_seg_df <- seg_df
        seg_df <- seg_df[-(1:nrow(seg_df)),]
    }
    # Append previous .csv file to master data frame
    final_seg_df <- final_seg_df %>% rbind(seg_df)
}

# Reads in all unique run IDs and stores them
run_IDs <- list.files(here("Data","Intermediate","HyTES"))

# Initiatlizes final run metadata dataframe
final_metadata <- data.frame(matrix(ncol = 11, nrow = 0))

# Creates final run metadata dataframe
for (run in run_IDs){
    print(paste("Start:", run))
    # Creates individual run metadata dataframe
    grouped_seg_df <- data.frame(matrix(ncol = 12, nrow = 0))
    colnames(grouped_seg_df) <- colnames(final_seg_df)
    # Appends all segments together in one dataframe
    grouped_seg_df <- grouped_seg_df %>% rbind(final_seg_df[final_seg_df$Run == run, ])
    # Adds Time column
    grouped_seg_df <- grouped_seg_df %>%
                    mutate(Time = make_datetime(Year, Month, Day, Hour, Minute, Second)) %>%
                    subset(select = c(10:13))
    Min_Time <- as.numeric(min(grouped_seg_df$Time))
    Max_Time <- as.numeric(max(grouped_seg_df$Time))
    County <- grouped_seg_df$County[1]
    State <- grouped_seg_df$State[1]
    Aircraft <- grouped_seg_df$Aircraft[1]
    Year <- substr(run, 1,4) %>% as.numeric()
    segment_min_lat <- c()
    segment_max_lat <- c()
    segment_min_long <- c()
    segment_max_long <- c()
    file_names <- list.files(here("Data", "Intermediate", "HyTES", run))
    if (substr(run, 4, 4) > '3' & substr(run, 4, 4) < '7'){
        lat_files <- dir(here("Data", "Intermediate", "HyTES", run), "geoLat.csv$", ignore.case = TRUE, all.files = TRUE)
        long_files <- dir(here("Data", "Intermediate", "HyTES", run), "geoLong.csv$", ignore.case = TRUE, all.files = TRUE)
        for (lat_file in lat_files){
            lat_table <- read.csv(here("Data", "Intermediate", "HyTES", run, lat_file)) 
            segment_min_lat <- segment_min_lat %>% append(min(lat_table))
            segment_max_lat <- segment_max_lat %>% append(max(lat_table))
        }
        for (long_file in long_files){
            long_table <- read.csv(here("Data", "Intermediate", "HyTES", run, long_file)) 
            segment_min_long <- segment_min_long %>% append(min(long_table))
            segment_max_long <- segment_max_long %>% append(max(long_table))
        }
    }
    # Read the L1 data for Lat/Lon information
    else {
        for (j in 1:length(file_names)){
            is_L1 <- grepl('_L1', file_names[j])
            if (is_L1) {
                L1_dir <- here("Data", "Intermediate", "HyTES", run, file_names[j])
                L1_data <- H5File$new(L1_dir, mode = "r")
                location_df <- data.frame(Latitude = as.vector(L1_data[['pixel_geolocation']][1,,]),
                                          Longitude = as.vector(L1_data[['pixel_geolocation']][2,,]),
                                          Height = as.vector(L1_data[['pixel_geolocation']][3,,]),
                                          Steps = as.vector(L1_data[['pixel_geolocation']][4,,]))
                segment_min_lat <- segment_min_lat %>% append(min(location_df$Latitude))
                segment_max_lat <- segment_max_lat %>% append(max(location_df$Latitude))
                segment_min_long <- segment_min_long %>% append(min(location_df$Longitude))
                segment_max_long <- segment_max_long %>% append(max(location_df$Longitude))
            }
        }
    } 
    Min_Lat <- min(segment_min_lat)
    Max_Lat <- max(segment_max_lat)
    Min_Lon <- min(segment_min_long)
    Max_Lon <- max(segment_max_long)
    run_metadata <- c(run, Min_Time, Max_Time, Min_Lat, Max_Lat, Min_Lon, Max_Lon, County, State, Aircraft, Year)
    final_metadata <- final_metadata %>% rbind(run_metadata)
}

# Assigns column names
colnames(final_metadata) <- c("Run_ID", "Earliest_Time", "Latest_Time", "Min_Latitude", "Max_Latitude", 
                              "Min_Longitude", "Max_Longitude", "County", "State", "Aircraft", "Year")

# Changes the display of time variables into proper date time format with UTC time zone
final_metadata$Earliest_Time <- as.numeric(final_metadata$Earliest_Time)
final_metadata$Latest_Time <- as.numeric(final_metadata$Latest_Time)
final_metadata$Earliest_Time <- as.POSIXct(final_metadata$Earliest_Time, origin = lubridate::origin, tz = "UTC")
final_metadata$Latest_Time <- as.POSIXct(final_metadata$Latest_Time, origin = lubridate::origin, tz = "UTC")

# Adds Polygon object to each run observation in order to create the GeoJSON file
#first = TRUE
#for (i in 1:nrow(final_metadata)){
#    lon <- c(final_metadata$Min_Longitude[i], final_metadata$Max_Longitude[i])
#    lat <- c(final_metadata$Min_Latitude[i], final_metadata$Max_Latitude[i])
#    Poly_Coord_df <- data.frame(lon, lat)
#    Polygon <- Poly_Coord_df %>% 
#        st_as_sf(coords = c("lon", "lat"), 
#           crs = "WGS84") %>% 
#        st_bbox() %>% 
#        st_as_sfc()
#    if (first){
#        run_metadata_sf <- st_sf(final_metadata[i,], Polygon)
#        first = FALSE
#    }
#    else{
#        new_metadata_sf <- st_sf(final_metadata[i,], Polygon)
#        run_metadata_sf <- run_metadata_sf %>% rbind(new_metadata_sf)
#    }
#}
 
# Write a shapefile of metadata information 
# st_write(run_metadata_sf, here("Data", "metadata.shp"))

# Write final run metadata to a .csv file
write.csv(final_metadata, here("Data","run_metadata.csv"))
