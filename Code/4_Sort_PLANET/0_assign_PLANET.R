### This file finds which PLANET files to use for each run

library(here)
library(raster)
library(sf)
library(lubridate)

#### NOTE: IF YOU REMOVE COLUMNS, BE SURE TO CORRECT THE CODE SO IT REFLECTS THOSE CHANGES (aka whenever the metadata file is read/used)
metadata <- read.csv(here("Data","Run_Metadata","run_metadata.csv"))
x_min <- read.csv(here("Data","Run_Metadata","run_metadata.csv"))[7]
x_max <- read.csv(here("Data","Run_Metadata","run_metadata.csv"))[8]
y_min <- read.csv(here("Data","Run_Metadata","run_metadata.csv"))[5]
y_max <- read.csv(here("Data","Run_Metadata","run_metadata.csv"))[6]

### Write .csv file containing PLANET name, date, and extent coords
PLANET <- data.frame(matrix(ncol = 6, nrow = 0))

# Converts PLANET .tiff files to 4 point rasters to get extent coords in proper form
for (k in 1:length(list.files(path = here("Data","Raw","PLANET")))){
    date <- list.files(path = here("Data","Raw","PLANET"))[k]
    for (i in 1:length(list.files(path = here("Data","Raw","PLANET", date)))){
        file_name <- list.files(path = here("Data","Raw","PLANET", date))[i]
        file_dir <- here("Data","Raw","PLANET", date, file_name)
        example <- brick(file_dir)
        four_point <- raster(extent(example))
        crs(four_point) <- "+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null +wktext +no_defs"
        projection <- projectRaster(four_point, crs = '+proj=longlat +datum=WGS84')
        x_min_pt <- as.numeric(extent(projection)[1])
        x_max_pt <- as.numeric(extent(projection)[2])
        y_min_pt <- as.numeric(extent(projection)[3])
        y_max_pt <- as.numeric(extent(projection)[4])
        observation <- c(file_name, date, x_min_pt, x_max_pt, y_min_pt, y_max_pt)
        PLANET <- PLANET %>% rbind(observation)
    }
}

# Finalize PLANET data frame
colnames(PLANET) <- c('file_name','date','x_min','x_max','y_min','y_max')
PLANET <- PLANET %>% transform(x_min = as.numeric(x_min),
                               x_max = as.numeric(x_max),
                               y_min = as.numeric(y_min),
                               y_max = as.numeric(y_max))

# Create a polygon of the extents for each HyTES observation
lst <- lapply(1:nrow(metadata), function(x){
  ## create a matrix of coordinates that also 'close' the polygon
  res <- matrix(c(metadata[x, 'Min_Longitude'], metadata[x, 'Min_Latitude'],
           metadata[x, 'Min_Longitude'], metadata[x, 'Max_Latitude'],
           metadata[x, 'Max_Longitude'], metadata[x, 'Max_Latitude'],
           metadata[x, 'Max_Longitude'], metadata[x, 'Min_Latitude'],
           metadata[x, 'Min_Longitude'], metadata[x, 'Min_Latitude'])  ## need to close the polygon
         , ncol = 2, byrow = T)
  ## create polygon objects
  st_polygon(list(res))
})

# Create a polygon of the extents for each PLANET observation
planet_lst <- lapply(1:nrow(PLANET), function(x){
  ## create a matrix of coordinates that also 'close' the polygon
  res <- matrix(c(PLANET[x, 'x_min'], PLANET[x, 'y_min'],
           PLANET[x, 'x_min'], PLANET[x, 'y_max'],
           PLANET[x, 'x_max'], PLANET[x, 'y_max'],
           PLANET[x, 'x_max'], PLANET[x, 'y_min'],
           PLANET[x, 'x_min'], PLANET[x, 'y_min'])  ## need to close the polygon
         , ncol = 2, byrow = T)
  ## create polygon objects
  st_polygon(list(res))
})

# Retrieve the date of HyTES run in the form of year-month
get_date <- sapply(1:nrow(metadata), function(x){
    yr_mo <- substr(metadata[x,4], 1, 7)
})

# Generate shapefiles for HyTES and PLANET with ID, date, and polygon of extent
HyTES_sf <- st_sf(Run_ID = metadata[,'Run_ID'], date = get_date, st_sfc(lst))
PLANET_sf <- st_sf(file_name = PLANET[,'file_name'], date = PLANET[,'date'], st_sfc(planet_lst))

# Initializes list
file_names <- list()

# Assigns correct PLANET files for each HyTES run
for (i in 1:nrow(HyTES_sf)){
    PLANET_files <- c()
    # Goes through all PLANET files for each run
    for (j in 1:nrow(PLANET_sf)){
        if(HyTES_sf$date[i] == PLANET_sf$date[j]){
            # If the HyTES run is contained in a single shape file, record name of shape file and break
            if (st_contains(st_geometry(PLANET_sf[j,]), st_geometry(HyTES_sf[i,]), sparse = FALSE) == TRUE){
                PLANET_files <- c(PLANET_sf$file_name[j])
                break
            }            
            # If there is any overlap between PLANET file and HyTES file, record name of PLANET file and append to vector
            if (st_overlaps(st_geometry(HyTES_sf[i,]), st_geometry(PLANET_sf[j,]), sparse = FALSE) == TRUE){
                PLANET_files <- c(PLANET_files, PLANET_sf$file_name[j])
            }
        }
    }
    print(PLANET_files)
    # If there are no PLANET files associated to run, label as NA
    if (length(PLANET_files) == 0){
        file_names[[i]] <- NA
    }
    # Set file_names to observation
    else{
        file_names[[i]] <- PLANET_files
    }
}

# Make a data frame containing run id, date in the form of year-month, and associated PLANET file names
new_run_keys <- data.frame('run_id' = metadata[,'Run_ID'], 'date' = get_date, 'file_names' = I(file_names))

# Write to .csv file
write.csv2(new_run_keys, here("Data","Intermediate","run_keys.csv"))