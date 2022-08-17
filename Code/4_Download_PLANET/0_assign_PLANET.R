### This file finds which PLANET files to use for each run

library(here)
library(raster)
library(sf)
library(lubridate)

#### NOTE: IF YOU REMOVE COLUMNS, BE SURE TO CORRECT THE CODE SO IT REFLECTS THOSE CHANGES (aka whenever the metadata file is read/used)
metadata <- read.csv(here("Data","Run_Metadata","run_metadata_final.csv"))
x_min <- read.csv(here("Data","Run_Metadata","run_metadata_final.csv"))[8]
x_max <- read.csv(here("Data","Run_Metadata","run_metadata_final.csv"))[9]
y_min <- read.csv(here("Data","Run_Metadata","run_metadata_final.csv"))[6]
y_max <- read.csv(here("Data","Run_Metadata","run_metadata_final.csv"))[7]


locate <- function(point, time){
    i = 1
    for(j in 1:nrow(PLANET)){
        within_x <- (as.numeric(point[1]) > as.numeric(PLANET[j,3])) & (as.numeric(point[1]) < as.numeric(PLANET[j,4]))
        within_y <-  (as.numeric(point[2]) > as.numeric(PLANET[j,5])) & (as.numeric(point[2]) < as.numeric(PLANET[j,6]))
        right_date <- time == PLANET[j, 2]
        if (within_x & within_y & right_date){
            break
        }
        else{
            i = i + 1 
        } 
    }
    return(i)
}

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

colnames(PLANET) <- c('file_name','date','x_min','x_max','y_min','y_max')
PLANET <- PLANET %>% transform(x_min = as.numeric(x_min),
                               x_max = as.numeric(x_max),
                               y_min = as.numeric(y_min),
                               y_max = as.numeric(y_max))

run_keys <- data.frame(matrix(ncol = 2, nrow = 0))

for(i in 1:nrow(metadata)){
    point_1 <- c(x_min[i,1], y_min[i,1])
    point_2 <- c(x_min[i,1], y_max[i,1])
    point_3 <- c(x_max[i,1], y_min[i,1])
    point_4 <- c(x_max[i,1], y_max[i,1])
    yr_mo <- substr(metadata[i,4], 1, 7)
    index_1 <- locate(point_1, yr_mo)
    index_2 <- locate(point_2, yr_mo)
    index_3 <- locate(point_3, yr_mo)
    index_4 <- locate(point_4, yr_mo)
    if (index_1 == (nrow(PLANET)+1) | index_2 == (nrow(PLANET)+1) | index_3 == (nrow(PLANET)+1) | index_4 == (nrow(PLANET)+1)){
        next
    }    
    quad_1 <- PLANET[index_1, 1]
    quad_2 <- PLANET[index_2, 1]
    quad_3 <- PLANET[index_3, 1]
    quad_4 <- PLANET[index_4, 1]
    quads <- list(quad_1, quad_2, quad_3, quad_4)
    observation <- c(metadata[i,3], yr_mo, as.list(quads))
    run_keys <- run_keys %>% rbind(observation)
}

colnames(run_keys) <- c('run_id', 'date', 'file_name_1', 'file_name_2', 'file_name_3', 'file_name_4')

write.csv(run_keys, here("Data","Intermediate","run_keys_2.csv"))






lst <- lapply(1:nrow(metadata), function(x){
  ## create a matrix of coordinates that also 'close' the polygon
  res <- matrix(c(metadata[x, 'Min.Longitude'], metadata[x, 'Min.Latitude'],
           metadata[x, 'Min.Longitude'], metadata[x, 'Max.Latitude'],
           metadata[x, 'Max.Longitude'], metadata[x, 'Max.Latitude'],
           metadata[x, 'Max.Longitude'], metadata[x, 'Min.Latitude'],
           metadata[x, 'Min.Longitude'], metadata[x, 'Min.Latitude'])  ## need to close the polygon
         , ncol = 2, byrow = T)
  ## create polygon objects
  st_polygon(list(res))
})

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

get_date <- sapply(1:nrow(metadata), function(x){
    yr_mo <- substr(metadata[x,4], 1, 7)
})

HyTES_sf <- st_sf(Run_ID = metadata[,'Run_ID'], date = get_date, st_sfc(lst))
PLANET_sf <- st_sf(file_name = PLANET[,'file_name'], date = PLANET[,'date'], st_sfc(planet_lst))


file_names <- list()

for (i in 1:nrow(HyTES_sf)){
    PLANET_files <- c()
    for(j in 1:nrow(PLANET_sf)){
        if(HyTES_sf$date[i] == PLANET_sf$date[j]){
            if (st_contains(st_geometry(PLANET_sf[j,]), st_geometry(HyTES_sf[i,]), sparse = FALSE) == TRUE){
                PLANET_files <- c(PLANET_sf$file_name[j])
                break
            }            
            if (st_overlaps(st_geometry(HyTES_sf[i,]), st_geometry(PLANET_sf[j,]), sparse = FALSE) == TRUE){
                PLANET_files <- c(PLANET_files, PLANET_sf$file_name[j])
            }
        }
    }
    print(PLANET_files)
    if (length(PLANET_files) == 0){
        file_names[[i]] <- NA
    }
    else{
        file_names[[i]] <- PLANET_files
    }
}

new_run_keys <- data.frame('run_id' = metadata[,'Run_ID'], 'date' = get_date, 'file_names' = I(file_names))

write.csv2(new_run_keys, here("Data","Intermediate","run_keys_2.csv"))
