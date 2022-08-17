library(here)
library(hdf5r)

# Debug for manually moving files over
files <- list.files(here("Data","Raw","HyTES"))

text_file <- file(here("Code","Scratch_Code", "text.txt"))

for (file in files){
    link <- paste0("https://hytes.jpl.nasa.gov/orders_complete/", file, ".zip")
    write(link, file = here("Code","Scratch_Code", "text.txt"), append = TRUE) 
    print(file)
}

close(text_file)



# Debug for old data
processFile = function(filepath) {
  lines <- list()
  con = file(filepath, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    lines <- append(lines, line)
  }
  close(con)
  return(lines)
}


links <- processFile(here("Code","Scratch_Code","text.txt"))

# CHANGE TO TEST DIFF LINKS
link <- links[4]
new_folder <- substr(link, 44, nchar(link)-4)
folder_dir <- file.path(here("Data","Raw", "HyTES", new_folder))
list.files(folder_dir)
test_lat <- read.csv(here(folder_dir,list.files(folder_dir)[4]))
test_lon <- read.csv(here(folder_dir, list.files(folder_dir)[5]))
L2_data <- H5File$new(here(folder_dir, list.files(folder_dir)[14]))

print(dim(test_lat))
print(dim(test_lon))
L2_data



# Adds Polygon object to each run observation in order to create the GeoJSON file
final_metadata <- read.csv(here("Data","Run_Metadata","run_metadata_2.csv"))

first = TRUE
for (i in 1:nrow(final_metadata)){
    lon <- c(final_metadata$Min_Longitude[i], final_metadata$Max_Longitude[i])
    lat <- c(final_metadata$Min_Latitude[i], final_metadata$Max_Latitude[i])
    Poly_Coord_df <- data.frame(lon, lat)
    Polygon <- Poly_Coord_df %>% 
        st_as_sf(coords = c("lon", "lat"), 
         crs = "WGS84") %>% 
        st_bbox() %>% 
        st_as_sfc()
    if (first){
        run_metadata_sf <- st_sf(final_metadata[i,], Polygon)
        first = FALSE
    }
    else{
        new_metadata_sf <- st_sf(final_metadata[i,], Polygon)
        run_metadata_sf <- run_metadata_sf %>% rbind(new_metadata_sf)
   }
}
 
# Write a shapefile of metadata information 
st_write(run_metadata_sf, here("Data", "Run_Metadata", "metadata.shp"))
