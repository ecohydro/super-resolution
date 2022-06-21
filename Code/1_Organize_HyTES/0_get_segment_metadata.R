# This file extracts the metadata of each segment!

library(here)
library(rhdf5)

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

links = processFile(here("Data","text.txt"))

time_types = list('UTC_Year', 'UTC_Month', 'UTC_Day','UTC_Hour','UTC_Minute')
time_metadata <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(time_metadata) <- c('Run','Segment','Year','Month','Day','Hour','Minute')
index <- 1

for (link in links){
  run <- "HELP!"
  new_folder <- substr(link, 44, nchar(link) - 4)
  folder_dir <- file.path(here("Data","Raw","HyTES"), new_folder)
  L1_HDF5_file <- list.files(folder_dir)[3]
  # This is how to extract the L2 file: (for later!)
  # L2_HDF5_file <- list.files(folder_dir)[10]
  time <- c(run, new_folder)
  for (t in time_types){
    time_path <- paste0('locational_metadata/', t)
    recorded_median <- median(h5read(here(folder_dir, L1_HDF5_file), time_path))
    time <- c(time, recorded_median)
  }
  time_metadata[index,] <- time
  index <- index + 1
}
