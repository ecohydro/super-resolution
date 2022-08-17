### This file moves the Raw HyTES data folders to the Intermediate HyTES data folders and sorts them based on run

library(here)

# Reads in part of segment_metadata.csv file to get keys between each run and segment
segment_keys <- read.csv(here("Data","segment_metadata.csv"))[2:3]

from <- here('Data','Raw','HyTES')                     #Current path of your folder
to   <- here('Data','Intermediate','HyTES')            #Path you want to move it.

# Transfers both L1 and L2 .hdf5 to intermediate folder grouped by run for further data processing
for (i in 1:nrow(segment_keys)){
    path1 <- paste(from, segment_keys[i,1], sep= '/')
    path2 <- paste(to, segment_keys[i,2], sep = '/')
    files_to_move <- dir(path1, "*.hdf5", ignore.case = TRUE, all.files = TRUE)
    more_files_to_move <- dir(path1, "geo(Alt|Lat|Long).csv$", ignore.case = TRUE, all.files = TRUE)
    dir.create(path2)
    file.copy(from = file.path(path1, files_to_move), to = path2,  recursive = TRUE)
    file.copy(from = file.path(path1, more_files_to_move), to = path2,  recursive = TRUE)
    unlink(path1, recursive = TRUE)
}