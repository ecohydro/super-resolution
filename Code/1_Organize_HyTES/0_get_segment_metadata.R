### This file extracts the metadata of each segment and writes it to a .csv file (segment_metadata.csv)

library(here)
library(dplyr)
library(readxl)
library(stringr)
library(tidyr)

# This function processes the file given filepath
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

# Opens and reads the text file with url links
# REPLACE WITH CORRECT DIRECTORY here("Data","text.txt")
links = processFile(here("Data","text.txt"))

# Creates empty data frame to store metadata
metadata <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(metadata) <- c('Folder_Name','Year','Month','Day','Hour','Minute','Second','County', 'State')
index <- 1

# Extracts relevant metadata information and appends it to data frame for 2017-2022 data
for (link in links){
  new_folder <- substr(link, 44, nchar(link) - 4)
  folder_dir <- file.path(here("Data","Raw","HyTES"), new_folder)
  L1_HDF5_file <- list.files(folder_dir)[3]
  if (grepl(".hdf5", L1_HDF5_file)){
    Year <- substr(L1_HDF5_file, 1, 4)
    Month <- substr(L1_HDF5_file, 5, 6)
    Day <- substr(L1_HDF5_file, 7, 8)
    Hour <- substr(L1_HDF5_file, 10, 11)
    Minute <- substr(L1_HDF5_file, 12, 13)
    Second <- substr(L1_HDF5_file, 14,15)
    Location <- gsub("([a-z])([A-Z])","\\1 \\2",sub(".*_","",sub("_L1.*","", L1_HDF5_file)))
    State <- sub(".* ","", Location)
    County <- sub(" [^ ]*$", "", Location)
  }
  else{
    L1_file_name <- list.files(folder_dir)[4]
    Year <- substr(L1_file_name, 1, 4)
    Month <- substr(L1_file_name, 6, 7)
    Day <- substr(L1_file_name, 9, 10)
    Hour <- substr(L1_file_name, 12, 13)
    Minute <- substr(L1_file_name, 14, 15)
    Second <- substr(L1_file_name, 16, 17)
    County <- str_sub(gsub("([a-z])([A-Z])","\\1 \\2",sub(".*_","",sub(".Line.*","", L1_file_name))), 19, -1)
    State <- NA
  }
  metadata_row <- c(new_folder, Year, Month, Day, Hour, Minute, Second, County, State)
  metadata[index,] <- metadata_row
  index <- index + 1
}

# Reads master sheet of metadata provided by Gerardo Rivera
HyTES_excel <- read_excel(here("Data",'Raw',"hytes_anna_boser_UCSB_22Jun2022.xlsx"))

# Filters data to remove ER2 flights and also create the correct time columns
HyTES_metadata <- HyTES_excel %>% filter(Aircraft != 'ER2') %>% 
                  extract("Acquistion DateTime", c("Date","Time"), "([^]+) ([^)]+) ([^]+) ([^)]+)$") %>% 
                  separate(Date, c("Year","Month","Day")) %>% 
                  separate(Time, c("Hour","Minute","Second")) %>%
                  subset(select = c(Year:Second, Run, Aircraft))

# There are a total of 1542 segments
# Extract out HyTES location and append to HyTES metadata data frame
HyTES_metadata_location <- HyTES_excel %>% filter(Aircraft != 'ER2') %>% subset(select = c('State/Country'))
HyTES_metadata <- HyTES_metadata %>% cbind(HyTES_metadata_location)

# Merges data frames together
metadata_merged <- metadata %>% merge(HyTES_metadata, by = c("Year","Month","Day","Hour","Minute","Second"))

# If State is NA, replace from HyTES_metadata master sheet and remove last column
metadata_merged$State[is.na(metadata_merged$State)] <- metadata_merged[,ncol(metadata_merged)]
metadata_merged <- metadata_merged[, -ncol(metadata_merged)]

# Gets run name ID
run_name <- data.frame(apply(metadata_merged[,c(1:3,8:10)], 2, paste0))
run_name <- paste0(run_name$Year, run_name$Month, run_name$Day, run_name$County, run_name$State, run_name$Run) 
run_name <- gsub(" ","", run_name)

# Appends run name ID to metadata data frame and reorders columns
metadata_final <- metadata_merged %>% cbind(run_name) %>%
                  subset(select = -c(10)) %>%
                  rename(Run = run_name) %>%
                  select(7,11,1:6,8:10)

# Writes final segment metadata sheet to a .csv file
write.csv(metadata_final, here("Data","segment_metadata.csv"))
