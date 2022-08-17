### This file downloads each .zip file and saves it to the Raw HyTES folder.
library(here)

# This function unzips the file and places it in a directory (exdir)
unzip <- function(file, exdir, .file_cache = FALSE) { #the normal unzip doesn't work well with large files. 
  
  if (.file_cache == TRUE) {
    print("decompression skipped")
  } else {
    
    # Set working directory for decompression
    # simplifies unzip directory location behavior
    wd <- getwd()
    setwd(exdir)
    
    # Run decompression
    decompression <-
      system2("unzip",
              args = c("-o", # include override flag
                       file),
              stdout = TRUE)
    
    # uncomment to delete archive once decompressed
    # file.remove(file) 
    
    # Reset working directory
    setwd(wd); rm(wd)
    
    # Test for success criteria
    # change the search depending on 
    # your implementation
    if (grepl("Warning message", tail(decompression, 1))) {
      print(decompression)
    }
  }
}   

# This funbction downloads the extracted data
download_extract <- function(url, directory, new_folder){
  dir.create(here(directory, new_folder))
  destination <- here(directory, new_folder, paste0(new_folder, ".zip"))
  download.file(url, destination)
  unzip(destination, exdir = here(directory, new_folder))
  unlink(destination)
}

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

options(timeout = max(30000, getOption("timeout"))) #timeout for download.file is 60s by default; not enough for large files

################################################################################
# Download all url links
################################################################################

links = processFile(here('Data','text.txt'))

for (link in links){
  new_folder <- substr(link, 44, nchar(link) - 4)
  download_extract(url = link,
                   directory = here("Data","Raw","HyTES"),
                   new_folder = new_folder)
}
