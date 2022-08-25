# Super Resolution Project
This repository is for processing JPL's HyTES data and PLANET's RGB data to create finalized raster datasets to be then fed into a model which will be able to heighten the resolution of LST (land surface temperature) images from 70m to 20m.
See data from JPL's website here: https://hytes.jpl.nasa.gov/


## Structure and Contents

This repository is organized into two folders: Data and Code. The Data folder contains raw and intermediate datasets, as well as final dataset used in future model processing. The code folder consists of code to (1) download data and (2) process the data into the intermediate and final datasets. The code folder is numbered in order of folders that are completed sequentially, and the numbers listed next to the data folders and data files correspond to the number of the code file used to create it.

* Code: All code for data download and processing

    * 0_Download_HyTES: Download all raw HyTES data. Note: Each individual HyTES product is sent via email, so these following files makes the downloading process more streamline.
        * `0_extract_emails.scpt`: An Apple Script used to extract all email content sent for each HyTES product.
        * `1_extract_links.py`: This file extracts the email contents from the apple script and writes the zip file links to a text file (`Data/Raw/links.txt`).
        * `2_download_zip.R`: This file downloads each .zip file and saves it to the Raw HyTES folder.
    * 1_Organize_HyTES: These files help organize each HyTES product by run.
        * `0_get_segment_metadata.R`: This file extracts the metadata of each segment and writes it to a `.csv` file located in the `Segment_Metadata` folder.
        * `1_organize_by_run.R`: This file moves the Raw HyTES data folders to the Intermediate HyTES data folders and sorts them based on run.
    * 2_Get_Metadata.R: This folder retrieves metadata information for each *run*
        * `0_run_metadata.R`: This file extracts all relevant metadata information for each run and writes information to `Run_Metadata` folder.
    * 3_Shapefile_HyTES: Creates shape files for each HyTES run.
        * `0_shapefile_HyTES.R`: This file converts the Intermediate HyTES `.hdf5` and `.csv` files into shape files and saves them in a separate folder. Also calculates surface_radiance for each pixel.
    * 4_Sort_PLANET: Sorts PLANET for each HyTES run
        * `0_assign_PLANET.R`: This file finds which PLANET files to use for each run and saves the keys to the Intermediate folder.
    * 5_Resample_HyTES: Resample HyTES data to ~5m PLANET resolution after cropping PLANET data to HyTES data.
        * `0_stitch_PLANET.R`: Stitch relevant PLANET files, crop the PLANET data to the extent of each run, and convert resampled HyTES data to rasters for model.
        * `1_big_runs.R`: For runs with excessive amounts of associated PLANET files, make an individual raster for each PLANET file.
    * 6_Make_Tiles: Makes tiles of 672 x 672 for each HyTES raster.
        * `0_tile.R`: Tile images of an arbitrary size into equal size tiles of 672 x 672.
        * `1_input_output_resample.R`: Take the target and basemap files and aggregate and resample them to the input and output size.
    * Scratch_Code: A folder containing scratch code that are either used for debugging or were simply not used for other Code files.


* Data (1,0): All data is listed in the `.gitignore` and is thus not stored on github. Note: All HyTES raw data has been removed due to space shortage. To retrieve HyTES data in its raw form, please order L2 products from https://hytes.jpl.nasa.gov/order. If operating on Mac, enable Rules on your Mail app and use the AppleScript and python script from `Code/0_Download_HyTES/0_extract_emails.scpt` and `Code/0_Download_HyTES/1_extract_links.py` respectively to extract all zip file url links.

    * Raw: All raw data, which can be found at given urls and/or retrieved using scripts in `code/0_Download_HyTES`. 
        * HyTES (0,2): This folder contains all product data from HyTES for each segment from 2014 onwards. Specifically, each segment folder contains geopixel information (L1 data) and LST/Emissivity data (L2 data). 
        * PLANET: This folder contains monthly composite RGB data from PLANET. Specifically, each PLANET `.tiff` file contain 4-bands for 4096 x 4096 pixels at 4.77 m GSD. 
        * `links.txt` (0,1): This `.txt` file has a list of the zip file links for each HyTES product.
        * `metadata.xlsx`: Finally, the Raw Folder also contains a `.xlsx` file which provides all metadata information for each individual HyTES product/segment.
    * Intermediate: All intermediate data including shape file versions of HyTES data
        * HyTES (1,1): HyTES data grouped by run. Also only contains relevant L1 and L2 data.
        * HyTES_sf (3,0): HyTES data in the form of a shapefile. Also contains surface radiance for each pixel.
        * `file_error.txt` (3,0): Contains run IDs for all files which failed to open correctly.
        * `run_keys.csv` (4,0): Data frame where each observation is each unique run ID followed by run date (year-month) and a list of all associated PLANET files.
        * `wrong_dim.txt` (3,0): A list of run IDs where the dimension sizes of the L1 and L2 data do not match correctly by a significant amount. 
    * Final: Contains finalized HyTES data in the form of a raster.
        * Emis (5,0/1): Contains Emissivity data for each HyTES raster.
        * LST (5,0/1): Contains LST data for each HyTES raster.
        * RGB (5,0/1): Contains RGB data for each HyTES raster.
        * Srfc_Rad (5,0/1): Contains Surface Radiance data for each HyTES raster.
        * `complete_run.txt` (5,0): A list of all completed HyTES runs for data processing.
        * `fat_runs.txt` (5,0): A list of HyTES runs which have more than 4 PLANET files associated to the run. 
    * For_CNN: 
    * Run_Metadata: Contains all metadata information for each HyTES *run*. Also contains a shape file version of all run metadata information.
    * Segment_Metadata: Contains all metadata information for each HyTES *segment* (multiple files since the process was done in batches)
    * Tiles:
