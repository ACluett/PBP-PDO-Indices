##### PBP-PDO-Indices ######
Supplementary code for Cluett et al. (2025).

This README file describes requirements and instructions for running 'ERSST-EOFs-calculation.sh' and 'ERSST-EOFs-plot.R'
Prepared by AA Cluett, 10/01/2025

-------------------------------------------------------------
----------- PART 1: EOF CALCULATIONS WITH CDO ---------------
-------------------------------------------------------------

The script 'ERSST-EOFs-calculation.sh' calculates the main manuscript results. Output is stored in the 'output' folder.

---------SYSTEM REQUIREMENTS ---------

This zsh script was developed using climate data operators (CDO) version 2.4.0 and netcdf operators (NCO) version 5.2.4 on macOS Monterey Version 12.0.

wget is required to download the NOAA Extended Reconstructed SST V5 using the script, although this file can be downloaded directly from: https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v5

This code has been tested on the following versions:
Climate Data Operators version 2.4.0 (https://mpimet.mpg.de/cdo)
System: aarch64-apple-darwin21.6.0
CXX Compiler: llvm_clang++ -std=gnu++20 -g -O2  -pthread
CXX version : Homebrew clang version 17.0.6
C Compiler: llvm_clang -g -O2  -pthread -pthread
C version : Homebrew clang version 17.0.6
F77 Compiler: gfortran -g -O2
F77 version : GNU Fortran (Homebrew GCC 13.2.0) 13.2.0
Features: 16GB 10threads c++20 Fortran pthreads HDF5 NC4/HDF5 OPeNDAP sz
Libraries: yac/3.1.0 NetCDF/4.9.2 HDF5/1.14.3
CDI data types: SizeType=size_t
CDI file types: srv ext ieg grb1 grb2 nc1 nc2 nc4 nc4c nc5 nczarr 
     CDI library version : 2.4.0
 cgribex library version : 2.2.0
 ecCodes library version : 2.35.0
  NetCDF library version : 4.9.2 of May  1 2023 09:07:13 $
    exse library version : 1.5.0
    FILE library version : 1.9.1

NCO netCDF Operators version 5.2.4 "Kamehameha I" built by brew on Monterey-arm64.local at Jun 13 2024 01:37:37
ncks version 5.2.4

No non-standard hardware required.

---------INSTALLATION GUIDE---------

To install CDO and NCO dependencies:
-For CDO, see: https://code.mpimet.mpg.de/projects/cdo/wiki 
-For NCO, see: https://nco.sourceforge.net/

To install wget, see: https://www.gnu.org/software/wget/

All three dependencies can be easily installed using homebrew: https://formulae.brew.sh/

Expected installation time: <5 min

---------DEMO---------

To run:

Execute the script in zsh command line: ./ERSST-EOFs-calculation.sh

User may have to adjust permissions using the following command: chmod +ux ERSST-EOFs-calculation.sh
User may have to remove default apple quarantine using the follow command: xattr -d com.apple.quarantine ERSST-EOFs-calculation.sh

Expected output:

This script will download the most recent NOAA Extended Reconstructed SST V5 file from the NOAA Physical Sciences Laboratory.

EOFs calculated for the 'total' SSTa (i.e., GMSST not removed) (eofs-total-'years'.nc) and GMSST-removed SSTa (eofs-gmsstr-'years'.nc) are produced for retrospective 30-year sliding windows from 1994-2023. The netcdf files can be viewed using Panoply (https://www.giss.nasa.gov/tools/panoply/).

The eigenvalues and percent variance for each are saved in csv files (e.g., eigval-total-'years'.nc, eigval-gmsstr-'years'.nc, var-total-'years'.nc, and var-gmsstr-'years'.nc).

Principal component timeseries for the first two EOFs (GMSST-removed and total SSTa) are exported in csv files (gmsstr-1994-2023-pc1.csv, gmsstr-1994-2023-pc2.csv, total-1994-2023-pc1.csv, and total-1994-2023-pc2.csv).

Expected run time: <5 min

---------INSTRUCTIONS FOR USE---------

This script reproduces the EOFs presented in this study. To reproduce the EOFs calculated on sliding windows of different lengths, change the year defined in 'yearstart' to the appropriate window.


-------------------------------------------------------------
---------------- PART 2: PLOT RESULTS IN R ------------------
-------------------------------------------------------------

The script 'ERSST-EOFs-plot.R' produces the main manuscript figures and exports the PBP values. Figures are stored in the 'figs' folder.


---------SYSTEM REQUIREMENTS ---------

This R script was developed using R version 4.3.1 on macOS Monterey Version 12.0.

This code requires the following R packages and has been tested with the listed versions:
data.table version 1.14.8
dplyr version 1.1.3
ggplot2 version 3.4.4
ggh4x version 0.2.8
lubridate version 1.9.2
metR version 0.14.1
ncdf4 version 1.21
patchwork version 1.1.3
purrr version 1.0.2
raster version 3.6.26
stringr version 1.5.0
tictoc version 1.2
tidyr version 1.3.0
zoo version 1.8.12

No non-standard hardware required.

The following files are required and included in the 'mapdata' folder:
- gridweights.nc
- LME_mask_ERSST.nc

---------INSTALLATION GUIDE---------

For instructions for installing R and R studio, see: https://rstudio-education.github.io/hopr/starting.html

All required packages can be installed from CRAN using the following command: install.packages(c("data.table", "dplyr", "ggplot2", "ggh4x", "lubridate", "metR", "ncdf4", "patchwork", "purrr", "raster", "stringr", " "tictoc",  "tidyr", "zoo"))

---------DEMO---------

To run:

Set the working direcotry to the repository folder ('PBP-PDO-Indices'). Execute the script in R studio.

Expected output:

This script will generate main text figures 1-4. Note that final figure formatting was performed in Adobe Illustrator.
PBP values are also exported in a csv (PBP_Index_Values.csv).

Expected run time: <5 min


