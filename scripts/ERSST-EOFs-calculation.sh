#!/bin/zsh

export CDO_WEIGHT_MODE=on
export MAX_JACOBI_ITER=100

mkdir -p -- 'code-example'

wget "https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v5/sst.mnmean.nc"
cdo -selgrid,2 sst.mnmean.nc grid.nc # select correct lat/lon grid

# Calculate anomalies
cdo -L -ymonsub grid.nc -ymonmean grid.nc total-anoms.nc # calculate monthly anoms
cdo fldmean total-anoms.nc ersst_globalmean.nc # calculate global mean ssta
cdo -L sub ersst_globalmean.nc -timmean ersst_globalmean.nc ersst_gAnom.nc # calculate gmsst anom
cdo enlarge,total-anoms.nc ersst_gAnom.nc tmpfile.nc # expand gmssta to full grid
ncdiff total-anoms.nc tmpfile.nc gmsstr-anoms.nc # subtract gmssta from each grid cell
ncks -C -H -v sst -s "%.10f\n" ersst_gAnom.nc > ersst_globalmean.csv # export gmssta

# Crop to PDO region
cdo sellonlatbox,105,260,20,65 total-anoms.nc total-crop.nc # crop total anoms to PDO region
cdo sellonlatbox,105,260,20,65 gmsstr-anoms.nc gmsstr-crop.nc # crop gmsstr anoms to PDO region

j=1 # Begin to loop through 30 year sliding windows
while [ $j -le 30 ]
do

	yearstart=$((1994 - ($j-1) )); echo $chunkstart
	yearend=$((2023- ($j-1))); echo $chunkend

	# total anoms
	cdo selyear,$yearstart/$yearend total-crop.nc total-crop-"$yearstart"_"$yearend".nc #subset years
	cdo -v -L eof,3 total-crop-"$yearstart"_"$yearend".nc eigval-total-"$yearstart"_"$yearend".nc eofs-total-"$yearstart"_"$yearend".nc #calculate eofs
	cdo -L div eigval-total-"$yearstart"_"$yearend".nc -timsum eigval-total-"$yearstart"_"$yearend".nc var-total-"$yearstart"_"$yearend".nc # calculate %var
	ncks -C -H -v sst -s "%.10f\n" var-total-"$yearstart"_"$yearend".nc > var-total-"$yearstart"_"$yearend".csv # export %var
	ncks -C -H -v sst -s "%.10f\n" eigval-total-"$yearstart"_"$yearend".nc > eigval-total-"$yearstart"_"$yearend".csv # export eigs

	# repeat for gmsstr anoms
	cdo selyear,$yearstart/$yearend gmsstr-crop.nc gmsstr-crop-"$yearstart"_"$yearend".nc
	cdo -v -L eof,3 gmsstr-crop-"$yearstart"_"$yearend".nc eigval-gmsstr-"$yearstart"_"$yearend".nc eofs-gmsstr-"$yearstart"_"$yearend".nc
	cdo -L div eigval-gmsstr-"$yearstart"_"$yearend".nc -timsum eigval-gmsstr-"$yearstart"_"$yearend".nc var-gmsstr-"$yearstart"_"$yearend".nc
	ncks -C -H -v sst -s "%.10f\n" var-gmsstr-"$yearstart"_"$yearend".nc > var-gmsstr-"$yearstart"_"$yearend".csv
	ncks -C -H -v sst -s "%.10f\n" eigval-gmsstr-"$yearstart"_"$yearend".nc > eigval-gmsstr-"$yearstart"_"$yearend".csv

    ((j++))
done

# Project 2023-1994 EOFs on full timeseries to produce indices

cdo eofcoeff eofs-gmsstr-1994_2023.nc gmsstr-crop.nc gmsstr-1994-2023-ts.nc
ncks -C -H -v sst -s "%.10f\n" gmsstr-1994-2023-ts.nc00000.nc > gmsstr-1994-2023-pc1.csv
ncks -C -H -v sst -s "%.10f\n" gmsstr-1994-2023-ts.nc00001.nc > gmsstr-1994-2023-pc2.csv

cdo eofcoeff eofs-total-1994_2023.nc total-crop.nc total-1994-2023-ts.nc
ncks -C -H -v sst -s "%.10f\n" total-1994-2023-ts.nc00000.nc > total-1994-2023-pc1.csv
ncks -C -H -v sst -s "%.10f\n" total-1994-2023-ts.nc00001.nc > total-1994-2023-pc2.csv


