# stormwater functions

Scripts and functions for stormwater project processing

merge_years.py
 - Helper script to merge wrf month variable files into yearly files
 - inpu t: -g to specify gcm [required]
 - input: -i to change input data directory [default /home/disk/rocinante/DATA/temp/WRF/var/]
 - input: -v to specify input variable name [default PREC]
 - input: -l to specify output variable name [default PREC_year]
 - note: does not change variable name within filename

mk_csv.py
 - Helper script to extract csv time-series data from netcdf variable file -> creates rawWRF data
 - input: -g to specify gcm [required]
 - input: -v to specify variable name [default PREC]
 - input: -i to specify input directory [default /home/disk/rocinante/DATA/temp/WRF/var/]
 - input: -o to specify output directory [default /home/disk/tsuga2/jswon11/workdir/2019_04_stormwater-newwrf/data/pub/]


Process_stormwater_stats.py
 - Run script to process time-series data and create metrics file
 - input: -l to specify a regex qualifier on locations to run [default all locations within input directory]
 - input: -g to speicfy a regex qualifier on gcms to run [default all gcms within input directory]
 - input: -i to sepeicfy input directory
 - note: does a line count to search for broken files in WRF inputs [1138799 line records]


mk_tableau.py
 - Create tableau files from the metric files within the pub directory

bc_nowin.py
 - Module to create training data and apply bias-correction on a time-series using the training data
 - create_training: obs timeseries [required], his timeseries [required]
 - apply_bc: timeseries [required], training [required]
 

merge_var.py
 - merged numerous other variables onto the precip file



