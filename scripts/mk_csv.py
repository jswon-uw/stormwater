#!/usr/bin/python3
import xarray as xr
import pandas as pd
import numpy as np
from glob import glob
import argparse
import os
import sys
from pprint import pprint as print


sys.path.append('/home/disk/tsuga2/jswon11/workdir/2019_05_wrf-netcdf/scripts/')
import wnp

parser = argparse.ArgumentParser()
parser.add_argument('-g', nargs='?', default=None, help='gcm')
parser.add_argument('-v', nargs='?', default='PREC', help='variable directory')
parser.add_argument('-i', nargs='?', default=None, help='input directory')
parser.add_argument('-o', nargs='?', default=None, help='output directory')
args = parser.parse_args()

if args.g:
    gcm = args.g
else:
    print('Provide a gcm')
    sys.exit()
 
bdir = args.i if args.i else '/home/disk/rocinante/DATA/temp/WRF/var/'
odir = args.o if args.o else '/home/disk/tsuga2/jswon11/workdir/2019_04_stormwater-newwrf/data/pub/'
var = args.v
    
# Get list of lat/lon coordinates to extract
def get_list():
    print('Getting lat/lon coordinates')
    mlist = '/home/disk/tsuga2/jswon11/workdir/2019_04_stormwater-newwrf/scripts/MASTER_Precip_Gauge_List.csv'
    #mlist = '/home/disk/tsuga2/jswon11/workdir/2019_04_stormwater-newwrf/scripts/old/short2.csv'
    df = pd.read_csv(mlist)
    df = df[['Source', 'ID', 'stn_lat', 'stn_lon']]
    df.columns = ['source', 'pid', 'lat', 'lon']
    return df

# Get datasets to given lat/lon coordinates from monthly netcdf files
def get_data(gcm, var, mdf, bdir):
    ddir = '{}/{}/2d/{}/*d02*'.format(bdir, gcm, var)
    data = sorted(glob(ddir))
    
    if len(data) <= 0:
        print('No data found.')
        sys.exit()
    
    print('Begin parsing netcdf:')
    #print(data)
    
    ds = wnp.arr_extract_latlon(data, mdf[['lat', 'lon']])
    return ds
    

# Save csv from the dataset for each of the coordinates
def save_csv(gcm, ds_arr, mdf):
    print('Saving csvs:')
        
    for key in ds_arr.keys():
        ds = ds_arr[key]
        lat = '{:0.5f}N'.format(abs(ds.lats.values.tolist()))
        lon = '{:0.5f}W'.format(abs(ds.lons.values.tolist()))
        row = mdf[(mdf['lat'] == key[0]) & (mdf['lon'] == key[1])]
        source = row.source.values.tolist()[0]
        pid = row.pid.values.tolist()[0]
        print('    Saving: {} {} - {}, {}'.format(source, pid, lat, lon))
        
        xdir = '{}/{}_{}/rawWRF/'.format(odir, source, pid)
        if not os.path.exists(xdir):
            os.makedirs(xdir)
 
        outpath = '{}/{}_{}_{}_{}_rawWRF_{}_rcp85.csv'.format(xdir, source, pid, lat, lon, gcm)

        print(outpath)
        wnp.ds2csv(ds, outpath, ['Precip (mm)'])
        



def main():        
    mdf = get_list()
    ds_arr = get_data(gcm, var, mdf, bdir)
    save_csv(gcm, ds_arr, mdf)
    
    print('All clear')
main()


## Format
## YYYY, MM, DD, HH, PRECIP(mm)(round 5)
## Naming
## [SOURCE]_[ID]_[WRF_LAT]N_[WRF_LON]W_[rawWRF]_[GCM]_[RCP].csv


#T2 - T 2m (C)
#Q2 - Specific Humidty (kg/kg)
#SWDOWN - Incoming Shortwave (W/m2)
#GLW - Incoming Longwave (W/m2)
#U10/V10 - Wind Speed(m/s)
#PREC - Precip (mm)
#PSFC - Pressure (Pa)
#Tdew - Tdew 2m (C)
#Cloud - Cloud Fraction (-)
#pet - Pot. Evap. (mm/day)
