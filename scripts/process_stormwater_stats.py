#!/usr/bin/env python
import stormwater_functions as sf
import pandas as pd
from glob import glob
import os.path as osp
from datetime import datetime 
import sys
import argparse
#from pprint import pprint as print

#srcs = ['OBS', 'rawWRF', 'bcWRF']
#srcs = ['rawWRF']
srcs = ['bcWRF']
#srcs = ['OBS']
#srcs = ['rawWRF', 'bcWRF']


parser = argparse.ArgumentParser()
parser.add_argument('-l', default='', help='location specifier')
parser.add_argument('-g', default='', help='gcm specifier')
parser.add_argument('-i', nargs='?', default=None, help='input directory')
parser.add_argument('-t', default=False, help='toggle logs')
parser.add_argument('--stats', action='store_true', help='Use existing max files to remake stats')
parser.add_argument('-m', type=int, default=-1, help='Min year threshold for stat calculations')
args = parser.parse_args()
lreg = args.l
greg = args.g
idir = args.i if args.i else '/home/disk/tsuga2/jswon11/workdir/2019_04_stormwater-newwrf/data/pub/'
myr = args.m
wc = 1138799
lname = datetime.now().strftime('logs/process-sf_%Y-%m-%d_%H-%M.log')
if args.t:
    logf = open(lname, 'w')


# Helper function to run aggregated max functions
def run_max(f, ra, out, com, msg):
    if osp.exists(out) & (args.stats):
        xdf = pd.read_csv(out, skiprows=[0])        
        print('{}: {}'.format(msg, out))        
    else:
        xdf = f(ra)
        sf.csv_comment(xdf, out, com)
        print(msg)

    return xdf
    

for s in srcs:
    data = glob('{}/*/{}/*{}*{}*_rcp85.csv'.format(idir, s, lreg, greg))
    data = sorted(data + glob('{}/*/{}/*{}*{}*W.csv'.format(idir, s, lreg, greg)))
    print('{}/*/{}/*{}*{}_rcp??.csv'.format(idir, s, lreg, greg))
    
    for d in data:
        basename = osp.splitext(d)[0]
        print(osp.basename(basename))
                
        # Filtered tables
        df = sf.read_ts(d)
        if (s != 'OBS') & (len(df) < wc):
            print('File appears to be broken: {}'.format(osp.basename(d)))
            if args.t:
                logf.write('BROKEN: {}\n'.format(osp.basename(d)))
            continue
        
        nz = sf.filter_nonzero(df)
        nz.to_csv('{}.non-zero.csv'.format(basename), index=False)

        if s != 'OBS':
            his = sf.filter_years(df, 1970, 2005)
            fut = sf.filter_years(df, 2006, 2099)
            his.to_csv('{}.1970-2005.csv'.format(basename), index=False)
            fut.to_csv('{}.2006-2099.csv'.format(basename), index=False)
            
                
        # Calculate Running Average
        ra_df = sf.calc_runsum(df)
        print('RUN SUM')

        
        # Calculate Monthly Maximum
        out = '{}.MOmax.csv'.format(basename)
        com = 'Raw WRF Monthly maximum precipitation values in mm\n'
        msg = 'MO MAX'
        mo_df = run_max(sf.calc_moagg, ra_df, out, com, msg)

        # Calculate Water Year Maximum
        out = '{}.WYmax.csv'.format(basename)
        com = 'Raw WRF Water Year maximum precipitation values in mm\n'
        msg = 'WY MAX'
        wy_df = run_max(sf.calc_wyagg, ra_df, out, com, msg)
        
        # Calculate Seasonal Maximum
        out = '{}.SNmax.csv'.format(basename)
        com = 'Raw WRF Seasonal maximum precipitation values in mm\n'
        msg = 'SN MAX'
        sn_df = run_max(sf.calc_snagg, ra_df, out, com, msg)

        # Calculate statistics
        stats = sf.calc_stats(ra_df, wy_df, sn_df, mo_df, myr, s=='OBS')        
        com = 'Raw WRF precipitation statistics (mm)\n'
        ftype = 'stats' if s == 'OBS' else 'stats-abs_vals'
        out = '{}.{}.csv'.format(basename, ftype)
        sf.csv_comment(stats, out, com)
        print('STATS')
        

        if s != 'OBS':
            # Calculate percent change statistics
            pchg = sf.calc_pchg(stats)
            com = 'Percent change in Raw WRF Precipitation statistics\n'
            out = '{}.stats-pct_chg.csv'.format(basename)
            sf.csv_comment(pchg, out, com)
            print('PCHG')    

if args.t:    
    logf.close()

