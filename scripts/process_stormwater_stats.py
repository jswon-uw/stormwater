#!/usr/bin/env python
import stormwater_functions as sf
import pandas as pd
from glob import glob
import os.path as osp
import sys



gcms = [
    'access1.0', 
    'access1.3',
    'bcc-csm1.1',
    'canesm2',
    'ccsm4',
    'csiro-mk3.6.0',
    'fgoals-g2',
    'gfdl-cm3',
    'giss-e2-h',
    'miroc5',
    'mri-cgcm3',
    'noresm1-m'
]

gcms = ['ccsm4', 'canesm2']
srcs = ['OBS', 'rawWRF', 'bcWRF']
srcs = ['rawWRF']

idir = '/home/disk/tsuga2/jswon11/workdir/2019_04_stormwater-newwrf/data/'


for s in srcs:
    data = glob('{}/pub/*/{}/*5.csv'.format(idir, s))
    for d in data:
        basename = osp.splitext(d)[0]
        print(osp.basename(basename))
        
        
        # Filtered tables
        df = sf.read_ts(d)
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
        mo_df = sf.calc_moagg(ra_df)
        com = 'Raw WRF Monthly maximum precipitation values in mm'
        out = '{}.MOmax.csv'.format(basename)
        sf.csv_comment(mo_df, out, com)
        print('MO MAX')
        
        
        # Calculate Water Year Maximum
        wy_df = sf.calc_wyagg(ra_df)
        com = 'Raw WRF Water Year maximum precipitation values in mm'
        out = '{}.WYmax.csv'.format(basename)
        sf.csv_comment(wy_df, out, com)
        print('WY MAX')
        
        
        # Calculate Seasonal Maximum
        sn_df = sf.calc_snagg(ra_df)
        com = 'Raw WRF Seasonal maximum precipitation values in mm'
        out = '{}.SNmax.csv'.format(basename)
        sf.csv_comment(sn_df, out, com)
        print('SN MAX')
        
        
        # Calculate statistics
        stats = sf.calc_stats(ra_df, wy_df, sn_df, mo_df, False)
        com = 'Raw WRF precipitation statistics (mm)'
        out = '{}.stats-abs_vals.csv'.format(basename)
        sf.csv_comment(stats, out, com)
        print('STATS')
        

        if s != 'OBS':
            # Calculate percent change statistics
            pchg = sf.calc_pchg(stats)
            com = 'Percent change in Raw WRF Precipitation statistics'
            out = '{}.stats-pct_chg.csv'.format(basename)
            sf.csv_comment(pchg, out, com)
            print('PCHG')    
