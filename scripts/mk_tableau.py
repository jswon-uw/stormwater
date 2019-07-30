#!/usr/bin/python3

import pandas as pd
from glob import glob
import numpy as np
import os
import sys
from pprint import pprint
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-l', default='', help='location specifier')
parser.add_argument('-g', default='', help='gcm specifier')
parser.add_argument('-i', nargs='?', default=None, help='input directory')
parser.add_argument('-o', nargs='?', default=None, help='output directory')

args = parser.parse_args()
lreg = args.l
greg = args.g
idir = args.i if args.i else '/home/disk/tsuga2/jswon11/workdir/2019_04_stormwater-newwrf/data/pub/'
odir = '/home/disk/tsuga2/jswon11/workdir/2019_04_stormwater-newwrf/data/tableau/'
odir = '/home/disk/tsuga2/jswon11/workdir/2019_04_stormwater-newwrf/scripts/minyr/'
odir = args.o if args.o else odir


mstf = '/home/disk/tsuga2/jswon11/workdir/2019_04_stormwater-newwrf/scripts/MASTER_Precip_Gauge_List.csv'
wrfs = [ 'rawWRF', 'bcWRF']
hformat = ['Obs Network', 'Station Name', 'Station ID','lat', 'lon',
           'OBS/rawWRF/bcWRF', 'Model', 'Scenario', 'Statistic', 'Metric',
           'Duration', 'N(obs)', 'Precip (mm)', 'Years']
vformat = ['Obs Network', 'Station Name', 'Station ID', 'lat', 'lon', 'Model',
           'Scenario', 'Statistic', 'Metric', 'Duration', 'N(obs)',
           'rawWRF (pct bias)', 'bcWRF (pct bias)']

fformat = ['Obs Network', 'Station Name', 'Station ID', 'lat', 'lon', 'Model',
           'Scenario', 'Future YRs', 'change_or_timeseries', 'Statistic', 'Metric',
           'Duration', 'rawWRF', 'bcWRF']

oformat = ['Obs Network', 'Station Name', 'Station ID','lat', 'lon',
           'Metric', 'Duration', 'N(obs)']

skey = ['Water Year or Month', 'Statistic', 'Years']

#exlist = ['Water Year', 'DJF', 'MAM', 'JJA', 'SON']
exlist = ['January', 'February', 'March', 'April', 'May', 'June', 'July',
          'August', 'September', 'October', 'November', 'December']
gcmrep = {'access1.0':'ACCESS 1.0','access1.3':'ACCESS 1.3', 'bcc-csm1.1':'BCC-CSM1.1',
          'canesm2':'CANESM2', 'ccsm4':'CCSM4', 'csiro-mk3.6.0':'CSIRO-MK3.6.0',
          'fgoals-g2':'FGOALS-G2', 'giss-e2-h':'GISS-E2-H', 'gfdl-cm3':'GFDL CM3',
          'miroc5':'MIROC5', 'mri-cgcm3':'MRI-CGCM3', 'noresm1-m':'NORESM1-M'}

rh = 'rawWRF'
bh = 'bcWRF'


# Create Nobs file -> prelim headers
def mk_nobs_file(out='nobs.csv'):
    print('Making Nobs file')
    mdf = pd.read_csv(mstf)
    c = mdf.columns 
    morder = ['Water Year', 'DJF', 'MAM', 'JJA','SON']
    flbl = ['Metric', 'Duration', 'N(obs)']
    durs = [1, 2, 3,6,12,24,48,72,120,240,360]    
    df = pd.DataFrame()

    # Filter
    mdf = mdf[(mdf[c[0]] + '_' + mdf[c[1]]).str.contains(lreg)]
    
    for idx,line in mdf.iterrows():
        src = line[0]
        name = line[3].replace('_', ' ')
        sid = line[1]
        lat = '{:0.5f}'.format(line[4])
        lon = '{:0.5f}'.format(line[5])
        
        wyfile = glob('{}/{}_{}/OBS/*.WYmax.csv'.format(idir, src, sid))
        snfile = glob('{}/{}_{}/OBS/*.SNmax.csv'.format(idir, src, sid))

        print('{} {}'.format(src, sid))
        if (len(wyfile) == 0) | (len(snfile) == 0):

            cdf = pd.DataFrame({'Metric':morder * len(durs),
                               'Duration':durs * len(morder),
                               'N(obs)':[np.nan]*len(morder) * len(durs)})
                    
        else:
            wyfile = wyfile[0]
            snfile = snfile[0]
            cdf = pd.DataFrame()
        
            wy_df = pd.read_csv(wyfile, skiprows=[0])
            wy_df.drop(wy_df.columns[0], axis=1, inplace=True)
            wy_df.columns = wy_df.columns.str.replace('-hr Precip', '')
            wy_df = pd.DataFrame(wy_df.count()).reset_index()
            wy_df.columns = ['Duration', 'N(obs)']
            wy_df.Duration = wy_df.Duration.astype(int)
            wy_df.insert(0, 'Metric', 'Water Year')
            cdf = cdf.append(wy_df)
            
            sn_df = pd.read_csv(snfile, skiprows=[0])
            sn_df.drop(sn_df.columns[0], axis=1, inplace=True)
            sn_df.columns = sn_df.columns.str.replace('-hr Precip', '')
            sn_df = pd.melt(sn_df, id_vars=['Season'])
            sn_df.variable = sn_df.variable.astype(int)
            
            sn_g = sn_df.groupby(['Season', 'variable']).count().reset_index()
            sn_g['Season'] = pd.Categorical(sn_g['Season'], morder)
            sn_g = sn_g.sort_values(['Season', 'variable'])
            sn_g.columns = flbl
            cdf = cdf.append(sn_g)

        cdf.insert(0, 'Obs Network', src)
        cdf.insert(1, 'Station Name', name)
        cdf.insert(2, 'Station ID', sid)
        cdf.insert(3, 'lat', lat)
        cdf.insert(4, 'lon', lon)
        df = df.append(cdf)
    df = df[oformat]
    df.to_csv(out, index=False)
    return df
        

# Helper function to format the dataframe for idf
def format_df(df, obs, src, fid, wrf, gcm, rcp, hformat):
    keys = ['Metric', 'Statistic', 'Years']
    cname = keys + [x.replace('-hr Precip', '') for x in df.columns if 'hr' in x]
    df.columns = cname
    df = pd.melt(df, id_vars=keys, var_name='Duration', value_name='Precip (mm)')
    df.Duration = df.Duration.astype(int)
    df['Obs Network'] = src
    df['Station ID'] = fid
    
    df = df.merge(obs, on=['Obs Network', 'Station ID', 'Metric', 'Duration'])
    df['OBS/rawWRF/bcWRF'] = wrf
    df['Model'] = gcm
    df['Scenario'] = rcp
    df = df[hformat]                
    
    print('{}_{} - {} | {}'.format(src, fid, wrf, gcm))
    return df


# Read file, filter out monthly data and melt to long form
def read_file(nfile, numid, lbl, ilist=exlist):
    df = pd.read_csv(nfile, skiprows=[0])
    name1 = df.columns[0]
    df = df[~df[name1].isin(ilist)]
    keys = df.columns[0:numid]
    df.columns = [x.replace('-hr Precip', '') for x in df.columns]
    df = df.melt(id_vars=keys, var_name='Duration', value_name=lbl)
    df.Duration = df.Duration.astype(int)    
    return df


# Read stats file and filter to earliest period only
def his_read(nfile, lbl):
    df = read_file(nfile, 3, lbl)
    df = df[df['Years'] == df['Years'].min()]
    df = df.drop(['Years'], axis=1)            
    df.columns = ['Metric', 'Statistic', 'Duration', lbl]    
    return df

    
def mk_idf_file(odf, out='idf.csv'):
    locs = odf[['Obs Network', 'Station ID']].drop_duplicates()
    #locs = locs[locs['Obs Network'].str.contains(lreg)]
    
    f = open(out, 'w')
    df = pd.DataFrame(columns = hformat)
    df.to_csv(f, index=False)
    
    for idx, loc in locs.iterrows():
        src = loc[0]
        fid = loc[1]
        codf = odf[(odf['Obs Network'] == src) & (odf['Station ID'].astype(str) == fid)]
        
        print('{} {}'.format(src, fid))
        
        # OBS merging
        data = glob('{}/{}_{}/OBS/*.stats.csv'.format(idir, src, fid))
        for d in data:
            s = os.path.basename(d).split('_')
            wrf = 'OBS'
            gcm = np.nan
            rcp = np.nan
            cdf = pd.read_csv(d, skiprows=[0])
            cdf = format_df(cdf, codf, src, fid, wrf, gcm, rcp, hformat)
            cdf.to_csv(f, index=False, header=False)
            
        # WRF merging
        data = glob('{}/{}_{}/*WRF/*rcp85.stats-abs_vals.csv'.format(idir, src, fid))
        ldf = pd.DataFrame()
        
        for wrf in wrfs:
            # Process gcms
            for d in data:
                s = os.path.basename(d).split('_')
                gcm = s[5].upper()
                rcp = s[6].split('.')[0].replace('rcp', 'RCP ').replace('5', '.5')
                
                cdf = pd.read_csv(d, skiprows=[0])
                cdf = format_df(cdf, codf, src, fid, wrf, gcm, rcp, hformat)
                ldf = ldf.append(cdf)
                
                
            # Calculate ensemble average
            print('-'*20)

            if ldf.empty:
                continue
            
            ldf['Precip (mm)'] = ldf['Precip (mm)'].astype(float).round(2)
            grp = [x for x in hformat if x not in ['Precip (mm)', 'Model']]
            gdf = ldf[grp + ['Precip (mm)']].groupby(grp).agg(np.nanmean).reset_index()
            gdf['Model'] = 'Ensemble'
            if gdf.empty:
                continue
            gdf = gdf[hformat]
            
            ldf.to_csv(f, index=False, header=False)
            gdf.to_csv(f, index=False, header=False)
        




# Helper function to read files for future pct changes
def fpc_reader(rfile, n):
    rdf = read_file(rfile, n, rh)
    bfile = rfile.replace(rh, bh)

    # Merge bcWRF
    if os.path.exists(bfile):
        bdf = read_file(bfile, n, bh)
        rdf = rdf.merge(bdf)
    else:
        rdf[bh] = np.nan

    return rdf

# Helper function to calculate percent change
def calc_pct(rdf, hdf):
    c = [x for x in hdf.columns if x not in [rh, bh]]
    hdf = hdf.groupby(c).mean().reset_index()
    hdf.columns = c + ['x', 'y']
    hdf = hdf.drop(['Future YRs'], axis=1)
    
    cdf = rdf.merge(hdf)
    cdf[rh] = ((cdf[rh] / cdf['x'] - 1)  * 100).round(2)
    cdf[bh] = ((cdf[bh] / cdf['y'] - 1)  * 100).round(2)
    cdf.drop(['x', 'y'], axis=1, inplace=True)

    return cdf

# Create future_pct_changes file
def mk_fpc_file(ndf, out):
    cfilter = [oformat[0], oformat[2]]
    
    f = open(out, 'w')
    df = pd.DataFrame(columns = fformat)
    df.to_csv(f, index=False)

    locs = ndf[~(ndf[cfilter].duplicated())].sort_values(by=cfilter)

    for idx, line in locs.iterrows():
        src = line[0]
        sid = line[2]
        print('{} {}'.format(src, sid))
        cdf = pd.DataFrame()

        # Get stats file
        rstats = glob('{}/{}_{}/rawWRF/*rcp85.stats-abs_vals.csv'.format(idir, src, sid))        
        for rfile in rstats:
            s = os.path.basename(rfile).split('_')
            gcm = s[5].upper()
            rcp = s[6].split('.')[0].replace('rcp', 'RCP ').replace('5', '.5')
            print('\t{} {}'.format(rcp, gcm))

            # Get stats values
            rdf = fpc_reader(rfile, 3)
            rdf.columns = ['Metric', 'Statistic', 'Future YRs', 'Duration', rh, bh]            
            hdf = rdf.copy()[rdf['Future YRs'] == rdf['Future YRs'].min()]
            rdf = rdf[rdf['Future YRs'] != rdf['Future YRs'].min()]
            rdf = calc_pct(rdf, hdf)
            rdf['change_or_timeseries'] = 'pct chg'
            
            # Get WYmax time series            
            rwy = rfile.replace('stats-abs_vals', 'WYmax')
            rwydf = fpc_reader(rwy, 1)
            rwydf.columns = ['Future YRs', 'Duration', rh, bh]
            rwydf.insert(1, 'Metric', 'Water Year')
            rwydf.insert(0, 'Statistic', 'Max')
            fy = rwydf['Future YRs']
            wyhis = rwydf.copy()[(fy >= 1970) & (fy <= 1999)]
            wyhis[['Future YRs']] = '1980s'
            rwydf = calc_pct(rwydf, wyhis)
            rwydf['change_or_timeseries'] = 'time series'
            
            
            # Get SNmax Time Series
            rsn = rfile.replace('stats-abs_vals', 'SNmax')
            rsndf = fpc_reader(rsn, 2)
            rsndf.columns = ['Future YRs', 'Metric', 'Duration', rh, bh]            
            rsndf.insert(0, 'Statistic', 'Max')
            fy = rsndf['Future YRs']
            snhis = rsndf.copy()[(fy >= 1970) & (fy <= 1999)]
            snhis[['Future YRs']] = '1980s'
            rsndf = calc_pct(rsndf, snhis)
            rsndf['change_or_timeseries'] = 'time series'
            
            
            rdf = rdf.append(rwydf, sort=True)
            rdf = rdf.append(rsndf, sort=True)
            rdf['Model'] = gcm
            rdf['Scenario'] = rcp
            cdf = cdf.append(rdf)
            
        c = [x for x in cdf.columns if x not in ['Model', rh, bh]]
        cg = cdf.groupby(c).mean().reset_index().round(2)
        cg['Model'] = 'Ensemble'
        cdf = cdf.append(cg, sort=True)

        
        # Merge header info and format
        for i in range(0, 5):
            cdf.insert(i, fformat[i], line[i])
        cdf = cdf[fformat]
        
        cdf.to_csv(f, index=False, header=False, na_rep='NaN', float_format='%.2f')
        
    


# Create validation_data file - pct chng vs obs
def mk_valid_file(ndf, out):
    cfilter = [oformat[0], oformat[2]]
    rh = 'rawWRF (pct bias)'
    bh = 'bcWRF (pct bias)'
    oh = 'obs'
    
    f = open(out, 'w')    
    df = pd.DataFrame(columns = vformat)
    df.to_csv(f, index=False)
    
    locs = ndf[~(ndf[cfilter].duplicated())].sort_values(by=cfilter)
    
    for idx, line in locs.iterrows():
        src = line[0]
        sid = line[2]
        print('{} {}'.format(src, sid))
        cdf = pd.DataFrame()
        
        # Get obs
        odata = glob('{}/{}_{}/OBS/*W.stats.csv'.format(idir, src, sid))
        odf = his_read(odata[0], oh) if odata else pd.DataFrame()
        
        # Get wrf data
        rdata = glob('{}/{}_{}/rawWRF/*rcp85.stats-abs_vals.csv'.format(idir, src, sid))        
        for rfile in rdata:
            s = os.path.basename(rfile).split('_')
            gcm = s[5] 
            rcp= s[6].split('.')[0].replace('rcp', 'RCP ').replace('5', '.5')
            
            rdf = his_read(rfile, rh)
            bfile = rfile.replace('rawWRF', 'bcWRF')

            # Merge bcwrf
            if os.path.exists(bfile):
                bdf = his_read(bfile, bh)
                rdf = rdf.merge(bdf)
            else:
                rdf[bh] = np.nan

            # Merge obs
            if odf.empty:
                rdf[oh] = np.nan                
            else:
                rdf = rdf.merge(odf)

            # Calculate pchnge
            rdf[rh] = ((rdf[rh] - rdf[oh]) / rdf[oh] * 100).round(2)
            rdf[bh] = ((rdf[bh] - rdf[oh]) / rdf[oh] * 100).round(2)
            
            rdf.drop(oh, inplace=True, axis=1)

            rdf['Model'] = gcm
            rdf['Scenario'] = rcp

            cdf = cdf.append(rdf)

            
        # Calculate ensemble
        gdf = cdf.groupby(['Metric', 'Statistic', 'Duration', 'Scenario']).agg(np.nanmean).reset_index()            
        gdf['Model'] = 'Ensemble'
        cdf = cdf.append(gdf)
        
        # Merge header info            
        for i in range(0, 5):
            cdf.insert(i, vformat[i], line[i])
        cdf = cdf.merge(ndf, how='left')
        cdf['Model'] = cdf['Model'].str.upper()
        
        # Clean-up and append
        cdf = cdf[vformat]
        
        cdf.to_csv(f, index=False, header=False, na_rep='NaN', float_format='%.2f')



obsf = '/home/disk/picea/mauger/2017_12_KingCounty_Stormwater/DATA/WRF_hourly/tableau/wrf_pext_tableau-Nobs.csv'
obsf = 'nobs.csv'
#odf = pd.read_csv(obsf)
odf = mk_nobs_file('{}/nobs.csv'.format(odir))
mk_idf_file(odf, '{}/idf.csv'.format(odir))
mk_fpc_file(odf, '{}/fpc.csv'.format(odir))
mk_valid_file(odf, '{}/vd.csv'.format(odir))























# get raw his -> 1975-2015
# get bc  his -> 
# get wymax raw
# get snmax raw
# get stats raw

# get wymax  bc
# get snmax  bc
# get stats bc


            
#loop through nobs (valid)
#-> get column data
#-> find obs
#-> read obs
#-> find rawWRF stats
#-> read rawWRF
#-> divide by obs
#-> find bcWRF stats
#-> read bcWRF
#-> combine raw and bc
#-> append to full df
#
#loop through nobs (idf)
#-> get column data
#-> find obs
#-> read obs
#
#loop through nobs (pct chg)
#-> get column data
#-> find rawWRF momax, wymax
#-> find bcwrf momax, wymax
#-> 

#keys = ['Metric', 'Statistic', 'Years']
#cname = keys + [x.replace('-hr Precip', '') for x in df.columns if 'hr' in x]
#df.columns = cname
#df = pd.melt(df, id_vars=keys, var_name='Duration', value_name='Precip (mm)')
#df.Duration = df.Duration.astype(int)
#df['Obs Network'] = src
#df['Station ID'] = fid
#
#df = df.merge(obs, on=['Obs Network', 'Station ID', 'Metric', 'Duration'])
#df['OBS/rawWRF/bcWRF'] = wrf
#df['Model'] = gcm
#df['Scenario'] = rcp
#df = df[hformat]
                                                
