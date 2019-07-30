#!/bin/usr/python3
import os.path
import sys
import argparse
import pandas as pd
import numpy as np
from datetime import datetime as dt
import calendar as cal

###############################################################################
# Creates a training table based on ratio of historical and observation data
# binned by 1% interval percentiles
###############################################################################
force_zero=0.025  # Force all value under value to 0

def print_full(x):
    pd.set_option('display.max_rows', len(x))
    print(x)
    pd.reset_option('display.max_rows')    



# Read and format dataframe for training
def format_df(fname, styr=1900, edyr=2099):
    hdr = ['Y', 'M', 'D', 'H', 'val']
    naList = ["NA", "NaN", ""]
    
    if(not os.path.isfile(fname)):
        print("File does not exist: " + fname)
        sys.exit()    
    
    df = pd.read_csv(fname, na_values=naList, names=hdr, low_memory=False, skiprows=[0])
    df.index = pd.to_datetime(df.Y.map(str) + '-' + df.M.map(str) + '-' + df.D.map(str) + ' ' + df.H.map(str) + ':00')

    # Filter dataframe to given start and end year
    df = df[(df.index.year >= styr) & (df.index.year <= edyr)]
    return df


# Grabs all data points that fall on given dofy within given window
def get_data(df):
    dfT = df[['val']].dropna().sort_values(by='val').reset_index(drop=True)
    return dfT

# Calculate percentile ranking
def pct_bin(series, n):

    # Calculate percentile
    ranked = np.ceil(series.rank(method='max', pct=1) * 100 / n * 100).astype(int)    
    df = pd.DataFrame({'val': series, 'rank':ranked})

    # Fit 0 values into their own bracket
    df['rank'] = df['rank'] - (df['val'] == 0) * 1

    # Aggregate bins
    df = df.groupby(['rank']).agg(['min', 'mean']).reset_index()
    df.columns = ['rank', 'floor', 'mean']

    # Return full frame
    df = df.set_index('rank').reindex(np.arange(1, n+1)).fillna(method='ffill')
    df = df.fillna(method='bfill').reset_index()
    return df
        
    
# Create training from observation and historical
def create_training(obs_file, his_file, wdw=False, obs_st=1950, obs_ed=2017, his_st=1900, his_ed=2015):
    # ----------------------------------------------------------------------------
    # Read OBS, HIS
    # ----------------------------------------------------------------------------
    print("---: " + obs_file.split('/')[-1])    

    obs = format_df(obs_file, obs_st, obs_ed)
    his = format_df(his_file, his_st, his_ed) 

    # Data length check
    if (len(obs) == 0) | (len(his) == 0):
        print('OBS: {}'.format(len(obs)))
        print('HIS: {}'.format(len(his)))
        print("Not enough data")
        sys.exit()
        

    # ----------------------------------------------------------------------------
    # Training Algorigthm
    # ----------------------------------------------------------------------------

    # Clean data
    obsT = get_data(obs)
    hisT = get_data(his)

    # Calculate percentile ranking
    obsR = pct_bin(obsT['val'], 100)
    hisR = pct_bin(hisT['val'], 100)

    

    # Calculate ratio
    hisR['ratio'] = obsR['mean'] / hisR['mean']
    hisR = hisR.set_index('rank').drop_duplicates()
    hisR = hisR[['floor', 'ratio']]
    hisR['ratio'].fillna('0', inplace=True)

    return hisR

def read_training(train_file):
    tdf = pd.read_csv(train_file)
    return tdf 
    

# apply bias-correction
def apply_bc(sim_file, tdf):
    sim = format_df(sim_file)
    tdf[tdf['floor'] < force_zero] = 0
    
    #print(tdf)

    n = round(sim.val.max()) * 100

    # Sort values into bins and evaluate multiplying factor
    bidx = np.digitize(sim['val'], tdf['floor'])
    factor = tdf.reset_index().loc[bidx -1, 'ratio'].values

    # Apply bias-correction
    sim['factor'] = factor
    sim['val'] = (sim['val'].astype(float) * factor.astype(float)).round(2)

    # Clean-up
    sim = sim[['Y', 'M', 'D', 'H', 'val']]
    sim.columns = ['YYYY', 'MM', 'DD', 'HH', 'Precip (mm)']        
    return sim

