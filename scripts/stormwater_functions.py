#!/bin/usr/python3

import pandas as pd
import numpy as np
import math
import sys


# Given a path to timeseries csv, reads and returns a formatted timeseries dataframe
def read_ts(infile):
    df = pd.read_csv(infile)
    df.columns = ['Y', 'M', 'D', 'H', '1hr']
    
    df['date'] = df.Y.map(str) + '-' + df.M.map(str) + '-' + df.D.map(str) + ' ' + df.H.map(str) + ':00'
    df.index = pd.to_datetime(df.date)

    return df
    
    
# Given a timeseries dataframe, create a running average timeseries dataset
def calc_runavg(df):
    df = df[['1hr']]
    df.columns = '1hr Precip'
    
    hrs = [2, 3, 6, 12, 24, 48, 72, 120, 240, 360]
    for h in hrs:
        df[str(h) + 'hr Precip'] = df['1hr Precip'].rolling(h).sum()
        
    return df

# Given a timeseries dataframe, calculate the monthly maximum time-series
def calc_momax(df):
    df = df.resample('M').max()

    # Formatting
    df.insert(0, 'Year', df.index.year)
    df.insert(1, 'Month', df.index.month)

    return df

# Given a timeseries dataframe, calculate the seasonal maximum time-series
def calc_snmax(df):
    df = df.shift(31*24).resample('3M').max()

    # Formatting
    df.insert(0, 'Year', df.index.year)
    df = df.reset_index(drop=True)

    seasons = ['DJF', 'MAM', 'JJA', 'SON']
    slist = [seasons[i % 3] for i in df.index]
    df.insert(1, 'Season', slist)

    return df

# Given a timeseries dataframe, calculate the water-year maximum time-series
def calc_wymax(df):
    df = df.shift(92 * 24).resample('Y').max()
    df.loc[0] = np.nan
    print(df)
    
    # Formatting
    df.insert(0, 'Water Year', df.index.year)

    return df
    
# Given a timeseries dataframe, return the dataset with only non-zero times
def filter_nonzero(df):
    df = df[df['1hr'] > 0]
    df.columns = 'Precip (mm)'

    df.insert(0, 'YYYY', df.index.year)
    df.insert(1, 'MM', df.index.month)
    df.insert(2, 'DD', df.index.day)
    df.insert(3, 'HH', df.index.hour)

    return df

# Given a timeseries dataframe and start/end time, return dataset with data within given time
def filter_years(df, styr=1678, edyr=2262):

    if styr > edyr:
        return pd.DataFrame()

    
    df = df.loc[str(styr)+'-01-01 00:00:00':str(edyr)+'-12-31 23:00:00']
    df.columns = 'Precip (mm)'

    df.insert(0, 'YYYY', df.index.year)
    df.insert(1, 'MM', df.index.month)
    df.insert(2, 'DD', df.index.day)
    df.insert(3, 'HH', df.index.hour)

    return df

# Given a yearly timeseries dataframe, calculates the extreme statistics
def calc_gev(df, yrs=[2,5,10,25,50,100]):
    df = pd.DataFrame(df)
    df.dropna(inplace=True)
    df.columns = ['value']    
    df.sort_values(by=['value'], ascending=False, inplace=True)
    df = df.reset_index(drop=True)
        
    # Calculate extreme moments
    df['rank'] = df.index + 1
    n = len(df)

    df['b1'] = (n - df['rank']) / (n * (n - 1)) * df['value']
    df['b2'] = df.b1 * (n - 1 - df['rank']) / (n - 2)
    df['b3'] = df.b2 * (n - 2 - df['rank']) / (n - 3)
    
    # Calculate extreme variables
    B = [0,0,0,0]
    L = [0,0,0,0]
    B[0] = df.value.mean()
    B[1] = df.b1.sum()
    B[2] = df.b2.sum()

    L[1] = B[0]
    L[2] = 2 * B[1] - B[0]
    L[3] = 6 * B[2] - 6 * B[1] + B[0]
    
    c = (2 * L[2]) / (L[3] + 3 * L[2]) - 0.630930
    kappa = (7.8590 * c) + (2.9554 * c * c)
    gamma = math.gamma(1 + kappa)
    alpha = kappa * L[2] / (gamma * (1 - math.pow(2, -kappa)))
    psi = L[1] + alpha * (gamma - 1) / kappa
    
    # Calculate extremes for given years
    edf = pd.DataFrame({'yrs':yrs})
    edf['prob'] = 1 - 1 / edf['yrs']
    edf['value'] = psi + alpha / kappa * (1 - (np.power(-np.log(edf['prob']), kappa)))
    #edf.ix[:, edf.columns != 'prob'] = (edf.ix[:, edf.columns != 'prob'].astype(np.double))

    return edf        
 
    
# Given water-year, seasonal, and monthly dataframes, returns extreme statistics dataframe
def calc_stats(wy_df, sn_df, mo_df, isHis=True):

    # Index dataframes
    wy_df.index = pd.to_datetime(wy_df['Water Year'].map(str) + '-01-01')
    sn_df.index = pd.to_datetime(sn_df['Year'].map(str) + '-01-01')
    mo_df.index = pd.to_datetime(mo_df['Year'].map(str) + '-' + mo_df['Month'].map(str) + '-01')
    
    
    df = pd.DataFrame()
    hrs = [2, 3, 6, 12, 24, 48, 72, 120, 240, 360]
    hrs = [str(h)+' hr Precip' for h in hrs]

    if isHis:
        styr = wy_df[['Water Year']][0]
        edyr = wy_df[['Water Year']][-1]
        pds = [(styr, edyr)]
    else:
        pds = [(1970, 1999), (1980, 2009), (2020, 2049), (2030, 2059)
               (2040, 2069), (2050, 2099), (2060, 2099), (2070, 2099)]
        pds = [(1970, 1999)]

    for styr,edyr in pds:
        # Process water years
        pd_df = wy_df.loc[str(styr)+'01-01':str(edyr)+'12-31']

        
        ddf = pd.DataFrame()
        
        # Calculate extreme statistics
        for col in pd_df.columns[1:]:
            hr_df = calc_gev(pd_df[[col]])            
            ddf = ddf.append(hr_df)

        # Flip dataframe
        
        
        # Combine total sum

        # Label metrics and period

        # Combine results to main dataframe

        df = df.append(ddf)
        
    
    print(df)
    return df
    
        # Process seasonal

        # Process monthly


# Given obs and fut statistical dataframe, returns percent change of statistical dataframe
def calc_pchg(xf, yf):
    print(1)








    
a = 'COOP_457781_47.43621N_121.49286W_rawWRF_access-1-0_rcp45.csv'
out = '../runavg.csv'


ra_df = pd.read_csv(out, index_col=0, parse_dates=[0], infer_datetime_format=True)
print(ra_df)

#ra_df = calc_runavg(a)
#ra_df.to_csv(out)


#mo_df = calc_momax(ra_df)
#print("MO MAX")
#print(mo_df)

wy_df = calc_wymax(ra_df)
print("WY MAX")
print(wy_df)

sn_df = calc_snmax(ra_df)
print('SN MAX:')
print(sn_df)

stats = calc_stats(wy_df, sn_df, mo_df)

print(stats)
