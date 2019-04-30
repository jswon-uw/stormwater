#!/usr/bin/env python
import pandas as pd
import numpy as np
import calendar
import math
import sys


# Given a path to timeseries csv, reads and returns a formatted timeseries dataframe
def read_ts(infile):
    df = pd.read_csv(infile)
    df.columns = ['Y', 'M', 'D', 'H', '1hr']
    
    df['date'] = df.Y.map(str) + '-' + df.M.map(str) + '-' + df.D.map(str) + ' ' + df.H.map(str) + ':00'
    df.index = pd.to_datetime(df.date)

    return df
    
    
# Given a timeseries dataframe, create a running sum timeseries dataset
def calc_runsum(df):
    print(df)
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
    slist = [seasons[i % 4] for i in df.index]
    df.insert(1, 'Season', slist)

    return df

# Given a timeseries dataframe, calculate the water-year maximum time-series
def calc_wymax(df):
    df = df.shift(92 * 24).resample('Y').max()
    df.iloc[0] = np.nan
    
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
    
    return edf[['yrs', 'value']]


# Helper function to process each period dataframe merging common functions in each calc_stats function
def process_stats_df(pdf, stats, period):
    pdf = pdf[[ x for x in pdf.columns if 'hr' in x]]
    df = pd.DataFrame()

    # Combine total sum
    ddf = pd.DataFrame(pdf.sum()).T
    ddf.insert(0, 'Water Year or Month', stats)
    ddf.insert(1, 'Statistic', 'Total')
    ddf.insert(2, 'Years', period)
    df = df.append(ddf)

    
    # Calculate extreme statistics
    ddf = pd.DataFrame()
    for col in pdf.columns:
        hr_df = calc_gev(pdf[[col]])
        hr_df.index = hr_df['yrs']
        hr_df = hr_df[['value']]
        hr_df.columns = [col]        
        ddf = pd.concat([ddf, hr_df], axis=1)
        
    ddf.insert(0, 'Water Year or Month', stats)    
    ddf.insert(1, 'Statistic', ddf.index.map(str) + '-yr')        
    ddf.insert(2, 'Years', period) 
    df = df.append(ddf)

    return df
    

# Given water-year dataframe and list of periods, return extreme statistics dataframe
def calc_stats_wy(wy_df, periods):
    wy_df.index = pd.to_datetime(wy_df['Water Year'].map(str) + '-01-01')

    df = pd.DataFrame()
    for styr, edyr in periods:
        pdf = wy_df.loc[str(styr)+'-01-01':str(edyr)+'-12-31']
        ddf = process_stats_df(pdf, 'Water Year', '{}-{}'.format(styr, edyr))
        df = df.append(ddf)
                
    # Sort and relabel results
    df['index'] = df.index
    df.sort_values(['index', 'Years'], inplace=True)
    df.drop(['index'], axis=1, inplace=True)
    
    return df


# Given seasonal dataframe and list of periods, return extreme statistics dataframe
def calc_stats_sn(sn_df, periods):
    sn_df.index = pd.to_datetime(sn_df['Year'].map(str) + '-01-01')
    sns = ['DJF', 'MAM', 'JJA', 'SON']    
    
    df = pd.DataFrame()
    for styr, edyr in periods:
        for sn in sns:
            pdf = sn_df[sn_df['Season'] == sn].loc[str(styr)+'-01-01':str(edyr)+'-12-31']
            ddf = process_stats_df(pdf, sn, '{}-{}'.format(styr, edyr))
            df = df.append(ddf)

    # Sort results
    df['index'] = df.index
    sorter = dict(zip(sns, range(len(sns))))
    df['sn_index'] = df['Water Year or Month'].map(sorter)    
    df.sort_values(['sn_index', 'index', 'Years'], inplace=True)
    df.drop(['index', 'sn_index'], axis=1, inplace=True)
    
    return df



# Given monthly dataframe and list of periods, return extreme statistics dataframe
def calc_stats_mo(mo_df, periods):
    mo_df.index = pd.to_datetime(mo_df['Year'].map(str) + '-' + mo_df['Month'].map(str) + '-01')
    df = pd.DataFrame()
    
    # Loop through each month for each period
    for (styr,edyr), mo in [(x, y) for x in periods for y in range(1,13)]:
        pdf = mo_df[mo_df['Month'] == mo].loc[str(styr)+'-01-01':str(edyr)+'-12-31']
        ddf = process_stats_df(pdf, mo, '{}-{}'.format(styr, edyr))
        df = df.append(ddf)

    # Sort results
    df['index'] = df.index
    df.sort_values(['Water Year or Month', 'index', 'Years'], inplace=True)
    df.drop(['index'], axis=1, inplace=True)
    
    # Replace Month values
    df['Water Year or Month'] = df['Water Year or Month'].apply(lambda x: calendar.month_name[x])
    
    return df
        

# Given water-year, seasonal, and monthly dataframes, returns combined extreme statistics dataframe
def calc_stats(wy_df, sn_df, mo_df, isHis=True):
    df = pd.DataFrame()
    
    if isHis:
        styr = int(wy_df.ix[0, 'Water Year'])
        edyr = int(wy_df.ix[-1, 'Water Year'])
        periods = [(styr, edyr)]
    else:
        periods = [(1970, 1999), (1980, 2009), (2020, 2049), (2030, 2059),
               (2040, 2069), (2050, 2099), (2060, 2099), (2070, 2099)]
        
    # Process water years
    df = df.append(calc_stats_wy(wy_df, periods))
        
    # Process seasonal
    df = df.append(calc_stats_sn(sn_df, periods))
        
    # Process monthly
    df = df.append(calc_stats_mo(mo_df, periods))

    return df

# Given obs and fut statistical dataframe, returns percent change of statistical dataframe
def calc_pchg(xf, yf):
    print(1)
