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
def calc_runsum(df, t=1):
    df = df[['1hr']]
    df.columns = ['1-hr Precip']
    
    hrs = [2, 3, 6, 12, 24, 48, 72, 120, 240, 360]
    for h in hrs:
        n = df.shape[1]
        roll = df['1-hr Precip'].rolling(h, min_periods=h*t).sum()
        df.insert(n, '{}-hr Precip'.format(h), roll)
    return df


# Given a timeseries dataframe, calculate the water-year aggregated time-series
def calc_wyagg(df, aggtype=max, t=0.95):
    ct = df.shift(92 * 24).resample('Y').count()
    df = df.shift(92 * 24).resample('Y').agg(aggtype)
    ff = np.repeat(ct.index.is_leap_year[:, np.newaxis], df.shape[1], axis=1)
    ff = (ff + 365) * 24 * t    
    df[ct < ff] = np.nan
    
    # Formatting
    df.insert(0, 'Water Year', df.index.year)

    return df

# Given a timeseries dataframe, calculate the seasonal aggregated time-series
def calc_snagg(df, aggtype=max, t=0.95):
    df = df.copy()
    seasons = ['DJF', 'MAM', 'JJA', 'SON']    
    
    df['Season'] = np.floor(df.index.month / 3).astype(int) % 4    
    df.index = df.index.shift(31*24, freq='H')                      # Shift to keep december first
    df.insert(0, 'Year', df.index.year)    
    ct = df.groupby(['Year', 'Season']).count().reset_index()
    df = df.groupby(['Year', 'Season']).agg(aggtype).reset_index()
    df.index = pd.to_datetime(df.apply(lambda x: '{}-{:02d}-01'.format(int(x.Year), int(x.Season*3+2)), axis=1))
    ff = ((df.index.is_leap_year)*1 + df.index.days_in_month + 61 + (df.index.month==2)*1)*24 * t
    ff = np.repeat(np.array(ff)[:, np.newaxis], df.shape[1], axis=1)
    mask = np.array(ct < ff)
    mask[:, 0:2] = False
    df[mask] = np.nan
    df['Season'] = df['Season'].apply(lambda x: seasons[x])
    
    return df
    

# Given a timeseries dataframe, calculate the monthly aggregated time-series
def calc_moagg(df, aggtype=max, t=0.95):
    ct = df.resample('M').count()    
    df = df.resample('M').agg(aggtype)
    ff = np.repeat(np.array(ct.index.days_in_month)[:, np.newaxis], df.shape[1], axis=1)
    ff = ff * 24 * t
    df[ct < ff] = np.nan
    
    # Formatting
    df.insert(0, 'Year', df.index.year)
    df.insert(1, 'Month', df.index.month)

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
def process_stats_df(sdf, pdf, stats, period):
    sdf = sdf[[ x for x in pdf.columns if 'hr' in x]]
    pdf = pdf[[ x for x in pdf.columns if 'hr' in x]]
    df = pd.DataFrame()

    # Combine total sum
    ddf = pd.DataFrame(sdf.sum()).T
    ddf = ddf / sdf.count()    
    ddf.iloc[0] = ddf.iloc[0,0].repeat(len(ddf.columns))
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
def calc_stats_wy(ra_df, wy_df, periods):
    wy_df.index = pd.to_datetime(wy_df['Water Year'].map(str) + '-01-01')
    sdf = calc_wyagg(ra_df, sum)
    
    df = pd.DataFrame()
    for styr, edyr in periods:
        psdf = sdf.loc[str(styr)+'-01-01':str(edyr)+'-12-31']
        pdf = wy_df.loc[str(styr)+'-01-01':str(edyr)+'-12-31']
        ddf = process_stats_df(psdf, pdf, 'Water Year', '{}-{}'.format(styr, edyr))
        df = df.append(ddf)
                
    # Sort and relabel results
    df['index'] = df.index
    df.sort_values(['index', 'Years'], inplace=True)
    df.drop(['index'], axis=1, inplace=True)

    return df


# Given seasonal dataframe and list of periods, return extreme statistics dataframe
def calc_stats_sn(ra_df, sn_df, periods):
    sn_df.index = pd.to_datetime(sn_df['Year'].map(str) + '-01-01')
    sdf = calc_snagg(ra_df, sum)    
    sns = ['DJF', 'MAM', 'JJA', 'SON']
    
    df = pd.DataFrame()
    for styr, edyr in periods:
        for sn in sns:
            psdf = sdf[sdf['Season'] == sn].loc[str(styr)+'-01-01':str(edyr)+'-12-31']
            pdf = sn_df[sn_df['Season'] == sn].loc[str(styr)+'-01-01':str(edyr)+'-12-31']
            ddf = process_stats_df(psdf, pdf, sn, '{}-{}'.format(styr, edyr))
            df = df.append(ddf)

    # Sort results
    df['index'] = df.index
    sorter = dict(zip(sns, range(len(sns))))
    df['sn_index'] = df['Water Year or Month'].map(sorter)    
    df.sort_values(['sn_index', 'index', 'Years'], inplace=True)
    df.drop(['index', 'sn_index'], axis=1, inplace=True)
    
    return df



# Given monthly dataframe and list of periods, return extreme statistics dataframe
def calc_stats_mo(ra_df, mo_df, periods):
    mo_df.index = pd.to_datetime(mo_df['Year'].map(str) + '-' + mo_df['Month'].map(str) + '-01')
    sdf = calc_moagg(ra_df, sum)
    df = pd.DataFrame()
    
    # Loop through each month for each period
    for (styr,edyr), mo in [(x, y) for x in periods for y in range(1,13)]:
        psdf = sdf[sdf['Month'] == mo].loc[str(styr)+'-01-01':str(edyr)+'-12-31']
        pdf = mo_df[mo_df['Month'] == mo].loc[str(styr)+'-01-01':str(edyr)+'-12-31']
        ddf = process_stats_df(psdf, pdf, mo, '{}-{}'.format(styr, edyr))
        df = df.append(ddf)

    # Sort results
    df['index'] = df.index
    df.sort_values(['Water Year or Month', 'index', 'Years'], inplace=True)
    df.drop(['index'], axis=1, inplace=True)
    
    # Replace Month values
    df['Water Year or Month'] = df['Water Year or Month'].apply(lambda x: calendar.month_name[x])
    
    return df
        

# Given water-year, seasonal, and monthly dataframes, returns combined extreme statistics dataframe
def calc_stats(ra_df, wy_df, sn_df, mo_df, isHis=True):
    df = pd.DataFrame()
    
    if isHis:
        styr = int(ra_df.ix[0, 'Water Year'])
        edyr = int(ra_df.ix[-1, 'Water Year'])
        periods = [(styr, edyr)]
    else:
        periods = [(1970, 1999), (1980, 2009), (2020, 2049), (2030, 2059),
               (2040, 2069), (2050, 2099), (2060, 2099), (2070, 2099)]
        
    # Process water years
    df = df.append(calc_stats_wy(ra_df, wy_df, periods))    
    # Process seasonal
    df = df.append(calc_stats_sn(ra_df, sn_df, periods))
    # Process monthly
    df = df.append(calc_stats_mo(ra_df, mo_df, periods))

    return df

# Given a statistical dataframe, returns percent change of statistical dataframe
def calc_pchg(st):
    his_period = ['1970-1999', '1980-2009']
    his = st[[x in his_period for x in st['Years']]]
    fut = st[[x not in his_period for x in st['Years']]]
    mon = list(his.columns[0:2])
    hrs = list(his.columns[3:])

    foj = his.merge(fut, on=mon)
    foj.rename(columns={'Years_y': 'Future Years', 'Years_x':'Historical Years'}, inplace=True)
    
    for h in hrs:
        n = foj.shape[1]
        fyr = foj[h + '_y']
        hyr = foj[h + '_x']
        pct = (fyr - hyr) / hyr * 100
        foj.insert(n, h, pct)

    foj = foj[mon + ['Future Years', 'Historical Years'] + hrs]
    return foj.round(1)



# Output dataframe as csv with comment above header
def csv_comment(df, out, comment, index=False):
    f = open(out, 'a')
    f.write(comment)
    df.to_csv(f, index=index)
    f.close()
