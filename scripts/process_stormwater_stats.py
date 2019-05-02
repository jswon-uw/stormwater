#!/usr/bin/env python
import stormwater_functions as sf
import pandas as pd
import sys

a = 'COOP_457781_47.43621N_121.49286W_rawWRF_access-1-0_rcp45.csv'
out = '../runavg.csv'


ra_df = pd.read_csv(out, index_col=0, parse_dates=[0], infer_datetime_format=True)
#ra_df = sf.calc_runsum(a)
#ra_df.to_csv(out)


df = pd.DataFrame()

mo_df = sf.calc_moagg(ra_df)
print("MO MAX")
#print(mo_df)

wy_df = sf.calc_wyagg(ra_df)
print("WY MAX")

sn_df = sf.calc_snagg(ra_df)
print('SN MAX')
#sn_df.to_csv('sn.csv', index=False, float_format='%.2f')

stats = sf.calc_stats(ra_df, wy_df, sn_df, mo_df, False)
#stats.to_csv('stats_check.csv', float_format='%.2f')

print(stats)



