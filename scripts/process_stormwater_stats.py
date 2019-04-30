#!/bin/usr/python3
import stormwater_functions as sf
import pandas as pd

a = 'COOP_457781_47.43621N_121.49286W_rawWRF_access-1-0_rcp45.csv'
out = '../runavg.csv'


ra_df = pd.read_csv(out, index_col=0, parse_dates=[0], infer_datetime_format=True)
#ra_df = sf.calc_runsum(a)
#ra_df.to_csv(out)


mo_df = sf.calc_momax(ra_df)
print("MO MAX")
print(mo_df)

wy_df = sf.calc_wymax(ra_df)
print("WY MAX")
print(wy_df)

sn_df = sf.calc_snmax(ra_df)
print('SN MAX')
print(sn_df)

stats = sf.calc_stats(wy_df, sn_df, mo_df, False)

print(stats)

