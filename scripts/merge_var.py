#!/bin/usr/python3

import pandas as pd
from glob import glob
import os.path as osp
from datetime import datetime
import sys
import argparse
from pprint import pprint

sys.path.append('/home/disk/tsuga2/jswon11/workdir/stormwater/scripts/')
#sys.path.append('/home/disk/tsuga2/jswon11/workdir/2019_04_stormwater-newwrf/scripts/pyeto/')
import humidity_lib as hlib
import hspf_cloud_fraction as hcf

parser = argparse.ArgumentParser()
parser.add_argument('-l', default='', help='location specifier')
parser.add_argument('-g', default='', help='gcm specifier')
parser.add_argument('-i', nargs='?', default=None, help='input directory')
args = parser.parse_args()
lreg = args.l
greg = args.g
idir = args.i if args.i else '/home/disk/tsuga2/jswon11/workdir/2019_04_stormwater-newwrf/data/pub/'


def read_ts(base, var):
    fname = base.replace('bcWRF', 'rawWRF').replace('pub', var)
    if not osp.exists(fname):
        print('{} does not exist: {}'.format(var, osp.basename(fname)))
        sys.exit()
    print('\tReading {}'.format(var))    
    df = pd.read_csv(fname)
    c = df.columns
    return df[c[-1]]


    
def main():        
    data = glob('{}/*{}*/bcWRF/*{}*rcp85.csv'.format(idir, lreg, greg))
    pprint(data)
    i = 0
    n = len(data)
    
    for d in data:
        i = i + 1
        print('{}/{}: {}'.format(i, n, d))
        prec = pd.read_csv(d)

        # Clear existing 
        if len(prec) > 5:
            prec = prec[['YYYY', 'MM', 'DD', 'HH', 'Precip (mm)']]
        
        
        t2 = read_ts(d, 'T2') - 273.15        
        q2 = read_ts(d, 'Q2')
        glw = read_ts(d, 'GLW')
        swd = read_ts(d, 'SWDOWN')
        u10 = read_ts(d, 'U10')
        v10 = read_ts(d, 'V10')
        wnd = (u10**2 + v10**2)**0.5
        psr = read_ts(d, 'PSFC')

        dew = hlib.dewpoint_approximation_vector(t2, q2, psr)
        pet = hlib.pet_fao56_vector(t2, q2, psr, wnd, swd, glw)
        cld = hcf.cloud_fraction_vector(glw, t2)
        
        print('\tSaving...')
        prec.insert(4, 'T 2m (C)', t2)
        prec.insert(6, 'Tdew 2m (C)', dew)
        prec.insert(7, '10m Wind Speed (m/s)', wnd)
        prec.insert(8, 'Incoming Shortwave (W/m2)', swd)
        prec.insert(9, 'Incoming Longwave (W/m2)', glw)
        prec.insert(10, 'Pot. Evap. (mm/day)', pet)
        prec.insert(11, 'Pressure (Pa)', psr)
        prec.insert(12, 'Specific Humidity (kg/kg)', q2)
        prec.insert(13, 'Cloud Fraction (-)', cld)

        prec.to_csv(d, index=False, float_format='%.5f')
        
main()




#for d in dd:
#    x = glob('{}/bcWRF/*rcp85.csv'.format(d))    
#    #pprint(x)
#    print(d)
#    print(len(x))


#YYYY,MM,DD,HH,
#T 2m (C),
#Precip (mm),
#Tdew 2m (C), -> humid: T, Q, PSFC
#10m Wind Speed (m/s),
#Incoming Shortwave (W/m2),
#Incoming Longwave (W/m2),
#Pot. Evap. (mm/day), -> humid: T, Q, PSFC, W, SWDOWN, GLW
#Pressure (Pa),
#Specific Humidity (kg/kg),
#Cloud Fraction (-) -> cloud: GLW, T2
