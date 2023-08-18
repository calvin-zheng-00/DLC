#!/usr/bin/env python3
import pandas as pd
import sys
fpath = sys.argv[1]
if len(sys.argv)>2:
    key = sys.argv[2]
    df = pd.read_hdf(fpath, key=key)
else:
    df = pd.read_hdf(fpath)

df.to_csv(sys.stdout, index=False)


# If your data is stored in hdf5 file called data.h5 and you have saved this script as hdf2df.py then:
#
# $ python hdf2df.py data.hf > data.csv