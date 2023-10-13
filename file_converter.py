#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys
import pickle
import itertools

#fpath = "e3v830e-20230918T105454-105508.h5"
#fout = "data2.csv"
#df = pd.read_hdf(fpath)
#df.to_csv(fout, index=False)

#objects = []
#with (open("detections.pickle", "rb")) as openfile:
#    while True:
#        try:
#            objects.append(pickle.load(openfile))
#        except EOFError:
#            break

df_list = pd.read_pickle(r'detections.pickle')
df = pd.DataFrame()
for i in range(len(df_list)):
    #df['i'] = df_list[i]
    additional = pd.DataFrame({
    str(i): df_list[i]
})
    df = pd.concat([df, additional], axis=1) 
df.to_csv("data3.csv", index=False)

#list_t = list(map(list, itertools.zip_longest(df_list, fillvalue=None)))
#, columns=['1','2','3','4','5','6','7']
#df = pd.DataFrame(df_list)
#df.to_csv("data3.csv", index=False)
