import pandas as pd
import os
import glob
import numpy as np

targetFolder = "P2_3"

path = os.getcwd()
path = path + "\\" + targetFolder + "\\trial1\\pose-3d"
csv_files = glob.glob(os.path.join(path, "*.csv"))

for f in csv_files:
	
	df = pd.read_csv(f)
	row_names = {'0': [],
        '1': [],
        '2': [],
		'3': [],
		'4': [],
		'5': [],
		'6': [],
		'7': [],
		'percent dropout': [],}
	out = pd.DataFrame({},row_names)
	fname = "out" + f.split("\\")[-1]
	
	df2 = df.filter(regex='_ncams$',axis=1)
	for col in df2.columns:
		ncams = [df2[col][df2[col]==0].count(), df2[col][df2[col]==1].count(), df2[col][df2[col]==2].count(), df2[col][df2[col]==3].count(),
		   df2[col][df2[col]==4].count(), df2[col][df2[col]==5].count(), df2[col][df2[col]==6].count(), df2[col][df2[col]==7].count(),
		   df2[col][(df2[col]==0) | (df2[col]==1)].count() / df2[col].count()]
		out[col] = ncams
	out.to_csv(fname)
