import os
import pandas as pd
import numpy as np

src = 'C:/Users/czhe0008/Documents/DLCprojects/openpose/data_conversion/confirmation/new/'
dest = 'C:/Users/czhe0008/Documents/DLCprojects/openpose/data_conversion/confirmation/new/'
# joints = ["W_x","W_y","W_z","TCMC_x","TCMC_y","TCMC_z","TMCP_x","TMCP_y","TMCP_z","TIP_x","TIP_y","TIP_z","TT_x","TT_y","TT_z","IMCP_x","IMCP_y","IMCP_z","IPIP_x","IPIP_y","IPIP_z","IDIP_x","IDIP_y","IDIP_z","IT_x","IT_y","IT_z","MMCP_x","MMCP_y","MMCP_z","MPIP_x","MPIP_y","MPIP_z","MDIP_x","MDIP_y","MDIP_z","MT_x","MT_y","MT_z","RMCP_x","RMCP_y","RMCP_z","RPIP_x","RPIP_y","RPIP_z","RDIP_x","RDIP_y","RDIP_z","RT_x","RT_y","RT_z","LMCP_x","LMCP_y","LMCP_z","LPIP_x","LPIP_y","LPIP_z","LDIP_x","LDIP_y","LDIP_z","LT_x","LT_y","LT_z"]
threshold = 3   # When to start flagging jumps
osc_threshold = 1   # When to start flagging oscillations
increment = 10  # How fast to expand acceptable range
window = 4

for foldername in os.listdir(src):
    trial = os.path.join(src, foldername)
    for file in os.listdir(trial):
        if file[-9:] == 'df_3d.csv':
            data = os.path.join(src, foldername, file)
            df = pd.read_csv(data)
            df2 = pd.DataFrame()
            print(file)
            #calculating flexion angles
            count = 0
            for joint in list(df):
                col = df.loc[:,joint]
                averaged = []
                for i in range(1,len(col), window):
                    averaged.append(np.mean(col[i:(i+window-1)]))
                df2[joint] = averaged
            df2.to_csv(dest+foldername+"_filtered.csv", index = False)