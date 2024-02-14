import os
import pandas as pd
import numpy as np

src = 'C:/Users/czhe0008/Documents/DLCprojects/openpose/data_conversion/p2_data/csv/'
dest = 'C:/Users/czhe0008/Documents/DLCprojects/openpose/data_conversion/p2_data/csv/'
# joints = ["W_x","W_y","W_z","TCMC_x","TCMC_y","TCMC_z","TMCP_x","TMCP_y","TMCP_z","TIP_x","TIP_y","TIP_z","TT_x","TT_y","TT_z","IMCP_x","IMCP_y","IMCP_z","IPIP_x","IPIP_y","IPIP_z","IDIP_x","IDIP_y","IDIP_z","IT_x","IT_y","IT_z","MMCP_x","MMCP_y","MMCP_z","MPIP_x","MPIP_y","MPIP_z","MDIP_x","MDIP_y","MDIP_z","MT_x","MT_y","MT_z","RMCP_x","RMCP_y","RMCP_z","RPIP_x","RPIP_y","RPIP_z","RDIP_x","RDIP_y","RDIP_z","RT_x","RT_y","RT_z","LMCP_x","LMCP_y","LMCP_z","LPIP_x","LPIP_y","LPIP_z","LDIP_x","LDIP_y","LDIP_z","LT_x","LT_y","LT_z"]
threshold = 5   # When to start flagging jumps
osc_threshold = 1   # When to start flagging oscillations
increment = 10  # How fast to expand acceptable range

for foldername in os.listdir(src):
    trial = os.path.join(src, foldername)
    for file in os.listdir(trial):
        if file[-9:] == 'df_3d.csv':
            data = os.path.join(src, foldername, file)
            df = pd.read_csv(data)
            print(file)
            #calculating flexion angles
            count = 0
            for joint in list(df):
                col = df.loc[:,joint]
                diff = threshold
                check = 0
                # Checking for large jumps
                for i in range(1,len(col)):
                    if abs(col[i] - col[check]) > diff:
                        diff = diff + increment
                    else:
                        grad = (col[i]-col[check])/(i-check)
                        for j in range(check+1,i):
                            col[j] = col[j-1] + grad
                        diff = threshold
                        check = i
                # Checking for oscillations
                for i in range(2,len(col)):
                    if (abs(col[i-1] - col[i-2]) > osc_threshold) & (abs(col[i] - col[i-1]) > osc_threshold) & (np.sign(col[i-1] - col[i-2])!=np.sign(col[i] - col[i-1])):
                        col[i-1] = col[i-2] + (col[i] - col[i-2])/2
                df.loc[:,joint] = col
            df.to_csv(dest+foldername+"_filtered.csv", index = False)
