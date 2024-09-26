import os
import pandas as pd
import numpy as np
from scipy import signal

src = 'C:/Users/czhe0008/Documents/DLCprojects/openpose/data_conversion/massive_3d/p16_3d'
dest = 'C:/Users/czhe0008/Documents/DLCprojects/openpose/data_conversion/massive_3d/p16_filter'
#joints = ["W_x","W_y","W_z","TCMC_x","TCMC_y","TCMC_z","TMCP_x","TMCP_y","TMCP_z","TIP_x","TIP_y","TIP_z","TT_x","TT_y","TT_z",
#          "IMCP_x","IMCP_y","IMCP_z","IPIP_x","IPIP_y","IPIP_z","IDIP_x","IDIP_y","IDIP_z","IT_x","IT_y","IT_z",
#          "MMCP_x","MMCP_y","MMCP_z","MPIP_x","MPIP_y","MPIP_z","MDIP_x","MDIP_y","MDIP_z","MT_x","MT_y","MT_z",
#          "RMCP_x","RMCP_y","RMCP_z","RPIP_x","RPIP_y","RPIP_z","RDIP_x","RDIP_y","RDIP_z","RT_x","RT_y","RT_z",
#          "LMCP_x","LMCP_y","LMCP_z","LPIP_x","LPIP_y","LPIP_z","LDIP_x","LDIP_y","LDIP_z","LT_x","LT_y","LT_z",
#          "N_x","N_y","N_z","C_x","C_y","C_z","RS_x","RS_y","RS_z","RE_x","RE_y","RE_z","RW_x","RW_y","RW_z",
#          "LS_x","LS_y","LS_z","LE_x","LE_y","LE_z","LR_x","LR_y","LR_z","P_x","P_y","P_z"]
window = 5 # For boxcar filtering. Currently unused
threshold = 30   # When to start flagging jumps
increment = 20  # How fast to expand acceptable range

for foldername in os.listdir(src):
    trial = os.path.join(src, foldername)
    for file in os.listdir(trial):
        if file[-9:] == 'df_3d.csv':
            data = os.path.join(src, foldername, file)
            df = pd.read_csv(data)
            #df = df.interpolate(method='linear', limit_direction='forward', axis=0)
            df2 = pd.DataFrame()
            print(file)
            count = 0
            for joint in list(df):
                col = df.loc[:,joint]
                col = np.asarray(col)
                diff = threshold
                check = 0
                #print(joint)
                # Filling blanks
                if np.isnan(col[-1]):
                    for i in range(1,len(col)):
                        if not np.isnan(col[i]):
                            col[-1] = col[i]
                            break
                for i in range(len(col)-2,1,-1):
                    if np.isnan(col[i]):
                        col[i] = col[i+1]
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
                # Butterworth
                b, a = signal.butter(2, 5, fs=40)
                if len(col) < 3 * max(len(a), len(b)):
                    filtered = signal.filtfilt(b, a, col, padlen=0)
                else:
                    filtered = signal.filtfilt(b, a, col)
                # Boxcar
                #win = signal.windows.boxcar(window)
                #filtered = signal.convolve(col, win, mode='same') / sum(win)
                df2[joint] = filtered
            df2 = df2.interpolate(method='linear', limit_direction='forward', axis=0)
            df2.to_csv(os.path.join(dest,foldername+"_filtered.csv"), index = False)