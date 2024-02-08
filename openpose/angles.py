import os
import pandas as pd
import numpy as np

src = 'C:/Users/czhe0008/Documents/DLCprojects/openpose/data_conversion/p2_pose3d'
dest = 'C:/Users/czhe0008/Documents/DLCprojects/openpose/data_conversion/p2_angles/'
finger = ["TCMC_x","TMCP_x","TIP_x","IMCP_x","IPIP_x","IDIP_x","MMCP_x","MPIP_x","MDIP_x","RMCP_x","RPIP_x","RDIP_x","LMCP_x","LPIP_x","LDIP_x"]
#abd = ["TMCP_x","IMCP_x","MMCP_x","RMCP_x","LMCP_x"]

for foldername in os.listdir(src):
    trial = os.path.join(src, foldername)
    for file in os.listdir(trial):
        if file[-9:] == 'df_3d.csv':
            data = os.path.join(src, foldername, file)
            df = pd.read_csv(data)
            print(file)
            #calculating flexion angles
            count = 0
            for joint in finger:
                col_name = joint[:-2] + "_flex"
                idx = df.columns.get_loc(joint)
                if count % 3 == 0:
                    v1 = [df.iloc[:, idx+3]-df.iloc[:, idx],df.iloc[:, idx+4]-df.iloc[:, idx+1],df.iloc[:, idx+5]-df.iloc[:, idx+2]]
                    v2 = [df.iloc[:, 0]-df.iloc[:, idx],df.iloc[:, 1]-df.iloc[:, idx+1],df.iloc[:, 2]-df.iloc[:, idx+2]]
                else:
                    v1 = [df.iloc[:, idx+3]-df.iloc[:, idx],df.iloc[:, idx+4]-df.iloc[:, idx+1],df.iloc[:, idx+5]-df.iloc[:, idx+2]]
                    v2 = [df.iloc[:, idx-3]-df.iloc[:, idx],df.iloc[:, idx-2]-df.iloc[:, idx+1],df.iloc[:, idx-1]-df.iloc[:, idx+2]]
                v1mag = np.sqrt(np.square(v1[0]) + np.square(v1[1]) + np.square(v1[2]))
                v1norm = [v1[0].div(v1mag), v1[1].div(v1mag), v1[2].div(v1mag)]
                v2mag = np.sqrt(np.square(v2[0]) + np.square(v2[1]) + np.square(v2[2]))
                v2norm = [v2[0].div(v2mag), v2[1].div(v2mag), v2[2].div(v2mag)]
                res = v1norm[0] * v2norm[0] + v1norm[1] * v2norm[1] + v1norm[2] * v2norm[2]
                df[col_name] = np.arccos(res)
                count += 1
            # Calculating angle of abduction

            df.to_csv(dest+foldername+".csv", index = False)

##Given points A,B and C, find the angle ABC
# v1 = {A.x - B.x, A.y - B.y, A.z - B.z}
# v2 = {C.x - B.x, C.y - B.y, C.z - B.z}

##The dot product of v1 and v2 is a function of the cosine of the angle between them
##(it's scaled by the product of their magnitudes). So first normalize v1 and v2:
# v1mag = sqrt(v1.x * v1.x + v1.y * v1.y + v1.z * v1.z)
# v1norm = {v1.x / v1mag, v1.y / v1mag, v1.z / v1mag}

# v2mag = sqrt(v2.x * v2.x + v2.y * v2.y + v2.z * v2.z)
# v2norm = {v2.x / v2mag, v2.y / v2mag, v2.z / v2mag}

##Dot product:
# res = v1norm.x * v2norm.x + v1norm.y * v2norm.y + v1norm.z * v2norm.z

# angle = acos(res)