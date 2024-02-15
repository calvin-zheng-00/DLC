import os
import pandas as pd
import numpy as np

src = 'C:/Users/czhe0008/Documents/DLCprojects/openpose/data_conversion/p2_data/filtered/'
dest = 'C:/Users/czhe0008/Documents/DLCprojects/openpose/data_conversion/p2_data/angles/'
finger = ["TCMC_x","TMCP_x","TIP_x","IMCP_x","IPIP_x","IDIP_x","MMCP_x","MPIP_x","MDIP_x","RMCP_x","RPIP_x","RDIP_x","LMCP_x","LPIP_x","LDIP_x"]
arm = ["RE_x"]
abduct = ["IMCP_x","MMCP_x","RMCP_x"]

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
            for joint in finger:
                col_name = joint[:-2] + "_flex"
                idx = df.columns.get_loc(joint)
                if count % 3 == 0:
                    v1 = [df.iloc[:, idx+3]-df.iloc[:, idx],df.iloc[:, idx+4]-df.iloc[:, idx+1],df.iloc[:, idx+5]-df.iloc[:, idx+2]]
                    v2 = [df.loc[:, 'RW_x']-df.iloc[:, idx],df.loc[:, 'RW_y']-df.iloc[:, idx+1],df.loc[:, 'RW_z']-df.iloc[:, idx+2]]
                else:
                    v1 = [df.iloc[:, idx+3]-df.iloc[:, idx],df.iloc[:, idx+4]-df.iloc[:, idx+1],df.iloc[:, idx+5]-df.iloc[:, idx+2]]
                    v2 = [df.iloc[:, idx-3]-df.iloc[:, idx],df.iloc[:, idx-2]-df.iloc[:, idx+1],df.iloc[:, idx-1]-df.iloc[:, idx+2]]
                v1mag = np.sqrt(np.square(v1[0]) + np.square(v1[1]) + np.square(v1[2]))
                v1norm = [v1[0].div(v1mag), v1[1].div(v1mag), v1[2].div(v1mag)]
                v2mag = np.sqrt(np.square(v2[0]) + np.square(v2[1]) + np.square(v2[2]))
                v2norm = [v2[0].div(v2mag), v2[1].div(v2mag), v2[2].div(v2mag)]
                res = v1norm[0] * v2norm[0] + v1norm[1] * v2norm[1] + v1norm[2] * v2norm[2]
                df2[col_name] = np.arccos(res)
                count += 1
            for joint in arm:
                col_name = joint[:-2] + "_flex"
                idx = df.columns.get_loc(joint)
                v1 = [df.iloc[:, idx+3]-df.iloc[:, idx],df.iloc[:, idx+4]-df.iloc[:, idx+1],df.iloc[:, idx+5]-df.iloc[:, idx+2]]
                v2 = [df.iloc[:, idx-3]-df.iloc[:, idx],df.iloc[:, idx-2]-df.iloc[:, idx+1],df.iloc[:, idx-1]-df.iloc[:, idx+2]]
                v1mag = np.sqrt(np.square(v1[0]) + np.square(v1[1]) + np.square(v1[2]))
                v1norm = [v1[0].div(v1mag), v1[1].div(v1mag), v1[2].div(v1mag)]
                v2mag = np.sqrt(np.square(v2[0]) + np.square(v2[1]) + np.square(v2[2]))
                v2norm = [v2[0].div(v2mag), v2[1].div(v2mag), v2[2].div(v2mag)]
                res = v1norm[0] * v2norm[0] + v1norm[1] * v2norm[1] + v1norm[2] * v2norm[2]
                df2[col_name] = np.arccos(res)

            ## Wrist angle
            #v1 = [df.loc[:, 'MMCP_x']-df.loc[:, 'RW_x'],df.loc[:, 'MMCP_y']-df.loc[:, 'RW_y'],df.loc[:, 'MMCP_z']-df.loc[:, 'RW_z']]
            #v2 = [df.loc[:, 'RE_x']-df.loc[:, 'TCMC_x'],df.loc[:, 'RE_y']-df.loc[:, 'TCMC_y'],df.loc[:, 'RE_z']-df.loc[:, 'TCMC_z']]
            #v1mag = np.sqrt(np.square(v1[0]) + np.square(v1[1]) + np.square(v1[2]))
            #v1norm = [v1[0].div(v1mag), v1[1].div(v1mag), v1[2].div(v1mag)]
            #v2mag = np.sqrt(np.square(v2[0]) + np.square(v2[1]) + np.square(v2[2]))
            #v2norm = [v2[0].div(v2mag), v2[1].div(v2mag), v2[2].div(v2mag)]
            #res = v1norm[0] * v2norm[0] + v1norm[1] * v2norm[1] + v1norm[2] * v2norm[2]
            #df2["W_flex"] = np.arccos(res)

            # Calculating the normal vector to the plane that is the palm (using the wrist, index knuckle and ring knuckle)
            v1 = [df.loc[:, 'IMCP_x']-df.loc[:, 'RW_x'],df.loc[:, 'IMCP_y']-df.loc[:, 'RW_y'],df.loc[:, 'IMCP_z']-df.loc[:, 'RW_z']]
            v2 = [df.loc[:, 'RMCP_x']-df.loc[:, 'RW_x'],df.loc[:, 'RMCP_y']-df.loc[:, 'RW_y'],df.loc[:, 'RMCP_z']-df.loc[:, 'RW_z']]
            v1 = np.transpose(v1)
            v2 = np.transpose(v2)
            norm = np.cross(v1,v2) # An array of groups of 3 coordinates. [[norm1],[norm2],[norm3],...]
            normt = np.transpose(norm)  # Normal to the palm plane
            normmag = np.sqrt(np.square(normt[0]) + np.square(normt[1]) + np.square(normt[2]))
            unitnorm = [np.divide(normt[0],normmag), np.divide(normt[1],normmag), np.divide(normt[2],normmag)]

            # Calculating angle of abduction
                # Alternative adjustment in the future: Calculate angle between MCP-PIP and MCP-Wrist instead of between two different fingers
            for joint in abduct:
                col_name = joint[:-2] + "_abd"
                idx = df.columns.get_loc(joint)
                v1 = np.array([df.iloc[:, idx+3]-df.iloc[:, idx],df.iloc[:, idx+4]-df.iloc[:, idx+1],df.iloc[:, idx+5]-df.iloc[:, idx+2]])
                v2 = np.array([df.iloc[:, idx+15]-df.iloc[:, idx+12],df.iloc[:, idx+16]-df.iloc[:, idx+13],df.iloc[:, idx+17]-df.iloc[:, idx+14]])
                v1_np = norm*((v1[0]*normt[0] + v1[1] * normt[1] + v1[2] * normt[2])/np.square(normmag))[:, None]   # Vector component normal to the plane
                v1_np = np.transpose(v1_np)
                v1_p = [v1[0]-v1_np[0],v1[1]-v1_np[1],v1[2]-v1_np[2]]                                       # Vector component on the plane
                v1mag = np.sqrt(np.square(v1_p[0]) + np.square(v1_p[1]) + np.square(v1_p[2]))
                v1norm = [np.divide(v1_p[0],v1mag), np.divide(v1_p[1],v1mag), np.divide(v1_p[2],v1mag)]
                v2_np = norm*((v2[0]*normt[0] + v2[1] * normt[1] + v2[2] * normt[2])/np.square(normmag))[:, None]
                v2_np = np.transpose(v2_np)
                v2_p = [v2[0]-v2_np[0],v2[1]-v2_np[1],v2[2]-v2_np[2]]
                v2mag = np.sqrt(np.square(v2_p[0]) + np.square(v2_p[1]) + np.square(v2_p[2]))
                v2norm = [np.divide(v2_p[0],v2mag), np.divide(v2_p[1],v2mag), np.divide(v2_p[2],v2mag)]
                res = v1norm[0] * v2norm[0] + v1norm[1] * v2norm[1] + v1norm[2] * v2norm[2]
                df2[col_name] = np.arccos(res)
            
            # Thumb abduction (from TMCP to TCMC to IMCP)
            v1 = [df.loc[:, 'TMCP_x']-df.loc[:, 'TCMC_x'],df.loc[:, 'TMCP_y']-df.loc[:, 'TCMC_y'],df.loc[:, 'TMCP_z']-df.loc[:, 'TCMC_z']]
            v2 = [df.loc[:, 'IMCP_x']-df.loc[:, 'TCMC_x'],df.loc[:, 'IMCP_y']-df.loc[:, 'TCMC_y'],df.loc[:, 'IMCP_z']-df.loc[:, 'TCMC_z']]
            v1mag = np.sqrt(np.square(v1[0]) + np.square(v1[1]) + np.square(v1[2]))
            v1norm = [v1[0].div(v1mag), v1[1].div(v1mag), v1[2].div(v1mag)]
            v2mag = np.sqrt(np.square(v2[0]) + np.square(v2[1]) + np.square(v2[2]))
            v2norm = [v2[0].div(v2mag), v2[1].div(v2mag), v2[2].div(v2mag)]
            res = v1norm[0] * v2norm[0] + v1norm[1] * v2norm[1] + v1norm[2] * v2norm[2]
            df2["TMCP_abd"] = np.arccos(res)

            # Find wrist flex vs abduction
            # We have the normal of the wrist plane and two points for the orthogonal plane to pass through. (MMCP and Wrist)
            # Find the cross product between the normal to the plane and the vector (MMCP-Wrist) to get the plane
            # along which we will be finding the flex angle. Subtract flex angle from total angle to get abduction angle.
            vp = [df.loc[:, 'MMCP_x']-df.loc[:, 'RW_x'],df.loc[:, 'MMCP_y']-df.loc[:, 'RW_y'],df.loc[:, 'MMCP_z']-df.loc[:, 'RW_z']]
            flex_plane = np.cross(norm,np.transpose(vp))
            flex_plane_t = np.transpose(flex_plane)
            flex_plane_mag = np.sqrt(np.square(flex_plane_t[0]) + np.square(flex_plane_t[1]) + np.square(flex_plane_t[2]))
            v1 = np.array([df.loc[:, 'RE_x']-df.loc[:, 'RW_x'],df.loc[:, 'RE_y']-df.loc[:, 'RW_y'],df.loc[:, 'RE_z']-df.loc[:, 'RW_z']])
            v1_np = flex_plane*((v1[0]*flex_plane_t[0] + v1[1] * flex_plane_t[1] + v1[2] * flex_plane_t[2])/np.square(flex_plane_mag))[:, None]   # Vector component normal to the plane
            v1_np = np.transpose(v1_np)
            v1_p = [v1[0]-v1_np[0],v1[1]-v1_np[1],v1[2]-v1_np[2]]                                       # Vector component on the plane
            v1pmag = np.sqrt(np.square(v1_p[0]) + np.square(v1_p[1]) + np.square(v1_p[2]))
            v1pnorm = [np.divide(v1_p[0],v1pmag), np.divide(v1_p[1],v1pmag), np.divide(v1_p[2],v1pmag)]
            vpmag = np.sqrt(np.square(vp[0]) + np.square(vp[1]) + np.square(vp[2]))
            vpnorm = [np.divide(vp[0],vpmag), np.divide(vp[1],vpmag), np.divide(vp[2],vpmag)]
            res = v1pnorm[0] * vpnorm[0] + v1pnorm[1] * vpnorm[1] + v1pnorm[2] * vpnorm[2]
            df2["W_flex"] = np.arccos(res)
            # Abuct = angle between original wrist angle and flex angle
            v1mag = np.sqrt(np.square(v1[0]) + np.square(v1[1]) + np.square(v1[2]))
            v1norm = [np.divide(v1[0],v1mag), np.divide(v1[1],v1mag), np.divide(v1[2],v1mag)]
            res = v1norm[0] * v1pnorm[0] + v1norm[1] * v1pnorm[1] + v1norm[2] * v1pnorm[2]
            df2["W_abd"] = np.arccos(res)
            
            # To calculate rotation of wrist, find the plane generated by the wrist, elbow and shoulder (arm plane).
            # Find the angle between the palm plane and the arm plane.
            v1 = [df.loc[:, 'RW_x']-df.loc[:, 'RE_x'],df.loc[:, 'RW_y']-df.loc[:, 'RE_y'],df.loc[:, 'RW_z']-df.loc[:, 'RE_z']]
            v2 = [df.loc[:, 'RS_x']-df.loc[:, 'RE_x'],df.loc[:, 'RS_y']-df.loc[:, 'RE_y'],df.loc[:, 'RS_z']-df.loc[:, 'RE_z']]
            v1 = np.transpose(v1)
            v2 = np.transpose(v2)
            arm_norm = np.cross(v1,v2) # An array of groups of 3 coordinates. [[norm1],[norm2],[norm3],...]
            arm_normt = np.transpose(arm_norm)
            arm_normmag = np.sqrt(np.square(arm_normt[0]) + np.square(arm_normt[1]) + np.square(arm_normt[2]))
            unit_arm_norm = [np.divide(arm_normt[0],arm_normmag), np.divide(arm_normt[1],arm_normmag), np.divide(arm_normt[2],arm_normmag)]
            # Find plane perpendicular to RE_RW vector; RE_RW is the normal to the plane
            EWnormt = np.array([df.loc[:, 'RW_x']-df.loc[:, 'RE_x'],df.loc[:, 'RW_y']-df.loc[:, 'RE_y'],df.loc[:, 'RW_z']-df.loc[:, 'RE_z']])
            EWnorm = np.transpose(EWnormt)
            EWnormmag = np.sqrt(np.square(EWnormt[0]) + np.square(EWnormt[1]) + np.square(EWnormt[2]))
            # Find angle between palm normal and arm normal
            palm_np = EWnorm*((unitnorm[0]*EWnormt[0] + unitnorm[1] * EWnormt[1] + unitnorm[2] * EWnormt[2])/np.square(EWnormmag))[:, None]
            palm_np = np.transpose(palm_np)
            palm_p = [unitnorm[0]-palm_np[0],unitnorm[1]-palm_np[1],unitnorm[2]-palm_np[2]]
            palmmag = np.sqrt(np.square(palm_p[0]) + np.square(palm_p[1]) + np.square(palm_p[2]))
            palmnorm = [np.divide(palm_p[0],palmmag), np.divide(palm_p[1],palmmag), np.divide(palm_p[2],palmmag)]
            res = palmnorm[0] * unit_arm_norm[0] + palmnorm[1] * unit_arm_norm[1] + palmnorm[2] * unit_arm_norm[2]
            df2["W_rot"] = np.arccos(res)

            # Find shoulder abduction by using Shoulder, Chest, Pelvis to find the plane
            v1 = [df.loc[:, 'P_x']-df.loc[:, 'C_x'],df.loc[:, 'P_y']-df.loc[:, 'C_y'],df.loc[:, 'P_z']-df.loc[:, 'C_z']]
            v2 = [df.loc[:, 'RS_x']-df.loc[:, 'C_x'],df.loc[:, 'RS_y']-df.loc[:, 'C_y'],df.loc[:, 'RS_z']-df.loc[:, 'C_z']]
            v1 = np.transpose(v1)
            v2 = np.transpose(v2)
            chest_norm = np.cross(v1,v2) # An array of groups of 3 coordinates. [[norm1],[norm2],[norm3],...]
            chest_normt = np.transpose(chest_norm)
            chest_normmag = np.sqrt(np.square(chest_normt[0]) + np.square(chest_normt[1]) + np.square(chest_normt[2]))
            
            v1 = np.array([df.loc[:, 'RE_x']-df.loc[:, 'RS_x'],df.loc[:, 'RE_y']-df.loc[:, 'RS_y'],df.loc[:, 'RE_z']-df.loc[:, 'RS_z']])
            v1_np = chest_norm*((v1[0]*chest_normt[0] + v1[1] * chest_normt[1] + v1[2] * chest_normt[2])/np.square(chest_normmag))[:, None]
            v1_np = np.transpose(v1_np)
            v1_p = [v1[0]-v1_np[0],v1[1]-v1_np[1],v1[2]-v1_np[2]]
            v1mag = np.sqrt(np.square(v1_p[0]) + np.square(v1_p[1]) + np.square(v1_p[2]))
            v1norm = [np.divide(v1_p[0],v1mag), np.divide(v1_p[1],v1mag), np.divide(v1_p[2],v1mag)]
            v2 = np.array([df.loc[:, 'C_x']-df.loc[:, 'RS_x'],df.loc[:, 'C_y']-df.loc[:, 'RS_y'],df.loc[:, 'C_z']-df.loc[:, 'RS_z']])
            v2_np = chest_norm*((v2[0]*chest_normt[0] + v2[1] * chest_normt[1] + v2[2] * chest_normt[2])/np.square(chest_normmag))[:, None]
            v2_np = np.transpose(v2_np)
            v2_p = [v2[0]-v2_np[0],v2[1]-v2_np[1],v2[2]-v2_np[2]]
            v2mag = np.sqrt(np.square(v2_p[0]) + np.square(v2_p[1]) + np.square(v2_p[2]))
            v2norm = [np.divide(v2_p[0],v2mag), np.divide(v2_p[1],v2mag), np.divide(v2_p[2],v2mag)]
            res = v1norm[0] * v2norm[0] + v1norm[1] * v2norm[1] + v1norm[2] * v2norm[2]
            df2["RS_flex"] = np.arccos(res)

            v0mag = np.sqrt(np.square(v1[0]) + np.square(v1[1]) + np.square(v1[2]))
            v0norm = [np.divide(v1[0],v0mag), np.divide(v1[1],v0mag), np.divide(v1[2],v0mag)]
            res = v1norm[0] * v0norm[0] + v1norm[1] * v0norm[1] + v1norm[2] * v0norm[2]
            df2["RS_rot"] = np.arccos(res)
            
            # Change src with dest
            df2.to_csv(dest+foldername+"_angles.csv", index = False)

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