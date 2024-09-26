import os
import shutil
import pandas as pd

src_folder = 'C:/Users/czhe0008/Documents/DLCprojects/openpose/data_conversion/massive_3d/3d/p19_3d_pre'
dest_folder = 'C:/Users/czhe0008/Documents/DLCprojects/openpose/data_conversion/massive_3d/3d/p19_3d_sorted'

df = pd.read_csv("C:/Users/czhe0008/Documents/DLCprojects/openpose/data_conversion/activity_order/csv/instructions_p19.csv")
instruct_list = df.videos
instruct_no = 0
count = -1

for foldername in sorted(os.listdir(src_folder)):
    if count == 4:
        count = 0
        instruct_no += 1
    else:
        count += 1

    ### Check for faulty trials here, and continue to next loop
    #if (instruct_no + 1 == 4) & (count >= 2):
    #    count = 0
    #    instruct_no += 1

    new_path = os.path.join(dest_folder, instruct_list[instruct_no])
    if not os.path.exists(new_path):
        os.makedirs(new_path)
    curr_path = os.path.join(src_folder, foldername)
    #shutil.copytree(curr_path, new_path)
    shutil.move(curr_path, new_path)