import os
import shutil

src_folder = r'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\openpose_output'
dest_folder = r'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\p2_sorted'
for foldername in os.listdir(src_folder):
    trial = foldername[8:]
    new_path = dest_folder + '/' + trial
    if not os.path.exists(new_path):
        os.makedirs(new_path)
        print(trial)
    shutil.move(src_folder +'/' + foldername, new_path)