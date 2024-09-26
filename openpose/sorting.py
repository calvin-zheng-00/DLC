import os
import shutil

src_folder = '/home/czhe0008/sz11_scratch/calvin/p21_json'
dest_folder = '/home/czhe0008/sz11_scratch/calvin/p21_json_s'
for foldername in os.listdir(src_folder):
    trial = foldername[8:]
    new_path = os.path.join(dest_folder, trial)
    new_path = dest_folder + '/' + trial
    if not os.path.exists(new_path):
        os.makedirs(new_path)
        print(trial)
    old_path = os.path.join(src_folder, foldername)
    shutil.move(old_path, new_path)