import os
import shutil

src_folder = '/home/czhe0008/sz11_scratch/calvin/p21_json'
dest_folder = '/home/czhe0008/sz11_scratch/calvin/p21_json_double'
prev_start = 0
prev_end = 0
prev_folder = 'temp'
for foldername in sorted(os.listdir(src_folder)):
    starttime = int(foldername[9:15])
    endtime = int(foldername[16:])
    if (abs(starttime - prev_start) <= 1) or (abs(endtime - prev_end) <= 1):
        #trial = foldername[8:]
        new_path = os.path.join(dest_folder, prev_folder)
        if not os.path.exists(new_path):
            os.makedirs(new_path)
            print(foldername)
        curr_path = os.path.join(src_folder, foldername)
        for camera in os.listdir(curr_path):
            cam_path = os.path.join(curr_path, camera)
            shutil.move(cam_path, new_path)
        prev_path = os.path.join(src_folder, prev_folder)
        for camera in os.listdir(prev_path):
            cam_path = os.path.join(prev_path, camera)
            shutil.move(cam_path, new_path)
    else:
        prev_start = starttime
        prev_end = endtime
        prev_folder = foldername
for foldername in sorted(os.listdir(src_folder)):
    new_path = os.path.join(src_folder, foldername)
    test = len(os.listdir(new_path))
    if len(os.listdir(new_path)) == 0: 
        os.rmdir(new_path)