import pandas as pd
import numpy as np
import json
import os
from aniposelib.boards import CharucoBoard, Checkerboard
from aniposelib.cameras import Camera, CameraGroup
from aniposelib.utils import load_pose2d_fnames
from freemocap.core_processes.capture_volume_calibration.triangulate_3d_data import triangulate_3d_data


cgroup = CameraGroup.load('calibration.toml')
score_threshold = 0.15
tracked_points = 21
frames = 0
main_directory = 'C:/Users/czhe0008/Documents/DLCprojects/openpose/data_conversion'
src = 'p2_fix'
dest = 'p2_pose3d'
#model = 'DLC_resnet50_V2Aug14shuffle1_650000'
hand_markers_names = ['W','TCMC','TMCP','TIP','TT','IMCP','IPIP','IDIP','IT','MMCP','MPIP','MDIP','MT','RMCP','RPIP','RDIP','RT','LMCP','LPIP','LDIP','LT']
#body_marker_names = ['N','C','RS','RE','RW','LS','LE','LR','P','RP','RK','RF','LP','LK','LF','REye','RC','LEye','LC','LH','LToe','LA','RH','RToe','RA']
#header = pd.MultiIndex.from_product([[model], hand_markers_names, ['x','y','likelihood']])

src_dir = os.path.join(main_directory, src)
# iterate over trials
for trial in os.listdir(src_dir):
    trial_dir = os.path.join(src_dir, trial)
    # Checking camera dropouts
    count = 0
    for foldername in os.listdir(trial_dir):
        count += 1
    if count != 7:
        print('\n\n\n\nError, wrong number of cameras in: ' + trial + '\n\n\n\n')
        continue
    out_dir = os.path.join(main_directory, dest, trial)
    os.makedirs(out_dir)
    print('\n\n\n\n' + trial + '\n\n\n\n')
    output_arr = np.empty((0,1))
    output_arr_flag = True
    # iterate cams in the trial
    for foldername in os.listdir(trial_dir):
        camera = os.path.join(trial_dir, foldername)
        if os.path.isfile(camera):
            continue
        hand_data = np.empty((0,1))
        hand_data_flag = True
        for filename in os.listdir(camera):
            file = os.path.join(camera, filename)
            # checking if it is a file
            if os.path.isfile(file):
                # Opening JSON file
                f = open(file)
                data = json.load(f)
                if len(data['people']) == 0:
                    if hand_data_flag:
                        hand_data = np.empty((0,tracked_points*2)) #Hardcoding tracked points here
                        hand_data_flag = False
                    hand_data = np.append(hand_data, [np.zeros(tracked_points*2)], axis=0)
                else:
                    #body_markers = data['people'][0]['pose_keypoints_2d']
                    hand_right_markers = data['people'][0]['hand_right_keypoints_2d']
                    thresholder = np.array(hand_right_markers).astype(np.float32)
                    thresholder = thresholder[2::3]
                    thresholder[thresholder < score_threshold] = np.nan # 0
                    thresholder[thresholder >= score_threshold] = 1
                    scores = [None]*len(thresholder)*2
                    scores[::2] = thresholder
                    scores[1::2] = thresholder
                    del hand_right_markers[2::3]
                    hand_right_thresholded = np.multiply(scores,hand_right_markers)
                    if hand_data_flag:
                        hand_data = np.empty((0,tracked_points*2))
                        hand_data_flag = False
                    hand_data = np.append(hand_data, [hand_right_thresholded], axis=0) # hand_right_thresholded
                # Closing file
                f.close()
        hand_data_x = np.copy(hand_data[:,::2])
        hand_data_y = np.copy(hand_data[:,1::2])
        
    # num_camsXnum_framesXnum_tracked_pointsXnum_spatial_dimensions
    # Order from back to front
        num_frames = np.empty((0,1))
        num_frames_flag = True
        for i in range(0,np.shape(hand_data_x)[0]):
            num_points = np.empty((0,2))
            for j in range(0,np.shape(hand_data_x)[1]):
                num_points = np.append(num_points, [[hand_data_x[i,j],hand_data_y[i,j]]], axis=0)
            if num_frames_flag:
                num_frames = np.empty((0,np.shape(hand_data_x)[1],2))
                num_frames_flag = False
            num_frames = np.append(num_frames,[num_points], axis=0)
        if output_arr_flag:
            output_arr = np.empty((0,np.shape(hand_data_x)[0],np.shape(hand_data_x)[1],2))
            frames = np.shape(hand_data_x)[0]
            output_arr_flag = False
        if np.shape(num_frames)[0] > frames:
            num_frames = num_frames[:(frames-np.shape(num_frames)[0]),:,:]
        else:
            for i in range(0,frames-np.shape(num_frames)[0]):
                num_frames = np.append(num_frames,[np.zeros((tracked_points, 2))], axis=0)
        
    #Check numpy dimensions
    # If dimensions don't match, add 0s until they do, or remove rows until they do.

        output_arr = np.append(output_arr,[num_frames], axis=0)
        #hand_df = pd.DataFrame(hand_data, columns = header)
        #hand_df.to_csv('test.csv', index=False)
        #out_dir = "pose-2d/" + foldername + ".h5"
        #hand_df.to_hdf(out_dir, key='data', mode='w')
        flag = 1
    
    #ncams, nframes, npoints, ndims = triangulate_3d_data(
    #    anipose_calibration_object,
    #    mediapipe_2d_data: np.ndarray,
    #    output_data_folder_path: Union[str, Path],
    #    mediapipe_confidence_cutoff_threshold: float,
    #    use_triangulate_ransac: bool = False,
    #    kill_event: multiprocessing.Event = None,
    #)
    #    number_of_cameras = mediapipe_2d_data.shape[0]
    #    number_of_frames = mediapipe_2d_data.shape[1]
    #    number_of_tracked_points = mediapipe_2d_data.shape[2]
    #    number_of_spatial_dimensions = mediapipe_2d_data.shape[3]


    output_points_3d, output_reprojection_error = triangulate_3d_data(
            anipose_calibration_object=cgroup,
            mediapipe_2d_data=output_arr,
            output_data_folder_path=out_dir,
            mediapipe_confidence_cutoff_threshold=0.6,
            use_triangulate_ransac=False
        )
    # mediapipe_confidence_cutoff_threshold currently does nothing
    pose_3d = np.empty((0,3*np.shape(output_points_3d)[1]))
    for i in range(0,np.shape(output_points_3d)[0]):
        pose_3d_temp = np.empty((0))
        for j in range(0,np.shape(output_points_3d)[1]):
            xyz = output_points_3d[i,j,:]
            xyz.flatten()
            pose_3d_temp = np.append(pose_3d_temp,xyz, axis=0)
        pose_3d = np.append(pose_3d,[pose_3d_temp], axis=0)

    #out_header = ['W','TCMC','TMCP','TIP','TT','IMCP','IPIP','IDIP','IT','MMCP','MPIP','MDIP','MT','RMCP','RPIP','RDIP','RT','LMCP','LPIP','LDIP','LT']
    out_suffix = ['_x','_y','_z']*tracked_points
    out_header = np.repeat(hand_markers_names, 3)
    out_header = np.char.add(out_header,out_suffix)

    df_3d = pd.DataFrame(pose_3d, columns=out_header)
    df_err = pd.DataFrame(output_reprojection_error, columns=hand_markers_names)
    df_3d.to_csv(out_dir+"/df_3d.csv", index = False)
    df_err.to_csv(out_dir+"/df_err.csv", index = False)


#fpath = 'C:/Users/czhe0008/Documents/DLCprojects/test/data_conversion/e3v830e-20231004T152604-152626.h5'

#df = pd.read_hdf(fpath)

#fname_dict = {
#    '0e': 'e3v830e-20231004T152604-152626.h5',
#    '1f': 'e3v831f-20231004T152604-152626.h5',
#    '3c': 'e3v833c-20231004T152604-152626.h5',
#    '4c': 'e3v834c-20231004T152604-152626.h5',
#    '37': 'e3v8337-20231004T152604-152626.h5',
#    '80': 'e3v8380-20231004T152604-152626.h5',
#    '89': 'e3v8389-20231004T152604-152626.h5',
#}

#d = load_pose2d_fnames(fname_dict, cam_names=cgroup.get_names())

#score_threshold = 0.5

#n_cams, n_points, n_joints, _ = d['points'].shape
#points = d['points']
#scores = d['scores']

#bodyparts = d['bodyparts']