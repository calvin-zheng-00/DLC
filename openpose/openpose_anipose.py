import pandas as pd
import numpy as np
import json
import os
from aniposelib.boards import CharucoBoard, Checkerboard
from aniposelib.cameras import Camera, CameraGroup
from aniposelib.utils import load_pose2d_fnames
from triangulate_3d_data import triangulate_3d_data
#from freemocap.core_processes.capture_volume_calibration.triangulate_3d_data import triangulate_3d_data


cgroup = CameraGroup.load('calibration.toml')
score_threshold = 0.3 # 0.15
tracked_points_hand = 21 # joints in hand
tracked_points_body = 25 # joints in body
tracked_points = tracked_points_hand + tracked_points_body
frames = 0
main_directory = 'C:/Users/czhe0008/Documents/DLCprojects/openpose/data_conversion/massive_3d/json_calib_test'
src = 'p18'
dest = '3d'
#model = 'DLC_resnet50_V2Aug14shuffle1_650000'
hand_markers_names = ['W','TCMC','TMCP','TIP','TT','IMCP','IPIP','IDIP','IT','MMCP','MPIP','MDIP','MT','RMCP','RPIP','RDIP','RT','LMCP','LPIP','LDIP','LT']
body_markers_names = ['N','C','RS','RE','RW','LS','LE','LR','P','RP','RK','RF','LP','LK','LF','REye','RC','LEye','LC','LH','LToe','LA','RH','RToe','RA']
marker_names = hand_markers_names + body_markers_names
#header = pd.MultiIndex.from_product([[model], hand_markers_names, ['x','y','likelihood']])

src_dir = os.path.join(main_directory, src)
# iterate over trials
for trial in sorted(os.listdir(src_dir)):
    trial_dir = os.path.join(src_dir, trial)
    out_dir = os.path.join(main_directory, dest, trial)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    print('\n\n\n\n' + trial + '\n\n\n\n')
    output_arr = np.empty((0,1))
    output_arr_flag = True
    cam_drop = np.empty((0,1))
    cam_drop_flag = 1
    # Counting tracking dropouts
    drop_count = []
    drop_count_flag1 = True
    # Tracking different camera names for repojection error output
    cam_names = []
    # iterate cams in the trial
    for foldername in sorted(os.listdir(trial_dir)):
        camera = os.path.join(trial_dir, foldername)
        if os.path.isfile(camera):
            continue
        # Checking for cameras to ignore. May remove in the future
        if foldername[0:7] == 'e3v82e4':
            continue
        # End of checking
        cam_names.append(foldername[0:7])
        hand_data = np.empty((0,1))
        body_data = np.empty((0,1))
        data_flag = True
        drop_count_loop = []
        drop_count_flag2 = True
        for filename in sorted(os.listdir(camera)):
            file = os.path.join(camera, filename)
            # checking if it is a file
            if os.path.isfile(file):
                #if os.path.getsize(file) == 0:
                #    if data_flag:
                #        hand_data = np.empty((0,tracked_points_hand*2)) #Hardcoding tracked points here
                #        body_data = np.empty((0,tracked_points_body*2))
                #        data_flag = False
                #    temp_arr = np.empty(tracked_points_hand*2)
                #    temp_arr[:] = np.nan
                #    hand_data = np.append(hand_data, [temp_arr], axis=0)
                #    body_data = np.append(body_data, [temp_arr], axis=0)
                #else:
                # Opening JSON file
                f = open(file)
                data = json.load(f)
                if len(data['people']) == 0:
                    if data_flag:
                        hand_data = np.empty((0,tracked_points_hand*2)) #Hardcoding tracked points here
                        body_data = np.empty((0,tracked_points_body*2))
                        data_flag = False
                    temp_arr = np.empty(tracked_points_hand*2)
                    temp_arr[:] = np.nan
                    temp_arr2 = np.empty(tracked_points_body*2)
                    temp_arr2[:] = np.nan
                    hand_data = np.append(hand_data, [temp_arr], axis=0)
                    body_data = np.append(body_data, [temp_arr2], axis=0)
                else:
                    # Getting right hand coordinates
                    hand_right_markers = data['people'][0]['hand_right_keypoints_2d']
                    thresholder = np.array(hand_right_markers).astype(np.float32)
                    thresholder = thresholder[2::3]
                    thresholder[thresholder < score_threshold] = np.nan # 0
                    thresholder[thresholder >= score_threshold] = 1
                    # Subthread for counting dropouts
                    #temp = np.nan_to_num(thresholder)
                    #if drop_count_flag2:
                    #    drop_count_loop = temp
                    #    drop_count_flag2 = False
                    #else:
                    #    drop_count_loop = np.append([drop_count_loop],[temp])
                    # end of subthread
                    scores = [None]*len(thresholder)*2
                    scores[::2] = thresholder
                    scores[1::2] = thresholder
                    del hand_right_markers[2::3]
                    hand_right_thresholded = np.multiply(scores,hand_right_markers)
                    if data_flag:
                        hand_data = np.empty((0,tracked_points_hand*2))
                        # data_flag = False # Don't set to false when it still needs to check body
                    hand_data = np.append(hand_data, [hand_right_thresholded], axis=0) # hand_right_thresholded

                    # Getting body coordinates
                    body_markers = data['people'][0]['pose_keypoints_2d']
                    thresholder = np.array(body_markers).astype(np.float32)
                    thresholder = thresholder[2::3]
                    thresholder[thresholder < score_threshold] = np.nan # 0
                    thresholder[thresholder >= score_threshold] = 1
                    scores = [None]*len(thresholder)*2
                    scores[::2] = thresholder
                    scores[1::2] = thresholder
                    del body_markers[2::3]
                    body_thresholded = np.multiply(scores,body_markers)
                    if data_flag:
                        body_data = np.empty((0,tracked_points_body*2))
                        data_flag = False
                    body_data = np.append(body_data, [body_thresholded], axis=0) # hand_right_thresholded
                # Closing file
                f.close()
        hand_data = np.concatenate((hand_data,body_data), axis = 1)
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
        #Check numpy dimensions. If dimensions don't match, add 0s until they do, or remove rows until they do.
        if np.shape(num_frames)[0] > frames:
            num_frames = num_frames[:(frames-np.shape(num_frames)[0]),:,:]
        else:
            for i in range(0,frames-np.shape(num_frames)[0]):
                num_frames = np.append(num_frames,[np.zeros((tracked_points, 2))], axis=0)

        output_arr = np.append(output_arr,[num_frames], axis=0)
    
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


    try:
        output_points_3d, output_reprojection_error, output_repoerr_3d = triangulate_3d_data(
                anipose_calibration_object=cgroup,
                mediapipe_2d_data=output_arr,
                use_triangulate_ransac=False
            )
        #Removed arguments:
                #output_data_folder_path=out_dir,
                #mediapipe_confidence_cutoff_threshold=0.6,
    except:
        print('\n\n\n\nMissing cameras in trial: ' + trial + '\n\n\n\n')
    else:
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
        out_header = np.repeat(marker_names, 3)
        out_header = np.char.add(out_header,out_suffix)

        df_3d = pd.DataFrame(pose_3d, columns=out_header)
        df_err = pd.DataFrame(output_reprojection_error, columns=marker_names)
        df_3d.to_csv(out_dir+"/df_3d.csv", index = False)
        df_err.to_csv(out_dir+"/df_err.csv", index = False)
        for i in range(len(cam_names)):
            df_drop = pd.DataFrame(output_repoerr_3d[i], columns=marker_names)
            df_drop.to_csv(out_dir+"/" + cam_names[i] + ".csv", index = False)


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