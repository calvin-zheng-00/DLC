import logging
import multiprocessing
#from pathlib import Path
#from typing import Union

import numpy as np

logger = logging.getLogger(__name__)


def triangulate_3d_data(
    anipose_calibration_object,
    mediapipe_2d_data: np.ndarray,
    use_triangulate_ransac: bool = False,
):
    number_of_cameras = mediapipe_2d_data.shape[0]
    number_of_frames = mediapipe_2d_data.shape[1]
    number_of_tracked_points = mediapipe_2d_data.shape[2]
    number_of_spatial_dimensions = mediapipe_2d_data.shape[3]

    if not number_of_spatial_dimensions == 2:
        logger.error(
            f"This is supposed to be 2D data but, number_of_spatial_dimensions: {number_of_spatial_dimensions}"
        )
        raise Exception

    # reshape data to collapse across 'frames' so it becomes [number_of_cameras,
    # number_of_2d_points(numFrames*numPoints), XY]
    data2d_flat = mediapipe_2d_data.reshape(number_of_cameras, -1, 2)

    logger.info(
        f"Reconstructing 3d points from 2d points with shape: \n"
        f"number_of_cameras: {number_of_cameras},\n"
        f"number_of_frames: {number_of_frames}, \n"
        f"number_of_tracked_points: {number_of_tracked_points},\n"
        f"number_of_spatial_dimensions: {number_of_spatial_dimensions}"
    )

    if use_triangulate_ransac:
        logger.info("Using `triangulate_ransac` method")
        #data3d_flat = anipose_calibration_object.triangulate_ransac(data2d_flat, progress=True, kill_event=kill_event)
        data3d_flat = anipose_calibration_object.triangulate_ransac(data2d_flat, progress=True)
    else:
        logger.info("Using simple `triangulate` method ")
        #data3d_flat = anipose_calibration_object.triangulate(data2d_flat, progress=True, kill_event=kill_event)
        data3d_flat = anipose_calibration_object.triangulate(data2d_flat, progress=True)

    spatial_data3d_numFrames_numTrackedPoints_XYZ = data3d_flat.reshape(
        number_of_frames, number_of_tracked_points, 3
    )

    data3d_reprojectionError_flat = anipose_calibration_object.reprojection_error(data3d_flat, data2d_flat, mean=True)
    #New
    data3d_reprojectionError_full = anipose_calibration_object.reprojection_error(data3d_flat, data2d_flat, mean=False)
    reprojectionError_cam_frame_marker = np.linalg.norm(data3d_reprojectionError_full, axis=2).reshape(
        number_of_cameras, number_of_frames, number_of_tracked_points
    )
    #New

    reprojection_error_data3d_numFrames_numTrackedPoints = data3d_reprojectionError_flat.reshape(
        number_of_frames, number_of_tracked_points
    )

    return (
        spatial_data3d_numFrames_numTrackedPoints_XYZ,
        reprojection_error_data3d_numFrames_numTrackedPoints,
        reprojectionError_cam_frame_marker
    )
