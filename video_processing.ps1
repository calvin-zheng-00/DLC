Get-ChildItem "C:\Users\czhe0008\Documents\DLCprojects\openpose\openpose\examples\media\test" -Filter *.avi |
Foreach-Object {
    $video = "C:\Users\czhe0008\Documents\DLCprojects\openpose\openpose\examples\media\test\"+$_.Name
    $output = "C:\Users\czhe0008\Documents\DLCprojects\openpose\openpose\examples\media\test2\"+$_.BaseName
    & "bin\OpenPoseDemo.exe" --video $video --hand --hand_render 1 --net_resolution -1x64 --display 0 --render_pose 0 --hand_render 0 --write_json $output
}