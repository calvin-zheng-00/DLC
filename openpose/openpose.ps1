$src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\test_data\video\p18\mouse'
$dest = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\test_data\json\p18_mouse\'
$dest2 = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\Zac\tracked_videos\'

Set-Location "C:\Users\czhe0008\Documents\DLCprojects\openpose\openpose"

Get-ChildItem $src -Filter *.avi | 
Foreach-Object {
    $argvid = $_.FullName
    $outname = $_.BaseName
    $argjson = "$($dest)$($outname)"
    $outvid = "$($dest2)$($outname).avi"
    & bin\OpenPoseDemo.exe --video $argvid --hand --hand_render 1 --net_resolution -1x64 --write_json $argjson
}