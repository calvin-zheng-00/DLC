data_set = 1;
loops_event = size(ALLEEG(data_set).event);
loops_event = loops_event(2);
event_time = [];
event_type = 'S  6';

flag = 0;
for i = 2:loops_event
    if strcmp(ALLEEG(data_set).event(i).type, event_type) %& (flag == 0) % 'R  8'
        temp = ALLEEG(data_set).event(i).latency;
        temp = round((temp*40)/1000);
        event_time = [event_time,temp];
    %     flag = 1;
    % elseif strcmp(ALLEEG(data_set).event(i).type, 'R  8')
    %     flag = 1;
    % else
    %     flag = 0;
    end
end
% syncFrame = 58*40 + 10; % 10_12
syncFrame = 35*40 + 17; % 11_12
% syncFrame = 904; % 13 10
startFrame = 1137 - syncFrame; % S  3 event 1
event_time = event_time - startFrame;

% vidObj = VideoReader("\\storage.erc.monash.edu.au\shares\MNHS-MoCap\Calvin\Ex2\eeg_raw_vid\11_12\e3v830e-20251211T143706-145917.avi");
% numFrames = 0;
% 
% while hasFrame(vidObj)
%     readFrame(vidObj); % Read the frame to advance the video object's state
%     numFrames = numFrames + 1;
% end
%numFrames = 53242;
% event_time = event_time - numFrames;
%e3v82e4, e3v830e
cam_list = ["e3v831f","e3v833c","e3v834c","e3v8337","e3v8364","e3v8380","e3v8389","e3v8393"];
src = '\\storage.erc.monash.edu.au\shares\MNHS-MoCap\Calvin\Ex2\eeg_raw_vid\11_12\';
vidtemp = dir(fullfile(src,'*.avi'));
vicfolder = {vidtemp(~[vidtemp.isdir]).name};
for k = 1:size(cam_list,2)
    cam_name = cam_list(k);
    camfolder = vicfolder(startsWith(vicfolder,cam_name));
    curr_vid = 1;
    v = VideoReader(char(fullfile(src,camfolder(curr_vid))));
    curr_event_time = event_time;
    
    for i = 1:size(curr_event_time,2)
        while (curr_event_time(i)/40) > v.Duration
            numFrames = floor(v.Duration*v.Framerate);
            curr_event_time = curr_event_time - numFrames;
            curr_vid = curr_vid + 1;
            v = VideoReader(char(fullfile(src,camfolder(curr_vid))));
        end
        v.CurrentTime=curr_event_time(i)/40;
        video_out_file = strcat("\\storage.erc.monash.edu.au\shares\MNHS-MoCap\Calvin\Ex2\eeg_raw_vid\11_12_edit2\",cam_name,"_",event_type,"_",int2str(i),".avi");
        vw = VideoWriter(video_out_file);
        vw.FrameRate = 40;
        open(vw);
        
        frameIdx = curr_event_time(i);
    
        next_vid_flag = true;
        
        while hasFrame(v)
            frame = readFrame(v);
            writeVideo(vw, frame);
            if frameIdx > (curr_event_time(i)+400)
                close(vw);
                next_vid_flag = false;
                break
            end
            frameIdx = frameIdx + 1;
        end
        if next_vid_flag
            curr_event_time = curr_event_time - frameIdx;
            frameIdx = 1;
            curr_vid = curr_vid + 1;
            v = VideoReader(char(fullfile(src,camfolder(curr_vid))));
            while hasFrame(v)
                frame = readFrame(v);
                writeVideo(vw, frame);
                if frameIdx > (curr_event_time(i)+400)
                    close(vw);
                    break
                end
                frameIdx = frameIdx + 1;
            end
        end
    end
end