% close all;% clear all; clc;
% Channel 2 is used as the reference node, so subtract 1 from all channel
% numbers except channel 1.

% ALLEEG(1).event(4).type = 'S  2';

data = ALLEEG(1).data';
% Fs = ALLEEG.srate;
loops = size(ALLEEG(1).event);
loops = loops(2);
events = [];
sync = [];
% main_event = readtable('C:\Users\czhe0008\Documents\EEG\sessions\10_07\events_usb.csv');

for i = 2:loops
    temp = ALLEEG(1).event(i).type;
    temp(1:2) = [];
    temp = str2double(temp);
    events = [events,temp];
    temp2 = ALLEEG(1).event(i).latency;
    sync = [sync,temp2];
end

%% If new setup
% audio_start = sync(events == 1);
start_trigger = find(events == 1);
% test_trigger = find(events == 2);
% reach_trigger = start_trigger(423:452);
reach_trigger = start_trigger(1:2:end);
% reach_trigger = reach_trigger(1:30);
% reach_trigger = reach_trigger(391:450);
% reach_trigger = reach_trigger(2:2:end);
reach_start = sync(reach_trigger);

% grasp_trigger = start_trigger(2:2:end);
% grasp_trigger = grasp_trigger(1:30);
% grasp_trigger = grasp_trigger(31:90);
% grasp_trigger = grasp_trigger(2:2:end);
% grasp_start = sync(grasp_trigger);
% start_trigger = find(events == 1);
% test_trigger = start_trigger(181:240);
% test_trigger = test_trigger(2:2:end);
% test_start = sync(test_trigger);

% If old setup
% start_trigger = find(events == 7);
% diff_trigger = diff(start_trigger);
% trigger_idx = find(diff_trigger ~= 1) + 1;
% sync_idx = start_trigger([1,trigger_idx]);
% sync_idx(29) = []; % Faulty trials in recording for 22_05 events_bottle motor
% sync_idx([1,28]) = []; % Faulty trials in recording for 22_05 usb video
% sync_idx([6,24]) = []; % Faulty trials in recording for 22_05 events_bottle motor 2
% sync_idx([9,13,19,22,36]) = []; % Faulty trials in recording for Lucy Bottle motor

% If motor
% reach_start_idx = sync_idx - main_event.Var1' + main_event.Var2';
% reach_start = sync(reach_start_idx);
% grasp_start_idx = sync_idx - main_event.Var1' + main_event.Var3';
% grasp_start = sync(grasp_start_idx);

% If video
% reach_start = sync(sync_idx)+76*25; % Bottle 68, usb 76
% grasp_start = sync(sync_idx)+105*25; % Bottle 101, usb 105

for k = 1:31
    save_flag = [];
    for i = 1:length(reach_start)
        data_points = (reach_start(i)-999):(reach_start(i)+2000);
        if (reach_start(i)+2000)>size(data,1)
            continue
        end
        save_flag = [save_flag,data(data_points,k)];
        
    end
    save_folder = 'C:\Users\czhe0008\Documents\EEG\raw_data\test2\per_channel\';
    if not(isfolder(save_folder))
        mkdir(save_folder)
    end
    save_loc = strcat(save_folder,'channel',int2str(k),'.mat');

    save(save_loc, "save_flag")
end