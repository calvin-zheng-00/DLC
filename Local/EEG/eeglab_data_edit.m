% close all; clear all; clc;
% Channel 2 is used as the reference node, so subtract 1 from all channel
% numbers except channel 1.

data_set = 1;
loops_event = size(ALLEEG(data_set).event);
loops_event = loops_event(2);

second_check = 0;
for i = 2:loops_event
    if ALLEEG(data_set).event(i).code == "Stimulus"
        % Assuming stimulus come in pairs
        if second_check == 0
            second_check = 1;
        else
            second_check = 0;
            ALLEEG(data_set).event(i).type = 'S  2';
            ALLEEG(data_set).urevent(i).type = 'S  2';
        end
    end
end

count_events = 0;
for i = 2:loops_event
    if ALLEEG(data_set).event(i).type == 'S  1'
        % if (count_events >=30) && (count_events < 60)
        if (count_events <30)
            ALLEEG(data_set).event(i).type = 'S  3';
            count_events = count_events + 1;
        elseif (count_events >=90) && (count_events < 120)
            ALLEEG(data_set).event(i).type = 'S  4';
            count_events = count_events + 1;
        else
            count_events = count_events + 1;
        end
    end
end

% events = [];
% sync = [];
% for i = 2:loops_event
%     temp = ALLEEG(data_set).event(i).type;
%     temp(1:2) = [];
%     temp = str2double(temp);
%     events = [events,temp];
%     temp2 = ALLEEG(data_set).event(i).latency;
%     sync = [sync,temp2];
% end

%% If new setup
% data = ALLEEG(data_set).data';
% reach_trigger = find(events == 1);
% return_trigger = find(events == 2);
% reach_start = sync(reach_trigger);
% 
% for k = 1:31
%     save_flag = [];
%     for i = 1:length(reach_start)
%         data_points = (reach_start(i)-999):(reach_start(i)+2000);
%         if (reach_start(i)+2000)>size(data,1)
%             continue
%         end
%         save_flag = [save_flag,data(data_points,k)];
% 
%     end
%     save_folder = 'C:\Users\czhe0008\Documents\EEG\raw_data\test2\per_channel\';
%     if not(isfolder(save_folder))
%         mkdir(save_folder)
%     end
%     save_loc = strcat(save_folder,'channel',int2str(k),'.mat');
% 
%     save(save_loc, "save_flag")
% end