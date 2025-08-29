% close all;% clear all; clc;
% Channel 2 is used as the reference node, so subtract 1 from all channel
% numbers except channel 1.
data = ALLEEG.data';
Fs = ALLEEG.srate;
% Left front: 1, 4
% Right motor: 24, 29
% Centre: 23
loops = size(ALLEEG.event);
loops = loops(2);
events = [];
sync = [];
main_event = readtable('C:\Users\czhe0008\Documents\EEG\29_05_2\events_bottle.csv');

for i = 2:loops
    temp = ALLEEG.event(i).type;
    temp(1:2) = [];
    temp = str2double(temp);
    events = [events,temp];
    temp2 = ALLEEG.event(i).latency;
    sync = [sync,temp2];
end

start_trigger = find(events > 1);
diff_trigger = diff(start_trigger);
trigger_idx = find(diff_trigger ~= 1) + 1;
sync_idx = start_trigger([1,trigger_idx]);
% sync_idx(29) = []; % Faulty trials in recording for events_bottle motor
% sync_idx([1,28]) = []; % Faulty trials in recording for usb video
% sync_idx([9,13,19,22,36]) = []; % Faulty trials in recording for Lucy Bottle video

% If motor
reach_start_idx = sync_idx - main_event.Var1' + main_event.Var2';
reach_start = sync(reach_start_idx);
grasp_start_idx = sync_idx - main_event.Var1' + main_event.Var3';
grasp_start = sync(grasp_start_idx);

% If video
% sync_idx = sync_idx(1:2:end);
% reach_start = sync(sync_idx)+76*25; % Bottle 68, usb 76
% grasp_start = sync(sync_idx)+105*25; % Bottle 101, usb 105

% If imagery
% sync_start = sync(sync_idx(1:2:end));
% sync_end = sync(sync_idx(2:2:end));
% sync_diff = floor(mean(sync_end-sync_start));

for k = 1:31
    channels = k;%[5,6,7,10,11]; % These channels are over the motor cortex.

    % For imagery
    % all_data = [];
    % for i = 1:length(sync_start)
    %     data_points = (sync_start(i)-2000):(sync_start(i)+sync_diff);
    %     temp = mean(data(data_points,channels),2);%-mean(data(data_points,channels));
    %     temp = temp - mean(temp);
    %     all_data = [all_data,temp];
    % end
    % all_mean = mean(all_data,2);
    % [s,f,t,ps] = spectrogram(all_mean,hamming(250),230,1000,1000,"power",'yaxis');
    % figure
    % subplot(2,1,1)
    % hold on
    % plot(linspace(-2,(sync_diff/1000),length(all_mean)),all_mean)
    % xline(0,'-',{'Imagine Start'});
    % xlabel("Time (s)")
    % ylabel("Amplitude (\muV)")
    % title('Imagination time series')
    % 
    % set(gcf, 'Position', [100 100 1200 800])
    % subplot(2,1,2)
    % waterplot(s,f,t-2);
    % xlim([0 50])
    % zlim([1 10^5])
    % imagesc(t-2,f,ps)
    % axis xy
    % xlabel("Time (s)")
    % ylabel("Frequency (Hz)")
    % colorbar('southoutside')
    % ylim([0 100])
    % set(gca, ColorScale="log")
    % title('Reach Spectrogram')

    % For everything else
    grasp_data = [];
    grasp_time = [];
    for i = 1:length(grasp_start)
        data_points = (grasp_start(i)-999):(grasp_start(i)+1000);
        temp = mean(data(data_points,channels),2);
        grasp_time = [grasp_time,temp];
        temp = temp - mean(temp);
        grasp_data = [grasp_data,temp];
    end
    grasp_mean = mean(grasp_data,2);
    grasp_time_mean = mean(grasp_time,2);
    reach_data = [];
    reach_time = [];
    for i = 1:length(reach_start)
        data_points = (reach_start(i)-999):(reach_start(i)+1000);
        temp = mean(data(data_points,channels),2);
        reach_time = [reach_time,temp];
        temp = temp - mean(temp);
        reach_data = [reach_data,temp];
    end
    reach_mean = mean(reach_data,2);
    reach_time_mean = mean(reach_time,2);
    
    [sg,fg,tg,psg] = spectrogram(grasp_time_mean,hamming(250),230,1000,1000,"power",'yaxis');
    [sr,fr,tr,psr] = spectrogram(reach_time_mean,hamming(250),230,1000,1000,"power",'yaxis');

    % Mean Normalize
    pretrigger_mean_g = mean(psg(:,1:floor(width(psg)/2)),2);
    psg = psg./pretrigger_mean_g;
    pretrigger_mean_r = mean(psr(:,1:floor(width(psr)/2)),2);
    psr = psr./pretrigger_mean_r;
    
    figure
    subplot(2,2,1)
    hold on
    plot(linspace(-1,1,length(reach_time_mean)),reach_time_mean)
    xline(0,'-',{'Reach'});
    xlabel("Time (s)")
    ylabel("Amplitude (\muV)")
    title('Reach time series')
    
    subplot(2,2,2)
    hold on
    plot(linspace(-1,1,length(grasp_time_mean)),grasp_time_mean)
    xline(0,'-',{'Grasp'});
    xlabel("Time (s)")
    ylabel("Amplitude (\muV)")
    title('Grasp time series')
    
    set(gcf, 'Position', [100 100 1200 800])
    subplot(2,2,3)
    imagesc(tr-1+0.115,fr,psr)
    axis xy
    xlabel("Time (s)")
    ylabel("Frequency (Hz)")
    colorbar('southoutside')
    ylim([0 60])
    set(gca, ColorScale="log")
    % clim([0.01 1000]); % All data
    % clim([0.0001 0.3]); % Alpha
    title('Reach Spectrogram')
    
    set(gcf, 'Position', [100 100 1200 800])
    subplot(2,2,4)
    imagesc(tg-1+0.115,fg,psg)
    axis xy
    xlabel("Time (s)")
    ylabel("Frequency (Hz)")
    colorbar('southoutside')
    ylim([0 60])
    set(gca, ColorScale="log")
    % clim([0.01 1000]); % All data
    % clim([0.0001 0.3]); % Alpha
    title('Grasp Spectrogram')
    colormap(jet)
    if k == 1
        channel_name = "1";
    else
        channel_name = int2str(k+1);
    end
    % titlestr = strcat('Bottle Motor NoFilt, All trials, Norm, M1avg');
    % figstr = strcat('Bottle_Motor_NoFilt_Alltrials_Norm_M1avg');
    titlestr = strcat('Bottle Motor, All trials, Ch',channel_name);
    figstr = strcat('Bottle_Motor_Alltrials_Norm_Ch',channel_name);
    sgtitle(titlestr)
    saveas(gcf,strcat('29_05_2\fig\',figstr,'.png'))
end