close all;% clear all; clc;
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
main_event = readtable('C:\Users\czhe0008\Documents\EEG\29_05_1\events_bottle.csv');

for i = 2:loops
    temp = ALLEEG.event(i).type;
    temp(1:2) = [];
    temp = str2double(temp);
    events = [events,temp];
    temp2 = ALLEEG.event(i).latency;
    sync = [sync,temp2];
end

start_trigger = find(events < 12);
diff_trigger = diff(start_trigger);
trigger_idx = find(diff_trigger ~= 1) + 1;
sync_idx = start_trigger([1,trigger_idx]);
sync_idx([9,13,19,22,36]) = []; % Faulty trials in recording for Lucy Bottle video
% sync_idx(29) = []; % Faulty trials in recording for events_bottle motor
% sync_idx(1) = []; % Faulty trials in recording for usb video
% sync_idx(28) = []; % Faulty trials in recording for usb video

% If motor
reach_start_idx = sync_idx - main_event.Var1' + main_event.Var2';
reach_start = sync(reach_start_idx);
grasp_start_idx = sync_idx - main_event.Var1' + main_event.Var3';
grasp_start = sync(grasp_start_idx);

% reach_start(1) = []; % Faulty trials in recording for usb video
% reach_start(28) = []; % Faulty trials in recording for usb video
% grasp_start(1) = []; % Faulty trials in recording for usb video
% grasp_start(28) = []; % Faulty trials in recording for usb video

% If video
% sync_idx = sync_idx(1:2:end);
% reach_start = sync(sync_idx+76); % Bottle 68, usb 76
% grasp_start = sync(sync_idx+105); % Bottle 101, usb 105

for i = 1:length(grasp_start)
    channels = 11;%[5,6,7,10,11]; % These channels are over the motor cortex.
    % For everything else
    figure
    data_points = (grasp_start(i)-499):(grasp_start(i)+500);
    temp = mean(data(data_points,channels),2);
    [sg,fg,tg,psg] = spectrogram(temp,hamming(250),230,1000,1000,"power",'yaxis');
    pretrigger_mean_g = mean(psg(:,1:floor(width(psg)/2)),2);
    psg = psg./pretrigger_mean_g;
    subplot(2,2,2)
    hold on
    plot(linspace(-0.5,0.5,length(temp)),temp)
    xline(0,'-',{'Grasp'});
    xlabel("Time (s)")
    ylabel("Amplitude (\muV)")
    title('Grasp time series')
    subplot(2,2,4)
    imagesc(tg-0.5+0.115,fg,psg)
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
    data_points = (reach_start(i)-499):(reach_start(i)+500);
    temp = mean(data(data_points,channels),2);
    [sr,fr,tr,psr] = spectrogram(temp,hamming(250),230,1000,1000,"power",'yaxis');
    pretrigger_mean_r = mean(psr(:,1:floor(width(psr)/2)),2);
    psr = psr./pretrigger_mean_r;
    subplot(2,2,1)
    hold on
    plot(linspace(-0.5,0.5,length(temp)),temp)
    xline(0,'-',{'Reach'});
    xlabel("Time (s)")
    ylabel("Amplitude (\muV)")
    title('Reach time series')
    subplot(2,2,3)
    imagesc(tr-0.5+0.115,fr,psr)
    axis xy
    xlabel("Time (s)")
    ylabel("Frequency (Hz)")
    colorbar('southoutside')
    ylim([0 60])
    set(gca, ColorScale="log")
    % clim([0.01 1000]); % All data
    % clim([0.0001 0.3]); % Alpha
    title('Reach Spectrogram')

    titlestr = strcat('Bottle Motor, Ch12, trial',int2str(i));
    figstr = strcat('Bottle_Motor_Ch12_Trial',int2str(i));
    sgtitle(titlestr)
    set(gcf, 'Position', [100 100 1200 800])
    saveas(gcf,strcat('29_05_1\fig\',figstr,'.png'))
end