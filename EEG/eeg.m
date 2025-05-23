% close all; clear all; clc;
% Channel 2 is used as the reference node, so subtract 1 from all channel
% numbers except channel 1.
data = ALLEEG.data';
% data = data - mean(data);
Fs = ALLEEG.srate;
% Left front: 1, 4
% Right motor: 24, 29
% Centre: 23
channels = [15,16,17];%[7,15,16,17]; % These channels are over the optical cortex.
loops = size(ALLEEG.event);
loops = loops(2);
events = [];
sync = [];
main_event = readtable('C:\Users\czhe0008\Documents\EEG\22_05\events_bottle.csv');

for i = 2:loops
    temp = ALLEEG.event(i).type;
    temp(1:2) = [];
    temp = str2double(temp);
    events = [events,temp];
    temp2 = ALLEEG.event(i).latency;
    sync = [sync,temp2];
end

start_trigger = find(events < 12);
light_trigger = find(events == 12);
diff_trigger = diff(start_trigger);
trigger_idx = find(diff_trigger ~= 1) + 1;
sync_idx = start_trigger([1,trigger_idx]);
light_trigger = light_trigger((light_trigger>start_trigger(1) & light_trigger<start_trigger(2)));
light_trigger(end-10:end) = [];
light_1sec = sync(light_trigger);

% data_points = (light_1sec(1)-999):(light_1sec(5)+1000);
% light_data = mean(data(data_points,channels),2);
% light_norm = light_data - mean(light_data);
light_data = [];
light_time = [];
for i = 1:(length(light_1sec)-1)
    data_points = (light_1sec(i)-500):(light_1sec(i+1)-500);
    temp = mean(data(data_points,channels),2);
    if i ~= 1
        if height(temp)>height(light_time)
            temp = temp(1:height(light_time));
            light_time = cat(2,light_time,temp);
            temp = temp - mean(temp);
            light_data = cat(2,light_data,temp);
        else
            light_time = light_time(1:height(temp),:);
            light_data = light_data(1:height(temp),:);
            light_time = cat(2,light_time,temp);
            temp = temp - mean(temp);
            light_data = cat(2,light_data,temp);
        end
    else
        light_time = cat(2,light_time,temp);
        temp = temp - mean(temp);
        light_data = cat(2,light_data,temp);
    end
end
light_mean = mean(light_data,2);
light_time_mean = mean(light_time,2);

[s,f,t,ps] = spectrogram(light_mean,hamming(250),230,1000,1000,"power",'yaxis');

figure
subplot(2,1,1)
hold on
plot((1:length(light_time_mean))/1000,light_time_mean)
xline(0.5,'-',{'strobe'});
% xline((light_5sec(2)-light_5sec(1))/1000,'-',{'strobe'});
% xline((light_5sec(3)-light_5sec(1))/1000,'-',{'strobe'});
% xline((light_5sec(4)-light_5sec(1))/1000,'-',{'strobe'});
% xline((light_5sec(5)-light_5sec(1))/1000,'-',{'strobe'});
xlim([0,length(light_time_mean)/1000])
xlabel("Time (s)")
ylabel("Amplitude (\muV)")
title('light time series')

set(gcf, 'Position', [100 100 1200 800])
subplot(2,1,2)
imagesc(t,f,ps)
axis xy
xlabel("Time (s)")
ylabel("Frequency (Hz)")
colorbar('southoutside')
ylim([0 60])
set(gca, ColorScale="log")
% clim([0.01 1000]); % All data
% clim([0.0001 0.3]); % Alpha
title('light Spectrogram')
colormap(jet)

% sgtitle('Bottle Motor NoFilt Ch7')
sgtitle('Light test ChO')