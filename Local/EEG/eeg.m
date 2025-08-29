%close all; %clear all; clc;
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
    
for k = 1:1
    channels = [15,16,17]; % These channels are over the optical cortex.
    light_data = [];
    light_time = [];
    for i = 1:(length(light_1sec)-1)
        data_points = (light_1sec(i)-750):(light_1sec(i+1)-500);
        temp = mean(data(data_points,channels),2);
        if i ~= 1
            if height(temp)>height(light_time)
                temp = temp(1:height(light_time));
                light_time = cat(2,light_time,temp);
                % temp = temp - mean(temp);
                % light_data = cat(2,light_data,temp);
            else
                light_time = light_time(1:height(temp),:);
                % light_data = light_data(1:height(temp),:);
                light_time = cat(2,light_time,temp);
                % temp = temp - mean(temp);
                % light_data = cat(2,light_data,temp);
            end
        else
            light_time = cat(2,light_time,temp);
            % temp = temp - mean(temp);
            % light_data = cat(2,light_data,temp);
        end
    end
    % light_mean = mean(light_data,2);
    light_time_mean = mean(light_time,2);
    
    [s,f,t,ps] = spectrogram(light_time_mean,hamming(250),230,1000,1000,"power",'yaxis');

    % Mean Normalize
    pretrigger_mean = mean(ps(:,1:floor(width(ps)/2)),2);
    ps = ps./pretrigger_mean;
    % ps = max(ps,0.000001);
    
    figure
    subplot(2,1,1)
    hold on
    plot(((1:length(light_time_mean(250:end)))/1000)-0.5,light_time_mean(250:end))
    xline(0,'-',{'strobe'});
    % xline((light_5sec(2)-light_5sec(1))/1000,'-',{'strobe'});
    % xline((light_5sec(3)-light_5sec(1))/1000,'-',{'strobe'});
    % xline((light_5sec(4)-light_5sec(1))/1000,'-',{'strobe'});
    % xline((light_5sec(5)-light_5sec(1))/1000,'-',{'strobe'});
    xlim([-0.5,0.5])
    xlabel("Time (s)")
    ylabel("Amplitude (\muV)")
    title('Strobe Time Series')

    set(gcf, 'Position', [100 100 1200 800])
    subplot(2,1,2)
    imagesc(t-0.75+0.125,f,ps)
    axis xy
    xlabel("Time (s)")
    ylabel("Frequency (Hz)")
    colorbar('southoutside')
    ylim([0 60])
    set(gca, 'YTick', 0:20:60)
    set(gca, ColorScale="log")
    % clim([-2.6 1.4]); % All data
    % clim([0.0001 0.3]); % Alpha
    title('Strobe Spectrogram')
    colormap(jet)

    % figure
    % plot(((1:length(light_time_mean(500:1000)))/1000)-0.25,light_time_mean(500:1000))
    % xline(0,'-',{'event'});
    % xlim([-0.25 0.25]);
    % xlabel("Time (s)")
    % ylabel("Amplitude (\muV)")

    fontsize(16,"points")
    if k == 1
        channel_name = "O";
    else
        channel_name = int2str(k+1);
    end
    titlestr = strcat('Strobe Test, Avg All Trials, Ch',channel_name);
    title(titlestr)
    figstr = strcat('Strobe_test_NoFilt_Alltrials_Norm_Ch',channel_name);
    saveas(gcf,strcat('22_05\fig\',figstr,'.png'))
end