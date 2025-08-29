% close all;% clear all; clc;
% Channel 2 is used as the reference node, so subtract 1 from all channel
% numbers except channel 1.
data = ALLEEG.data';
Fs = ALLEEG.srate;
loops = size(ALLEEG.event);
loops = loops(2);
events = [];
sync = [];

for i = 2:loops
    temp = ALLEEG.event(i).type;
    temp(1:2) = [];
    temp = str2double(temp);
    events = [events,temp];
    temp2 = ALLEEG.event(i).latency;
    sync = [sync,temp2];
end

start_trigger = find(events < 8);
reach_trigger = start_trigger(1:2:end);
grasp_trigger = start_trigger(2:2:end);
% reach_trigger = reach_trigger(1:30);
reach_trigger = reach_trigger(121:180);
reach_trigger = reach_trigger(2:2:end);

reach_start = sync(reach_trigger);
for k = 1:31
    channels = k;%[5,6,7,10,11]; % These channels are over the motor cortex.
    reach_data = [];
    for i = 1:length(reach_start)
        data_points = (reach_start(i)-1499):(reach_start(i)+1000);
        temp = mean(data(data_points,channels),2);
        reach_data = [reach_data,temp];
    end
    reach_mean = mean(reach_data,2);
    
    [sr,fr,tr,psr] = spectrogram(reach_mean,500,480,1000,1000,"power",'yaxis');

    % Mean Normalize
    pretrigger_mean_r = mean(psr(:,1:floor(width(psr)/2)),2);
    psr = psr./pretrigger_mean_r;
    
    figure(2)
    subplot(2,1,1)
    hold on
    plot(linspace(-1,1,length(reach_mean(500:end))),reach_mean(500:end))
    xline(0,'-',{'Reach'});
    xlabel("Time (s)")
    ylabel("Amplitude (\muV)")
    title('Reach time series')
    
    set(gcf, 'Position', [100 100 1200 800])
    subplot(2,1,2)
    imagesc(tr-1.5+0.25,fr,psr)
    axis xy
    xlabel("Time (s)")
    ylabel("Frequency (Hz)")
    colorbar('southoutside')
    ylim([0 60])
    set(gca, ColorScale="log")
    % clim([0.01 1000]); % All data
    % clim([0.0001 0.3]); % Alpha
    title('Reach Spectrogram')
    colormap(jet)
    if k == 1
        channel_name = "1";
    else
        channel_name = int2str(k+1);
    end
    % titlestr = strcat('Bottle Motor NoFilt, All trials, Norm, M1avg');
    % figstr = strcat('Bottle_Motor_NoFilt_Alltrials_Norm_M1avg');
    titlestr = strcat('USB Imagery, All trials, Ch',channel_name);
    figstr = strcat('USB_Imagery_Alltrials_Norm_Ch',channel_name);
    sgtitle(titlestr)
    saveas(gcf,strcat('sessions\17_06\fig\',figstr,'.png'))
    close(2)
end