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

for i = 2:loops
    temp = ALLEEG.event(i).type;
    temp(1:2) = [];
    temp = str2double(temp);
    events = [events,temp];
    temp2 = ALLEEG.event(i).latency;
    sync = [sync,temp2];
end

noise_start = sync(1:2:end);
noise_end = sync(2:2:end);
first_set_start = sync(events == 7);
first_set_end = first_set_start(2:2:end);
first_set_end(end) = [];
first_set_start = first_set_start(1:2:end);
first_set_start(end) = [];
chan_list = [2,7,11,12,24];
% test_list = ["Nothing On, No Shield, Grounded Frame", "Nothing On, Shield, Grounded Frame",...
%     "Lights, Shield, Grounded Frame", "Nothing On, No Shield, No Grounding",...
%     "Nothing On, Shield, No Grounding", "Lights, Shield, No Grounding", "Lights, No Shield, Grounded Frame",...
%     "Lights, No Shield, No Grounding", "Lights, Fan, Shield, Grounded Frame", "Lights, Fan, Monitor, Shield, Grounded Frame",...
%     "Lights, Fan, Monitor, No Shield, Grounded Frame and EEG", "Lights, Fan, Monitor, Cameras, No Shield, Grounded Frame",...
%     "Lights, Fan, Monitor, Cameras, Shield, Grounded Frame", "Participant Moving Wires", "Observer Moving Wires"];
test_list = ["10 blinks", "closed eyes", "Looking left then transition to right",...
    "looking left 10 times","looking right 10 times","10 swallows","10 thumb presses",...
    "10 arm lifts"];
range = [];
range2 = [];

for k = 1:1
    channels = 11;%chan_list(k);

    for i = 1:length(first_set_start)
        data_points = (first_set_start(i)):(first_set_end(i));
        noise_data = mean(data(data_points,channels),2);
        noise_data = noise_data - mean(noise_data);
        [s,f,t,ps] = spectrogram(noise_data,500,450,1000,1000,"power",'yaxis');
    
        % Mean Normalize
        % pretrigger_mean = mean(ps(:,1:floor(width(ps)/2)),2);
        % ps = ps./pretrigger_mean;

        figure
        subplot(2,1,1)
        hold on
        temp = (1:length(noise_data))/1000;
        plot(temp,noise_data)
        xlabel("Time (s)")
        ylabel("Amplitude (\muV)")
        title('Noise time series')
        xlim([0 temp(end)])

        % upper_set = findpeaks(noise_data);
        % lower_set = findpeaks(noise_data*-1)*-1;
        % upper_lim = mean(maxk(upper_set,floor(length(upper_set)/2)));
        % lower_lim = mean(mink(lower_set,floor(length(lower_set)/2)));
        % upper_set = upper_set(upper_set>upper_lim);
        % lower_set = lower_set(lower_set<lower_lim);
        % upper_lim = mean(maxk(upper_set,floor(length(upper_set)/2)));
        % lower_lim = mean(mink(lower_set,floor(length(lower_set)/2)));
        % range = [range,upper_lim - lower_lim];
        range = [range,min(noise_data)];
        range2 = [range2,max(noise_data)];
        
        set(gcf, 'Position', [100 100 1200 800])
        subplot(2,1,2)
        imagesc(t+0.115,f,ps)
        axis xy
        xlabel("Time (s)")
        ylabel("Frequency (Hz)")
        colorbar('southoutside')
        ylim([0 60])
        set(gca, ColorScale="log")
        % clim([0.01 1000]); % All data
        % clim([0.0001 0.3]); % Alpha
        title('Noise Spectrogram')

        colormap(jet)
        channel_name = int2str(channels+1);
        titlestr = strcat("Noise test, ", test_list(i),", ch",channel_name);
        figstr = strcat("Noise_test_", test_list(i),"_ch",channel_name);
        sgtitle(titlestr)
        saveas(gcf,strcat('sessions\17_06\fig\',figstr,'.png'))
        close(2)
    end
    
    % titlestr = strcat('Bottle Motor NoFilt, All trials, Norm, M1avg');
    % figstr = strcat('Bottle_Motor_NoFilt_Alltrials_Norm_M1avg');
end