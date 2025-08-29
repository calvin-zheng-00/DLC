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

start_trigger = find(events < 8);
sync = sync(start_trigger);
noise_start = sync(1:2:end);
noise_end = sync(2:2:end);

% audio_start = sync(events == 1);
% audio_start1 = audio_start(1:29);
% audio_start2 = audio_start(30:end);
chan_list = [2,7,11,12,24];
% test_list = ["Nothing On, No Shield, Grounded Frame", "Nothing On, Shield, Grounded Frame",...
%     "Lights, Shield, Grounded Frame", "Nothing On, No Shield, No Grounding",...
%     "Nothing On, Shield, No Grounding", "Lights, Shield, No Grounding", "Lights, No Shield, Grounded Frame",...
%     "Lights, No Shield, No Grounding", "Lights, Fan, Shield, Grounded Frame", "Lights, Fan, Monitor, Shield, Grounded Frame",...
%     "Lights, Fan, Monitor, No Shield, Grounded Frame and EEG", "Lights, Fan, Monitor, Cameras, No Shield, Grounded Frame",...
%     "Lights, Fan, Monitor, Cameras, Shield, Grounded Frame", "Participant Moving Wires", "Observer Moving Wires"];
% test_list = ["10 blinks", "closed eyes", "Looking left then transition to right",...
%     "looking left 10 times","looking right 10 times","10 swallows","10 thumb presses",...
%     "10 arm lifts"];
range = [];
range2 = [];

for k = 1:31
    channels = k;%chan_list(k);

    noise_data = [];
    for i = 1:length(audio_start)
        data_points = (audio_start(i)-2499):(audio_start(i)+2000);
        temp = mean(data(data_points,channels),2);
        noise_data = [noise_data,temp];
    end
    noise_mean = mean(noise_data,2);
    [s,f,t,ps] = spectrogram(noise_mean,500,480,1000,1000,"power",'yaxis');

    % Mean Normalize
    pretrigger_mean = mean(ps(:,1:floor(width(ps)/2)),2);
    ps = ps./pretrigger_mean;

    figure(2)
    subplot(2,1,1)
    hold on
    plot(linspace(-1,1,length(noise_mean(500:end))),noise_mean(500:end))
    xline(0,'-',{'Reach'});
    xlabel("Time (s)")
    ylabel("Amplitude (\muV)")
    title('Reach time series')

    % upper_set = findpeaks(noise_data);
    % lower_set = findpeaks(noise_data*-1)*-1;
    % upper_lim = mean(maxk(upper_set,floor(length(upper_set)/2)));
    % lower_lim = mean(mink(lower_set,floor(length(lower_set)/2)));
    % upper_set = upper_set(upper_set>upper_lim);
    % lower_set = lower_set(lower_set<lower_lim);
    % upper_lim = mean(maxk(upper_set,floor(length(upper_set)/2)));
    % lower_lim = mean(mink(lower_set,floor(length(lower_set)/2)));
    % range = [range,upper_lim - lower_lim];
    % range = [range,min(noise_data)];
    % range2 = [range2,max(noise_data)];
    
    set(gcf, 'Position', [100 100 1200 800])
    subplot(2,1,2)
    imagesc(t-1.5+0.25,f,ps)
    axis xy
    xlabel("Time (s)")
    ylabel("Frequency (Hz)")
    colorbar('southoutside')
    ylim([4 80])
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
    titlestr = strcat('Audio Cue, All trials, Ch',channel_name);
    figstr = strcat('Audio_Cue_Alltrials_Norm_Ch',channel_name);
    sgtitle(titlestr)
    saveas(gcf,strcat('sessions\24_06\fig\audio cue\',figstr,'.png'))
    close(2)
end