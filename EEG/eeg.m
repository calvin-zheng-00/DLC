% close all; clear all; clc;
data = ALLEEG.data';
channels = 7; %[5,6,7,10,11];
data_window = 1:7000;
loops = size(ALLEEG.event);
loops = loops(2);
events = [];
sync = [];

for i = 2:loops
    % temp = str2double(ALLEEG.event(i).type(3:4));
    temp = ALLEEG.event(i).type;
    temp(1:2) = [];
    temp = str2double(temp);
    events = [events,temp];
    temp2 = ALLEEG.event(i).latency;
    sync = [sync,temp2];
end

start_trigger = find(events < 12);
% Reach start times: 222, 86, 46, 32, 68
% Reach end times: 668, 528, 432, 392, 448
% Reach triggers: 1, 17, 30, 44, 57
reach_start = [sync(start_trigger(1)+222),sync(start_trigger(17)+86),sync(start_trigger(30)+46),sync(start_trigger(44)+32),sync(start_trigger(57)+68)];
reach_end = [sync(start_trigger(1)+668),sync(start_trigger(17)+528),sync(start_trigger(30)+432),sync(start_trigger(44)+392),sync(start_trigger(57)+448)];
% Video start triggers: 68, 80, 95, 108, 120
% Video + imagery triggers: 132, 147, 159, 171, 182
% Video start and end: 2s, 20f to 9s, 20f
video_start = [sync(start_trigger(68)+200),sync(start_trigger(80)+200),sync(start_trigger(95)+200),sync(start_trigger(108)+200),sync(start_trigger(120)+200)];
video_end = [sync(start_trigger(68)+760),sync(start_trigger(80)+760),sync(start_trigger(95)+760),sync(start_trigger(108)+760),sync(start_trigger(120)+760)];
% Imagery triggers: 204, 231, 254, 279, 302
imagery_start = [sync(start_trigger(204)),sync(start_trigger(231)),sync(start_trigger(254)),sync(start_trigger(279)),sync(start_trigger(302))];
imagery_end = [sync(start_trigger(205)),sync(start_trigger(232)),sync(start_trigger(255)),sync(start_trigger(280)),sync(start_trigger(303))];
Fs = ALLEEG.srate;
% d = designfilt('bandstopiir','FilterOrder',2, ...
%                'HalfPowerFrequency1',49,'HalfPowerFrequency2',51, ...
%                'DesignMethod','butter','SampleRate',Fs);

%% Eye data
% data_points = data_window + reach_start;
% reach1 = data(data_points,channels)-mean(data(data_points,channels));
% wpass = [8,12];% * ((2*pi)/Fs); %8-12 for alpha, 18-26 for beta
% reach_filt = bandpass(reach1,wpass,Fs);

data_points1 = (reach_start(1)-1000):(reach_end(1)+1000);
data_points2 = (reach_start(2)-1000):(reach_end(2)+1000);
data_points3 = (reach_start(3)-1000):(reach_end(3)+1000);
data_points4 = (reach_start(4)-1000):(reach_end(4)+1000);
data_points5 = (reach_start(5)-1000):(reach_end(5)+1000);
temp1 = data(data_points1,channels)-mean(data(data_points1,channels));
temp2 = data(data_points2,channels)-mean(data(data_points2,channels));
temp3 = data(data_points3,channels)-mean(data(data_points3,channels));
temp4 = data(data_points4,channels)-mean(data(data_points4,channels));
temp5 = data(data_points5,channels)-mean(data(data_points5,channels));
resamp2 = resample(temp2,length(temp1),length(temp2));
resamp3 = resample(temp3,length(temp1),length(temp3));
resamp4 = resample(temp4,length(temp1),length(temp4));
resamp5 = resample(temp5,length(temp1),length(temp5));
temp = [temp1,resamp2,resamp3,resamp4,resamp5];
reach_mean = mean([temp1,resamp2,resamp3,resamp4,resamp5],2);
% reach = data(data_points_mean,7)-mean(data(data_points_mean,7));
figure;
s1 = spectrogram(reach_mean);
spectrogram(reach_mean,100,80,1000,1000,'yaxis')
set(gca, ColorScale="log")

data_points1 = (video_start(1)-1000):(video_end(1)+1000);
data_points2 = (video_start(2)-1000):(video_end(2)+1000);
data_points3 = (video_start(3)-1000):(video_end(3)+1000);
data_points4 = (video_start(4)-1000):(video_end(4)+1000);
data_points5 = (video_start(5)-1000):(video_end(5)+1000);
temp1 = data(data_points1,channels)-mean(data(data_points1,channels));
temp2 = data(data_points2,channels)-mean(data(data_points2,channels));
temp3 = data(data_points3,channels)-mean(data(data_points3,channels));
temp4 = data(data_points4,channels)-mean(data(data_points4,channels));
temp5 = data(data_points5,channels)-mean(data(data_points5,channels));
resamp2 = resample(temp2,length(temp1),length(temp2));
resamp3 = resample(temp3,length(temp1),length(temp3));
resamp4 = resample(temp4,length(temp1),length(temp4));
resamp5 = resample(temp5,length(temp1),length(temp5));
temp = [temp1,resamp2,resamp3,resamp4,resamp5];
video_mean = mean([temp1,resamp2,resamp3,resamp4,resamp5],2);
% reach2 = data(data_points,7)-mean(data(data_points,7));
figure;
s2 = spectrogram(video_mean);
spectrogram(video_mean,100,80,1000,1000,'yaxis')
set(gca, ColorScale="log")

data_points1 = (imagery_start(1)):(imagery_end(1));
data_points2 = (imagery_start(2)):(imagery_end(2));
data_points3 = (imagery_start(3)):(imagery_end(3));
data_points4 = (imagery_start(4)):(imagery_end(4));
data_points5 = (imagery_start(5)):(imagery_end(5));
temp1 = data(data_points1,channels)-mean(data(data_points1,channels));
temp2 = data(data_points2,channels)-mean(data(data_points2,channels));
temp3 = data(data_points3,channels)-mean(data(data_points3,channels));
temp4 = data(data_points4,channels)-mean(data(data_points4,channels));
temp5 = data(data_points5,channels)-mean(data(data_points5,channels));
resamp2 = resample(temp2,length(temp1),length(temp2));
resamp3 = resample(temp3,length(temp1),length(temp3));
resamp4 = resample(temp4,length(temp1),length(temp4));
resamp5 = resample(temp5,length(temp1),length(temp5));
temp = [temp1,resamp2,resamp3,resamp4,resamp5];
imagery_mean = mean([temp1,resamp2,resamp3,resamp4,resamp5],2);
% reach3 = data(data_points,7)-mean(data(data_points,7));
figure;
s3 = spectrogram(imagery_mean);
spectrogram(imagery_mean,100,80,1000,1000,'yaxis')
set(gca, ColorScale="log")

%132
%192-204, 205-217


%% Reach data
% data_points = data_window + reach_start;
% reach = data(data_points,channels)-mean(data(data_points,channels));
% % reach_filt = filtfilt(d,reach);
% wpass = [8,12];% * ((2*pi)/Fs); %8-12 for alpha, 18-26 for beta
% reach_filt = bandpass(reach,wpass,Fs);
% % reach = mean(reach,2);
% 
% reach_start = sync(temp(13)+12);
% data_points = data_window + reach_start;
% reach = data(data_points,channels)-mean(data(data_points,channels));
% wpass = [8,12];
% reach_filt2 = bandpass(reach,wpass,Fs);
% 
% reach_start = sync(temp(27)+9);
% data_points = data_window + reach_start;
% reach = data(data_points,channels)-mean(data(data_points,channels));
% wpass = [8,12];
% reach_filt3 = bandpass(reach,wpass,Fs);
% 
% reach_start = sync(temp(41)+7);
% data_points = data_window + reach_start;
% reach = data(data_points,channels)-mean(data(data_points,channels));
% wpass = [8,12];
% reach_filt4 = bandpass(reach,wpass,Fs);
% 
% reach_start = sync(temp(55)+6);
% data_points = data_window + reach_start;
% reach = data(data_points,channels)-mean(data(data_points,channels));
% wpass = [8,12];
% reach_filt5 = bandpass(reach,wpass,Fs);

%% Plotting
% figure
% plot(reach_filt)
% xline(1000)
% xline(5125)
% legend({'FC5','FC1','C3','CP5','CP1'},"Location",'northeastoutside') % 'FC5'
% box off
% xlabel('Time (ms)')
% ylabel('Amplitude (\muV)')
% title('EEG recording over the motor cortex from reach start (alpha)')
% fontsize(12,"points")
% set(gcf,'position',[200,200,1200,500])
% ylim([-10,10])

% data_points = data_window + eye_start(1);
% temp = data(data_points,channels)-mean(data(data_points,channels));
% temp = temp(:,1);
% fourier = fft(temp);
% 
% Fs = 1000;            % Sampling frequency                    
% T = 1/Fs;             % Sampling period       
% L = 1000;             % Length of signal
% t = (0:L-1)*T;        % Time vector
% figure
% plot(Fs/L*(0:L-1),abs(fourier))

%% Reach info
% For ball grab activity, trial 1, sync signal is sent 6 seconds and 10
% frames into the video (button is released after frame 15)
% Button is pressed during the event idx
% [1092,1093,1094,1095,1096,1097,1098,1099,1100,1101,1102,1103]
% These are 12 events, which aligns with the 6 frames (each frame has an up
% and down edge)
% The button signal overwrites on top of the camera sync signal, keeping
% the 40Hz, but changing the event ID.
%
% 7 seconds, frame 34 has the hand lifting off the table
% 8 seconds, frame 28 is when grasp is complete
%
% Fz (channel 2) is not included in the data. Most likely used as a
% reference node.
%
% For ball trial 1
% Button press at 6 seconds, 10 frames
% reach start: 7 seconds, 34 frames
% grasp moment: 8 seconds, 28 frames
% temp(1) + 44 (+76)
%
% For ball trial 2
% Button press at 4 seconds, 25 frames
% reach start: 5 seconds, 17 frames
% grasp moment: 6 seconds, 8 frames
% temp(13) + 12 (+43)
%
% For ball trial 3
% Button press at 3 seconds, 7 frames
% reach start: 3 seconds, 36 frames
% grasp moment: 4 seconds, 30 frames
% temp(27) + 9 (+43)
%
% For ball trial 4
% Button press at 5 seconds, 3 frames
% reach start: 5 seconds, 30 frames
% grasp moment: 6 seconds, 25 frames
% temp(41) + 7 (+42)
%
% For ball trial 5
% Button press at 5 seconds, 15 frames
% reach start: 6 seconds, 1 frames
% grasp moment: 6 seconds, 35 frames
% temp(55) + 6 (+40)