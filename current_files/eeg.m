% close all; clear all; clc;
data = ALLEEG.data';
channels = [6,7,10,11]; %[5,6,7,10,11];
data_points = 1:2000;
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

temp = find(events < 12);
reach_start = sync(temp(1)+44); % I am graphing 1 second around the start of the reach (starts 0.5 seconds before)
grasp_start = sync(temp(1)+76); % I am graphing 1 second around the moment of grasp (starts 0.5 seconds before)
Fs = ALLEEG.srate;
% d = designfilt('bandstopiir','FilterOrder',2, ...
%                'HalfPowerFrequency1',49,'HalfPowerFrequency2',51, ...
%                'DesignMethod','butter','SampleRate',Fs);


data_points = data_points + reach_start;
reach = data(data_points,channels)-mean(data(data_points,channels));
% reach_filt = filtfilt(d,reach);
wpass = [8,12];% * ((2*pi)/Fs); %8-12 for alpha, 18-26 for beta
reach_filt = bandpass(reach,wpass,Fs);
% reach = mean(reach,2);

reach_start = sync(temp(13)+12);
data_points = (1:2000) + reach_start;
reach = data(data_points,channels)-mean(data(data_points,channels));
wpass = [8,12];
reach_filt2 = bandpass(reach,wpass,Fs);

reach_start = sync(temp(27)+9);
data_points = (1:2000) + reach_start;
reach = data(data_points,channels)-mean(data(data_points,channels));
wpass = [8,12];
reach_filt3 = bandpass(reach,wpass,Fs);

reach_start = sync(temp(41)+7);
data_points = (1:2000) + reach_start;
reach = data(data_points,channels)-mean(data(data_points,channels));
wpass = [8,12];
reach_filt4 = bandpass(reach,wpass,Fs);

reach_start = sync(temp(55)+6);
data_points = (1:2000) + reach_start;
reach = data(data_points,channels)-mean(data(data_points,channels));
wpass = [8,12];
reach_filt5 = bandpass(reach,wpass,Fs);

figure
plot(reach_filt5)
xline(500)
legend({'FC1','C3','CP5','CP1'},"Location",'northeastoutside') % 'FC5'
box off
xlabel('Time (ms)')
ylabel('Amplitude (\muV)')
title('EEG recording over the motor cortex from reach start (alpha)')
fontsize(12,"points")
set(gcf,'position',[200,200,1200,500])
ylim([-10,10])

% temp = data(data_points+data_start,channels)-mean(data(data_points+data_start,channels));
% temp = temp(:,1);
% fourier = fft(temp);
% 
% Fs = 1000;            % Sampling frequency                    
% T = 1/Fs;             % Sampling period       
% L = 1000;             % Length of signal
% t = (0:L-1)*T;        % Time vector
% figure
% plot(Fs/L*(0:L-1),abs(fourier))

% For ball grab activity, trial 1, sync signal is sent 6 seconds and 10
% frames into the video (button is released after frame 15)
% Button is pressed during the event idx
% [1092,1093,1094,1095,1096,1097,1098,1099,1100,1101,1102,1103]
% These are 12 events, which aligns with the 6 frames (each frame has an up
% and down edge)
% The button signal overwrites on top of the camera sync signal, keeping
% the 40Hz, but changing the event ID.

% 7 seconds, frame 34 has the hand lifting off the table
% 8 seconds, frame 28 is when grasp is complete

% Fz (channel 2) is not included in the data. Most likely used as a
% reference node.

% For ball trial 1
% Button press at 6 seconds, 10 frames
% reach start: 7 seconds, 34 frames
% grasp moment: 8 seconds, 28 frames
% temp(1) + 44 (+76)

% For ball trial 2
% Button press at 4 seconds, 25 frames
% reach start: 5 seconds, 17 frames
% grasp moment: 6 seconds, 8 frames
% temp(13) + 12 (+43)

% For ball trial 3
% Button press at 3 seconds, 7 frames
% reach start: 3 seconds, 36 frames
% grasp moment: 4 seconds, 30 frames
% temp(27) + 9 (+43)

% For ball trial 4
% Button press at 5 seconds, 3 frames
% reach start: 5 seconds, 30 frames
% grasp moment: 6 seconds, 25 frames
% temp(41) + 7 (+42)

% For ball trial 5
% Button press at 5 seconds, 15 frames
% reach start: 6 seconds, 1 frames
% grasp moment: 6 seconds, 35 frames
% temp(55) + 6 (+40)