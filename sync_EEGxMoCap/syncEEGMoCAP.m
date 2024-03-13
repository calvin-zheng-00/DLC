%% EEG & MoCap synchronisation

% Set up recording object
% You should only need to change the date-time part of the filepath
recordObj.path = 'Open Ephys\2024-03-06_13-16-52\Record Node 102\experiment1\recording1\continuous\Acquisition_Board-100.Rhythm Data\continuous.dat';
recordObj.channels = 64; % 2 headstages, 32 electrodes each
recordObj.data = fread(fopen(recordObj.path,'rb'),[recordObj.channels, Inf], 'int16');
recordObj.fs = 30000; % 30kHz sampling rate

% Set up timestamp object
timesObj.path = 'Open Ephys\2024-03-06_13-16-52\Record Node 102\experiment1\recording1\continuous\Acquisition_Board-100.Rhythm Data\timestamps.npy';
timesObj.data = readNPY(timesObj.path);

% clean EEG signal - advice from Hannah: downsample then lowpass
% Downsample first
fd = 300; % downsample by factor of 100 - arbitrary choice for now

recordObj.downSamp = resample_array(recordObj.data, fd, recordObj.fs);

indices = ceil(linspace(1,size(recordObj.data,2),size(recordObj.downSamp,2)));

for i = 1:size(recordObj.downSamp,1)
   recordObj.reSamp(i,indices) = recordObj.downSamp(i,:);
end

% Lowpass filter
f_cutoff = 15; 
for k = 1:recordObj.channels
    recordObj.filt(k,:) = lowpass(recordObj.reSamp(k,:),f_cutoff,fd);
end

% LATER PROBLEM: may need to shorten to get uniform recording length if comparing sessions

%% plotting 
% LATER PROBLEM: time shift 
temp = ones(size(timesObj.data,1),1);

for i = 1:recordObj.channels
    figure ()
    plot(recordObj.filt(i,1:1000));
    hold on
    plot(timesObj.data(1:1000),temp(1:1000),'ro');
end


