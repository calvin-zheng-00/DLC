close all; clear all; clc;
mode = 3;
src = "C:\Users\czhe0008\Documents\EEG\raw_data\reach\";
src_A = "C:\Users\czhe0008\Documents\EEG\raw_data\audio\";
maintemp = dir(fullfile(src,'*'));
mainfolder = setdiff({maintemp([maintemp.isdir]).name},{'.','..'});
audiotemp = dir(fullfile(src_A,'*'));
audiofolder = setdiff({audiotemp([audiotemp.isdir]).name},{'.','..'});
for k = 1:31
    %% Reach
    reach_data = [];
    for i = 1:numel(mainfolder)
        temp = cell2mat(mainfolder(i));
        % Determine which folders to grab data from
        if (temp(1:2) == "Ju") && (temp((end-5):end) == "tle_MI")
            subtemp = dir(fullfile(src,mainfolder{i},'*.mat'));
            subfolder = {subtemp(~[subtemp.isdir]).name};
            for j = 1:numel(subfolder)
                T = fullfile(src,mainfolder{i},subfolder{j});
                load(T)
                reach_data = [reach_data,mean(save_flag(1:end,k),2)]; %[6,7,8,11,12]
            end
        end
    end
    reach_data = reach_data(1001:end,:);

    % Noise removal
    amp = max(reach_data)-min(reach_data);
    amp_std = std(amp);
    amp_mean = mean(amp);
    amp_min = amp_mean - 2*amp_std;
    amp_max = amp_mean + 2*amp_std;
    amp_thresh = (amp<amp_max); %&(amp>amp_min);
    reach_data = reach_data(:,amp_thresh);
    reach_mean = mean(reach_data,2);

    %% Audio
    audio_data = [];
    for i = 1:numel(audiofolder)
        temp = cell2mat(audiofolder(i));
        % Determine which folders to grab data from
        if 1==1%(temp(1:4) == "July")
            subtemp = dir(fullfile(src_A,audiofolder{i},'*.mat'));
            subfolder = {subtemp(~[subtemp.isdir]).name};
            for j = 1:numel(subfolder)
                T = fullfile(src_A,audiofolder{i},subfolder{j});
                load(T)
                audio_data = [audio_data,mean(save_flag(1:end,k),2)]; %[6,7,8,11,12]
            end
        end
    end
    audio_data = audio_data(1001:end,:);

    % Noise removal
    amp = max(audio_data)-min(audio_data);
    amp_std = std(amp);
    amp_mean = mean(amp);
    amp_min = amp_mean - 2*amp_std;
    amp_max = amp_mean + 2*amp_std;
    amp_thresh = (amp<amp_max); %&(amp>amp_min);
    audio_data = audio_data(:,amp_thresh);
    audio_mean = mean(audio_data,2);


    reach_mean = reach_mean - audio_mean;
    [sr,fr,tr,psr] = spectrogram(reach_mean,500,480,1000,1000,"power",'yaxis');

    % Mean Normalize
    pretrigger_mean_r = mean(psr(:,1:50),2);
    psr = psr./pretrigger_mean_r;
    
    figure(2)
    subplot(2,1,1)
    hold on
    reach_mean_norm = reach_mean - mean(reach_mean);
    plot(linspace(-1,2,length(reach_mean_norm(500:end))),reach_mean_norm(500:end))
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
    if mode == 1
        titlestr = strcat('Bottle Motor, All trials, Ch',channel_name);
        figstr = strcat('Bottle_Motor_Alltrials_Norm_Ch',channel_name);
        save_folder = 'C:\Users\czhe0008\Documents\EEG\population_no_audio_all_sessions\Bottle_ME\';
    elseif mode == 2
        titlestr = strcat('Bottle Observation, All trials, Ch',channel_name);
        figstr = strcat('Bottle_Observation_Alltrials_Norm_Ch',channel_name);
        save_folder = 'C:\Users\czhe0008\Documents\EEG\population_no_audio_all_sessions\Bottle_AO\';
    elseif mode == 3
        titlestr = strcat('Bottle Imagery, All trials, Ch',channel_name);
        figstr = strcat('Bottle_Imagery_Alltrials_Norm_Ch',channel_name);
        save_folder = 'C:\Users\czhe0008\Documents\EEG\population_no_audio_all_sessions\Bottle_MI\';
    elseif mode == 4
        titlestr = strcat('USB Motor, All trials, Ch',channel_name);
        figstr = strcat('USB_Motor_Alltrials_Norm_Ch',channel_name);
        save_folder = 'C:\Users\czhe0008\Documents\EEG\population_no_audio_all_sessions\USB_ME\';
    elseif mode == 5
        titlestr = strcat('USB Observation, All trials, Ch',channel_name);
        figstr = strcat('USB_Observation_Alltrials_Norm_Ch',channel_name);
        save_folder = 'C:\Users\czhe0008\Documents\EEG\population_no_audio_all_sessions\USB_AO\';
    elseif mode == 6
        titlestr = strcat('USB Imagery, All trials, Ch',channel_name);
        figstr = strcat('USB_Imagery_Alltrials_Norm_Ch',channel_name);
        save_folder = 'C:\Users\czhe0008\Documents\EEG\population_no_audio_all_sessions\USB_MI\';
    end
    sgtitle(titlestr)
    if not(isfolder(save_folder))
        mkdir(save_folder)
    end
    saveas(gcf,strcat(save_folder,figstr,'.png'))
    close(2)
end