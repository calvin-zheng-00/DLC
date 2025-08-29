% close all; clear all; clc;
mode = 4;
src = "C:\Users\czhe0008\Documents\EEG\raw_data\test2\";
maintemp = dir(fullfile(src,'*'));
mainfolder = setdiff({maintemp([maintemp.isdir]).name},{'.','..'});
for k = 1:31
    reach_data = [];
    for i = 1:numel(mainfolder)
        temp = cell2mat(mainfolder(i));
        % Determine which folders to grab data from
        % if (temp(1:7) == "July_10") && (temp((end-5):end) == "USB_ME") %||(temp(1:7) == "June_12")||(temp(1:7) == "July_02")) 
        if (temp == "Aug_19_pinch_40")%(temp(1:2) == "Ju") && (temp((end-5):end) == "tle_ME") % temp(7:end) == "Right_Grasp"
            subtemp = dir(fullfile(src,mainfolder{i},'*.mat'));
            subfolder = {subtemp(~[subtemp.isdir]).name};
            for j = 1:numel(subfolder)
                T = fullfile(src,mainfolder{i},subfolder{j});
                load(T)
                % reach_data_prep = highpass(save_flag(1:end,k),1,1000);
                reach_data_prep = save_flag(1:end,k);
                reach_data = [reach_data,reach_data_prep]; %[6,7,8,11,12]
            end
        end
    end
    reach_data = reach_data(1001:(end-2000),:);
    %% syncing data based on moment of grasp
    % 4000, 3725, 3725, 3825, 3625, 3650, 3900, 3650, 3625, 3775
    % reach_data2 = [reach_data(3001:5000,2),reach_data(2726:4725,3),reach_data(2726:4725,4),...
    %     reach_data(2826:4825,5),reach_data(2626:4625,6),reach_data(2651:4650,7),...
    %     reach_data(2901:4900,8),reach_data(2651:4650,9),reach_data(2626:4625,10),...
    %     reach_data(2776:4775,11)];

    %% Noise removal
    % amp = max(reach_data)-min(reach_data);
    % amp_std = std(amp);
    % amp_mean = mean(amp);
    % amp_min = amp_mean - 2*amp_std;
    % amp_max = amp_mean + 2*amp_std;
    % amp_thresh = (amp<amp_max); %&(amp>amp_min);
    % reach_data = reach_data(:,amp_thresh);

    %% Trial rejection
    % trial_mean = mean(reach_data);
    % reach_norm = reach_data - trial_mean;
    % amp = [];
    % for i = 1:width(reach_norm)
    %     [pks,locs] = findpeaks(reach_norm(1:2500,i));
    %     [negpks,neglocs] = findpeaks(-reach_norm(1:2500,i));
    %     amp = [amp,mean(pks) - mean(-negpks)];
    % end
    % 
    % lower_bound = pks > reach_min;
    % lower_bound = sum(lower_bound,1);
    % lower_bound = lower_bound==length(reach_data);
    % upper_bound = pks < reach_max;
    % upper_bound = sum(upper_bound,1);
    % upper_bound = upper_bound==length(reach_data);
    % accepted = (upper_bound + lower_bound) == 2;
    % out = sum(accepted);

    reach_mean = mean(reach_data,2);
    
    % Fs = 1000;                   % samples per second
    % dt = 1/Fs;                   % seconds per sample
    % StopTime = 3.5;             % seconds
    % t = (0:dt:StopTime-dt)';     % seconds
    % siny = sin(2*pi*60*t)*0.3;
    % siny(1:1500) = 0;
    % reach_mean = reach_mean + siny;

    [sr,fr,tr,psr] = spectrogram(reach_mean,500,480,1000,1000,"power",'yaxis');

    % Mean Normalize
    pretrigger_mean_r = mean(psr(:,1:floor(size(psr,2)/3)),2);
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
    fr = fr(4:40);
    psr = psr(4:40,:);
    imagesc(tr-1,fr,psr)
    axis xy
    xlabel("Time (s)")
    ylabel("Frequency (Hz)")
    colorbar('southoutside')
    % ylim([4 80])
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
        save_folder = 'C:\Users\czhe0008\Documents\EEG\individual\Bottle_ME\';
    elseif mode == 2
        titlestr = strcat('Bottle Observation, All trials, Ch',channel_name);
        figstr = strcat('Bottle_Observation_Alltrials_Norm_Ch',channel_name);
        save_folder = 'C:\Users\czhe0008\Documents\EEG\individual\Bottle_AO\';
    elseif mode == 3
        titlestr = strcat('Bottle Imagery, All trials, Ch',channel_name);
        figstr = strcat('Bottle_Imagery_Alltrials_Norm_Ch',channel_name);
        save_folder = 'C:\Users\czhe0008\Documents\EEG\individual\Bottle_MI\';
    elseif mode == 4
        titlestr = strcat('USB Motor old, All trials, Ch',channel_name);
        figstr = strcat('USB_Motor_old_Alltrials_Norm_Ch',channel_name);
        save_folder = 'C:\Users\czhe0008\Documents\EEG\old_pipe\19_08_test_fig\';
    elseif mode == 5
        titlestr = strcat('USB Observation, All trials, Ch',channel_name);
        figstr = strcat('USB_Observation_Alltrials_Norm_Ch',channel_name);
        save_folder = 'C:\Users\czhe0008\Documents\EEG\individual\USB_AO\';
    elseif mode == 6
        titlestr = strcat('USB Imagery, All trials, Ch',channel_name);
        figstr = strcat('USB_Imagery_Alltrials_Norm_Ch',channel_name);
        save_folder = 'C:\Users\czhe0008\Documents\EEG\individual\USB_MI\';
    end
    sgtitle(titlestr)
    if not(isfolder(save_folder))
        mkdir(save_folder)
    end
    saveas(gcf,strcat(save_folder,figstr,'.png'))
    close(2)
end