close all; clear all; clc;
src = "C:\Users\czhe0008\Documents\EEG\raw_data\test\";
maintemp = dir(fullfile(src,'*'));
mainfolder = setdiff({maintemp([maintemp.isdir]).name},{'.','..'});
for k = 1:31
    reach_data = [];
    for i = 1:numel(mainfolder)
        temp = cell2mat(mainfolder(i));
        % Determine which folders to grab data from
        if (temp == "Right_Pinch_Release")%(temp(1:4) == "June") && (temp((end-5):end) == "USB_AO")
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
    [sr,fr,tr,psr] = spectrogram(reach_mean,500,480,1000,1000,"power",'yaxis');

    % Mean Normalize
    pretrigger_mean_r = mean(psr(:,1:floor(width(psr)/2)),2);
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
    titlestr = strcat('Right Pinch Release, All trials, Ch',channel_name);
    figstr = strcat('Right_Pinch_Release_Alltrials_Norm_Ch',channel_name);
    sgtitle(titlestr)
    save_folder = 'C:\Users\czhe0008\Documents\EEG\sessions\30_06\fig\Right_Pinch_Release\';
    if not(isfolder(save_folder))
        mkdir(save_folder)
    end
    saveas(gcf,strcat(save_folder,figstr,'.png'))
    close(2)
end