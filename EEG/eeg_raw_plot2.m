% close all; clear all; clc;
mode = 1;
src = "C:\Users\czhe0008\Documents\EEG\raw_data\test2\";
maintemp = dir(fullfile(src,'*'));
mainfolder = setdiff({maintemp([maintemp.isdir]).name},{'.','..'});
reach_data = [];
for i = 1:numel(mainfolder)
    temp = cell2mat(mainfolder(i));
    % Determine which folders to grab data from
    % if ((temp(1:7) == "June_24")||(temp(1:7) == "June_12")||(temp(1:7) == "July_02")) && (temp((end-5):end) == "USB_MI")
    if (temp == "per_channel")%(temp(1:2) == "Ju") && (temp((end-5):end) == "tle_ME") % temp(7:end) == "Right_Grasp"
        subtemp = dir(fullfile(src,mainfolder{i},'*.mat'));
        subfolder = {subtemp(~[subtemp.isdir]).name};
        for j = 1:numel(subfolder)
            T = fullfile(src,mainfolder{i},subfolder{j});
            load(T)
            reach_mean = mean(save_flag,2);
            [sr,fr,tr,psr] = spectrogram(reach_mean,500,480,1000,1000,"power",'yaxis');
            
            % Mean Normalize
            pretrigger_mean_r = mean(psr(:,1:42),2);
            psr = psr./pretrigger_mean_r;
            
            figure(2)
            subplot(2,1,1)
            hold on
            reach_mean_norm = reach_mean - mean(reach_mean);
            plot(linspace(-1,2,length(reach_mean_norm)),reach_mean_norm)
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
            
            if j == 1
                channel_name = "1";
            else
                channel_name = int2str(j+1);
            end
            titlestr = strcat('Channel data norm, All trials, Ch',channel_name);
            figstr = strcat('Channel_Alltrials_Ch',channel_name);
            save_folder = 'C:\Users\czhe0008\Documents\EEG\process2_test\';
            sgtitle(titlestr)
            if not(isfolder(save_folder))
                mkdir(save_folder)
            end
            saveas(gcf,strcat(save_folder,figstr,'.png'))
            close(2)
        end
    end
end