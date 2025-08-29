% close all; clear all; clc;
src = "C:\Users\czhe0008\Documents\EEG\raw_data\reach\";
maintemp = dir(fullfile(src,'*'));
mainfolder = setdiff({maintemp([maintemp.isdir]).name},{'.','..'});
for new_ref = 1:31
    save_folder = strcat('C:\Users\czhe0008\Documents\EEG\LDA\rereference\ref',int2str(new_ref));
    if not(isfolder(save_folder))
        mkdir(save_folder)
    end
    for k = 1:31
        if k ~= new_ref
            continue
        end
        grasp_data = [];
        for i = 1:numel(mainfolder)
            temp = cell2mat(mainfolder(i));
            % Determine which folders to grab data from
            if (temp(1:7) == "July_02") && (temp((end-5):end) == "tle_ME") % temp(7:end) == "Right_Grasp"
                subtemp = dir(fullfile(src,mainfolder{i},'*.mat'));
                subfolder = {subtemp(~[subtemp.isdir]).name};
                for j = 1:numel(subfolder)
                    T = fullfile(src,mainfolder{i},subfolder{j});
                    load(T)
                    % reref = save_flag(1501:end,k) - save_flag(1501:end,new_ref);
                    % [s,f,t,ps] = spectrogram(reref,hanning(250),225,1000,1000,"power",'yaxis');
                    % trigger = ceil(size(ps,2)/3);
                    % pretrigger_mean = mean(ps(:,1:trigger),2);
                    % ps = ps./pretrigger_mean;
                    % 
                    % temp1 = [];
                    % for m = trigger:size(ps,2)
                    %     alpha = 8:12;
                    %     beta = 13:30;
                    %     gamma = 31:80;
                    %     temp1 = [temp1,mean(ps(alpha,m)),mean(ps(beta,m)),mean(ps(gamma,m))];
                    % end

                    reref = save_flag(1001:end,k) - save_flag(1001:end,new_ref);
                    [s,f,t,ps] = spectrogram(reref,500,480,1000,1000,"power",'yaxis');
                    pretrigger_mean = mean(ps(:,1:50),2);
                    ps = ps./pretrigger_mean;

                    figure
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

                    temp1 = [];
                    for m = 51:10:150
                        alpha = 8:12;
                        beta = 13:30;
                        gamma = 31:80;
                        temp1 = [temp1,mean(mean(ps(alpha,m:m+9))),mean(mean(ps(beta,m:m+9))),mean(mean(ps(gamma,m:m+9)))];
                    end
                    grasp_data = [grasp_data;temp1];
                end
            end
        end
        dest = strcat(save_folder,'\GMotor',int2str(k),'.mat');
        save(dest, "grasp_data")
    end
end