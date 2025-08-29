close all; clear all; clc;
src = "C:\Users\czhe0008\Documents\EEG\raw_data\test2\";
maintemp = dir(fullfile(src,'*'));
mainfolder = setdiff({maintemp([maintemp.isdir]).name},{'.','..'});
for k = 1:31
    pinch_std_data = [];
    for i = 1:numel(mainfolder)
        temp = cell2mat(mainfolder(i));
        % Determine which folders to grab data from
        if (temp == "Aug_14_audio_std") %(temp(1:7) == "July_02") && (temp((end-5):end) == "USB_ME") % temp(7:end) == "Right_Grasp" % (temp == "Aug_14_grasp_std")
            subtemp = dir(fullfile(src,mainfolder{i},'*.mat'));
            subfolder = {subtemp(~[subtemp.isdir]).name};
            for j = 1:numel(subfolder)
                T = fullfile(src,mainfolder{i},subfolder{j});
                load(T)
                [s,f,t,ps] = spectrogram(save_flag(1001:(end-2000),k),500,480,1000,1000,"power",'yaxis');
                pretrigger_mean = mean(ps(:,1:50),2);
                ps = ps./pretrigger_mean;
                temp1 = [];
                for m = 51:10:150
                    alpha = 8:12;
                    beta = 13:30;
                    gamma = 31:80;
                    temp1 = [temp1,mean(mean(ps(alpha,m:m+9))),mean(mean(ps(beta,m:m+9))),mean(mean(ps(gamma,m:m+9)))];
                end
                % for m = 1:size(ps,2)
                %     alpha = 8:12;
                %     beta = 13:30;
                %     gamma = 31:80;
                %     temp1 = [temp1,mean(mean(ps(alpha,m))),mean(mean(ps(beta,m))),mean(mean(ps(gamma,m)))];
                % end
                pinch_std_data = [pinch_std_data;temp1];
            end
        end
    end
    dest = strcat('C:\Users\czhe0008\Documents\EEG\LDA\new_window\PMotor',int2str(k),'.mat');
    save(dest, "pinch_std_data")
end
