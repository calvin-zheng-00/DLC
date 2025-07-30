close all; clear all; clc;
src = "C:\Users\czhe0008\Documents\EEG\raw_data\reach\";
maintemp = dir(fullfile(src,'*'));
mainfolder = setdiff({maintemp([maintemp.isdir]).name},{'.','..'});
for k = 1:31
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
                [s,f,t,ps] = spectrogram(save_flag(1001:end,k),500,480,1000,1000,"power",'yaxis');
                pretrigger_mean = mean(ps(:,1:50),2);
                ps = ps./pretrigger_mean;
                temp1 = [];
                for m = 51:10:150
                    % for n = 5:15:80
                    %     temp1 = [temp1,mean(mean(ps(n:n+9,m:m+9)))];
                    % end
                    alpha = 8:12;
                    beta = 13:30;
                    gamma = 31:80;
                    temp1 = [temp1,mean(mean(ps(alpha,m:m+9))),mean(mean(ps(beta,m:m+9))),mean(mean(ps(gamma,m:m+9)))];
                end
                grasp_data = [grasp_data;temp1];
            end
        end
    end
    dest = strcat('C:\Users\czhe0008\Documents\EEG\LDA\move\GMotor',int2str(k),'.mat');
    save(dest, "grasp_data")
end
