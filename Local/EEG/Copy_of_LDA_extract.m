close all; clear all; clc;
src = "C:\Users\czhe0008\Documents\EEG\raw_data\reach\";
maintemp = dir(fullfile(src,'*'));
mainfolder = setdiff({maintemp([maintemp.isdir]).name},{'.','..'});
for k = 1:31
    grasp_data = [];
    for i = 1:numel(mainfolder)
        temp = cell2mat(mainfolder(i));
        % Determine which folders to grab data from
        if (temp(1:7) == "July_10") && (temp((end-5):end) == "tle_ME") % temp(7:end) == "Right_Grasp"
            subtemp = dir(fullfile(src,mainfolder{i},'*.mat'));
            subfolder = {subtemp(~[subtemp.isdir]).name};
            for j = 1:numel(subfolder)
                T = fullfile(src,mainfolder{i},subfolder{j});
                load(T)
                [s,f,t,ps] = spectrogram(save_flag(1001:end,k),hanning(500),400,1000,1000,"psd",'yaxis');
                pretrigger_mean = mean(ps(:,1:floor(size(ps,2)/3)),2);
                ps = ps./pretrigger_mean;
                
                [peak_amplitude, peak_location] = max(ps(4:100,11:end));
                temp1 = [peak_amplitude',f(peak_location),t(11:end)'];
                grasp_data = [grasp_data;temp1];
            end
        end
    end
    % grasp_data = grasp_data(1:609,:);
    dest = strcat('C:\Users\czhe0008\Documents\EEG\LDA\spectral\G2Motor',int2str(k),'.mat');
    save(dest, "grasp_data")
end
