% close all; clear all; clc;

data_set = 1;
data = ALLEEG(data_set).data';
loops = size(ALLEEG(data_set).event);
loops = loops(2);
events = [];
sync = [];

for i = 1:loops
    temp = ALLEEG(data_set).event(i).type;
    temp(1:2) = [];
    temp = str2double(temp);
    events = [events,temp];
    temp2 = ALLEEG(data_set).event(i).latency;
    sync = [sync,temp2];
end
grasp_trigger = find(events == 3);
grasp_start = sync(grasp_trigger);
pinch_trigger = find(events == 4);
pinch_start = sync(pinch_trigger);

grasp_data = [];
for j = 1:31
    grasp_data_temp = [];
    for i = 1:length(grasp_start)
        data_points = (grasp_start(i)-499):(grasp_start(i)+2000);
        % if (grasp_start(i)+2000)>size(data,1)
        %     continue
        % end
        grasp_data_temp = [grasp_data_temp,data(data_points,j)];
    end
    grasp_data = [grasp_data,mean(grasp_data_temp,2)];
end
pinch_data = [];
for j = 1:31
    pinch_data_temp = [];
    for i = 1:length(grasp_start)
        data_points = (grasp_start(i)-499):(grasp_start(i)+2000);
        % if (grasp_start(i)+2000)>size(data,1)
        %     continue
        % end
        pinch_data_temp = [pinch_data_temp,data(data_points,j)];
    end
    pinch_data = [pinch_data,mean(pinch_data_temp,2)];
end

for i = 1:31
    [s,f,t,ps] = spectrogram(grasp_data(:,i),500,480,1000,1000,"psd",'yaxis');
    pretrigger_mean = mean(ps(:,1:floor(size(ps,2)/5)),2);
    psn = ps./pretrigger_mean;

    figure(4)
    subplot(3,1,1)
    hold on
    reach_mean_norm = grasp_data(:,i) - mean(grasp_data(:,i));
    plot(linspace(-0.5,2,length(reach_mean_norm)),reach_mean_norm)
    xline(0,'-',{'Reach'});
    xlabel("Time (s)")
    ylabel("Amplitude (\muV)")
    title('Time Domain')

    subplot(3,1,2)
    f = f(4:40);
    ps = ps(4:40,:);
    imagesc(t-0.5,f,ps)
    axis xy
    xlabel("Time (s)")
    ylabel("Frequency (Hz)")
    colorbar('southoutside')
    set(gca, ColorScale="log")
    title('Unnormalized Spectrogram')
    % clim([0.0001 1]);

    set(gcf, 'Position', [100 0 1200 1000])
    subplot(3,1,3)
    psn = psn(4:40,:);
    imagesc(t-0.5,f,psn)
    axis xy
    xlabel("Time (s)")
    ylabel("Frequency (Hz)")
    colorbar('southoutside')
    set(gca, ColorScale="log")
    title('Normalized Spectrogram')
    % clim([0.002 20]);
    colormap(jet)

    if i == 1
        channel_name = "1";
    else
        channel_name = int2str(i+1);
    end
    titlestr = strcat('New baseline 24 06 grasp, Ch',channel_name);
    figstr = strcat('New_baseline_24_06_grasp_Ch',channel_name);
    save_folder = 'C:\Users\czhe0008\Documents\EEG\comparison\';
    sgtitle(titlestr)
    if not(isfolder(save_folder))
        mkdir(save_folder)
    end
    saveas(gcf,strcat(save_folder,figstr,'.png'))
    close(4)

    % temp1 = [];
    % for m = 42:12:126
    %     alpha = 8:12;
    %     beta = 13:30;
    %     gamma = 31:40;
    %     temp1 = [temp1,mean(mean(ps(alpha,m:m+11))),mean(mean(ps(beta,m:m+11))),mean(mean(ps(gamma,m:m+11)))];
    % end
    % dest = strcat('C:\Users\czhe0008\Documents\EEG\LDA\new_window\PMotor',int2str(k),'.mat');
    % save(dest, "pinch_std_data")
end