%close all; clear all; clc;

function rawdata_extract(save_folder,filename,ALLEEG,data_set,trigger,start_bound,end_bound)
% This function extracts the data opened by EEGLAB and saves it into a file
% to be used by other functions.
% filename (input): the name of the output file
% save_folder (input): The location the output file will be saved to
% triggers (input): An array identifying the triggers that represent the
% data to be saved into the output file.
% start_bound (input): The amount of time (ms) to include before the
% trigger
% end_bound (input): The amount of time (ms) to include after the trigger
% File output format: There will be multiple output files, each being a 2D
% matrix containing the raw EEG data between the start bound and end bound
% times across trials (triggers), with each channel having its own file.

    data = ALLEEG(data_set).data'; % This contains all of the raw EEG data.
    loops = size(ALLEEG(data_set).event);
    loops = loops(2);
    events = [];
    sync = [];
    
    % This for loop finds all the triggers and their associated time stamps
    for i = 2:loops
        temp = ALLEEG(data_set).event(i).type;
        temp(1:2) = [];
        temp = str2double(temp);
        events = [events,temp];
        temp2 = ALLEEG(data_set).event(i).latency;
        sync = [sync,temp2];
    end
    
    % This section selects only the triggers chosen by the inputs
    start_trigger = find(events == trigger);
    reach_start = sync(start_trigger);
    
    % This section formats all the data and saves it.
    for i = 1:length(reach_start)
        data_points = (reach_start(i)-start_bound):(reach_start(i)+end_bound);
        save_flag = data(data_points,:);
        if not(isfolder(save_folder))
            mkdir(save_folder)
        end
        save_loc = strcat(save_folder,filename,int2str(i),'.mat');
        save(save_loc, "save_flag")
    end
end

function plotting_old(src,session,start_bound,end_bound,out_loc,fname)
% This function plots the spectrogram of the data prepared by the
% rawdata_extract function.
% src (input): the location of the EEG data
% session (input): The name of the session to grab data from
% start_bound (input): The lower bound time index within the EEG data
% matrix that the function will grab data from.
% end_bound (input): The upper bound time index within the EEG data
% matrix that the function will grab data from.
% out_loc (input): The output location to save the LDA features to.
% File output format: A figure for each channel within the chosen session.
% Each figure contains the time domain data, an unnormalized spectrogram,
% and a normalized spectrogram.

    maintemp = dir(fullfile(src,'*'));
    mainfolder = setdiff({maintemp([maintemp.isdir]).name},{'.','..'});
    % Looping through channels
    for k = 1:31
        reach_data = [];
        for i = 1:numel(mainfolder)
            temp = cell2mat(mainfolder(i));
            % Determine which session folders to grab data from
            if (temp == session)
                subtemp = dir(fullfile(src,mainfolder{i},'*.mat'));
                subfolder = {subtemp(~[subtemp.isdir]).name};
                for j = 1:numel(subfolder)
                    T = fullfile(src,mainfolder{i},subfolder{j});
                    load(T)
                    reach_data_prep = save_flag(start_bound:end_bound,k);
                    reach_data = [reach_data,reach_data_prep];
                end
            end
        end
        % Average across trials
        reach_mean = mean(reach_data,2);
        % Create spectrogram
        [~,f,t,ps] = spectrogram(reach_mean,500,480,1000,1000,"power",'yaxis');
        % Normalize spectrogram
        pretrigger_mean = mean(ps(:,1:50),2);
        psn = ps./pretrigger_mean;
        % Plotting spectrogram
        figure(4)
        subplot(3,1,1)
        hold on
        reach_mean_norm = reach_mean - mean(reach_mean);
        plot(linspace(-1,2,length(reach_mean_norm(500:end))),reach_mean_norm(500:end))
        xline(0,'-',{'Reach'});
        xlabel("Time (s)")
        ylabel("Amplitude (\muV)")
        title('Time Domain')
    
        subplot(3,1,2)
        f = f(4:40);
        ps = ps(4:40,:);
        imagesc(t-1,f,ps)
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
        imagesc(t-1,f,psn)
        axis xy
        xlabel("Time (s)")
        ylabel("Frequency (Hz)")
        colorbar('southoutside')
        set(gca, ColorScale="log")
        title('Normalized Spectrogram')
        % clim([0.002 20]);
        colormap(jet)
    
        if k == 1
            channel_name = "1";
        else
            channel_name = int2str(k+1);
        end
        fname2 = fname;
        fname2(isspace(fname)) = '_';
        titlestr = strcat(fname,channel_name);
        figstr = strcat(fname2,channel_name);
        sgtitle(titlestr)
        if not(isfolder(out_loc))
            mkdir(out_loc)
        end
        saveas(gcf,strcat(out_loc,figstr,'.png'))
        close(4)
    end
end

function plotting_new(src,session,start_bound,end_bound,out_loc,fname)
% This function plots the spectrogram of the data prepared by the
% rawdata_extract function.
% src (input): the location of the EEG data
% session (input): The name of the session to grab data from
% start_bound (input): The lower bound time index within the EEG data
% matrix that the function will grab data from.
% end_bound (input): The upper bound time index within the EEG data
% matrix that the function will grab data from.
% out_loc (input): The output location to save the LDA features to.
% File output format: A figure for each channel within the chosen session.
% Each figure contains the time domain data, an unnormalized spectrogram,
% and a normalized spectrogram.

    maintemp = dir(fullfile(src,'*'));
    mainfolder = setdiff({maintemp([maintemp.isdir]).name},{'.','..'});
    % Looping through channels
    for k = 1:31
        reach_data = [];
        for i = 1:numel(mainfolder)
            temp = cell2mat(mainfolder(i));
            % Determine which session folders to grab data from
            if (temp == session)
                subtemp = dir(fullfile(src,mainfolder{i},'*.mat'));
                subfolder = {subtemp(~[subtemp.isdir]).name};
                for j = 1:numel(subfolder)
                    T = fullfile(src,mainfolder{i},subfolder{j});
                    load(T)
                    reach_data_prep = save_flag(start_bound:end_bound,k);
                    reach_data = [reach_data,reach_data_prep];
                end
            end
        end
        % Average across trials
        reach_mean = mean(reach_data,2);
        % Create spectrogram
        [~,f,t,ps] = spectrogram(reach_mean,500,480,1000,1000,"power",'yaxis');
        % Normalize spectrogram
        pretrigger_mean = mean(ps(:,1:floor(size(ps,2)/3)),2);
        psn = ps./pretrigger_mean;
        % Plotting spectrogram
        figure(4)
        subplot(3,1,1)
        hold on
        reach_mean_norm = reach_mean - mean(reach_mean);
        plot(linspace(-1,2,length(reach_mean_norm)),reach_mean_norm)
        xline(0,'-',{'Reach'});
        xlabel("Time (s)")
        ylabel("Amplitude (\muV)")
        title('Time Domain')
    
        subplot(3,1,2)
        f = f(4:40);
        ps = ps(4:40,:);
        imagesc(t-1,f,ps)
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
        imagesc(t-1,f,psn)
        axis xy
        xlabel("Time (s)")
        ylabel("Frequency (Hz)")
        colorbar('southoutside')
        set(gca, ColorScale="log")
        title('Normalized Spectrogram')
        % clim([0.002 20]);
        colormap(jet)
    
        if k == 1
            channel_name = "1";
        else
            channel_name = int2str(k+1);
        end
        fname2 = fname;
        fname2(isspace(fname)) = '_';
        titlestr = strcat(fname,channel_name);
        figstr = strcat(fname2,channel_name);
        sgtitle(titlestr)
        if not(isfolder(out_loc))
            mkdir(out_loc)
        end
        saveas(gcf,strcat(out_loc,figstr,'.png'))
        close(4)
    end
end


%% Main code
%% EEG data extract code
rawdata_extract("C:\Users\czhe0008\Documents\EEG\raw_data\comparison\19_08_grasp\","19_08_grasp",ALLEEG,1,3,2499,4000)

%% Plotting EEG data using old method
plotting_old("C:\Users\czhe0008\Documents\EEG\raw_data\comparison\","19_08_grasp",1001,4500,"C:\Users\czhe0008\Documents\EEG\comparison\",'Old 19 08 grasp Ch')

%% Plotting EEG data using new method
plotting_new("C:\Users\czhe0008\Documents\EEG\raw_data\comparison\","New_24_06_grasp",1501,4500,"C:\Users\czhe0008\Documents\EEG\comparison\",'New 24 06 grasp Ch')