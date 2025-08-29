%close all; clear all; clc;

function rawdata_extract(save_folder,filename,ALLEEG,triggers,start_bound,end_bound)
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

    data = ALLEEG.data'; % This contains all of the raw EEG data.
    loops = size(ALLEEG.event);
    loops = loops(2);
    events = [];
    sync = [];
    
    % This for loop finds all the triggers and their associated time stamps
    for i = 2:loops
        temp = ALLEEG.event(i).type;
        temp(1:2) = [];
        temp = str2double(temp);
        events = [events,temp];
        temp2 = ALLEEG.event(i).latency;
        sync = [sync,temp2];
    end
    
    % This section selects only the triggers chosen by the inputs
    start_trigger = find(events == 1);
    reach_trigger = start_trigger(triggers);
    reach_start = sync(reach_trigger);
    
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

function LDA_extract(src,session,start_bound,end_bound,out_loc)
% This function converts the data prepared by the rawdata_extract function
% into a feature matrix to be used in an LDA
% src (input): the location of the EEG data
% session (input): The name of the session to grab data from
% start_bound (input): The lower bound time index within the EEG data
% matrix that the function will grab data from.
% end_bound (input): The upper bound time index within the EEG data
% matrix that the function will grab data from.
% out_loc (input): The output location to save the LDA features to.
% File output format: There will be multiple output files, each being a 2D
% matrix containing the features between the start bound and end bound
% times across all trials, with each channel having its own file.

    maintemp = dir(fullfile(src,'*'));
    mainfolder = setdiff({maintemp([maintemp.isdir]).name},{'.','..'});
    for k = 1:31
        data = [];
        for i = 1:numel(mainfolder)
            temp = cell2mat(mainfolder(i));
            % Determine which session folders to grab data from
            if (temp == session)
                subtemp = dir(fullfile(src,mainfolder{i},'*.mat'));
                subfolder = {subtemp(~[subtemp.isdir]).name};
                for j = 1:numel(subfolder)
                    T = fullfile(src,mainfolder{i},subfolder{j});
                    load(T)
                    % Create spectrogram
                    [~,~,~,ps] = spectrogram(save_flag(start_bound:end_bound,k),500,480,1000,1000,"power",'yaxis');
                    % Normalize spectrogram
                    pretrigger_mean = mean(ps(:,1:50),2);
                    ps = ps./pretrigger_mean;
                    % Select features from spectrogram
                    temp1 = [];
                    for m = 51:10:150
                        alpha = 8:12; % Alpha band
                        beta = 13:30; % Beta band
                        gamma = 31:80; % Gamma band
                        temp1 = [temp1,mean(mean(ps(alpha,m:m+9))),mean(mean(ps(beta,m:m+9))),mean(mean(ps(gamma,m:m+9)))];
                    end
                    data = [data;temp1];
                end
            end
        end
        dest = strcat(out_loc,int2str(k),'.mat');
        save(dest, "data")
    end
end


function [grasp_guess,grasp_true] = LDA(grasp_data,pinch_data)
% This function creates a linear discriminator to classify two classes.
% grasp_data (input): An array of features from a single session of data
% defining a single class
% pinch_data (input): An array of features from a single session of data
% defining a single class
% grasp_guess (output): The total cross validated estimations on the
% classes of the input data
% grasp_true (output): The true class of all the cross validated input data

    grasp_guess = [];
    grasp_true = [];
    for i = 1:size(grasp_data,1)
        % Splitting training and testing sets.
        testdata = i;
        traindata = 1:size(grasp_data,1);
        traindata(testdata) = [];
        % Creating classifier
        
        labels = [repmat("Grasp",size(grasp_data(traindata,:),1),1);repmat("Pinch",size(pinch_data(traindata,:),1),1)];
        motor = [grasp_data(traindata,:);pinch_data(traindata,:)];
        MdlLinear = fitcdiscr(motor,labels);%,'DiscrimType','pseudolinear');
    
        % Predicting data outcomes
        true_out = [repmat("Grasp",size(grasp_data(testdata,:),1),1);repmat("Pinch",size(pinch_data(testdata,:),1),1)];
        testmotor = [grasp_data(testdata,:);pinch_data(testdata,:)];
        guess_out = string(predict(MdlLinear,testmotor));
        
        grasp_guess = [grasp_guess;guess_out];
        grasp_true = [grasp_true;true_out];
    end
end



%% Main code
%% EEG data extract code
rawdata_extract("C:\Users\czhe0008\Documents\EEG\raw_data\test2\Aug_19_pinch_40\","Aug_19_pinch_40",ALLEEG,31:60,2499,4000)

%% LDA feature selection code
LDA_extract("C:\Users\czhe0008\Documents\EEG\raw_data\test2\","Aug_19_pinch_40",1001,4500,"C:\Users\czhe0008\Documents\EEG\LDA\sync_test1\PMotor")

%% LDA code
all_grasp_guess = [];
all_grasp_true = [];
for j = 1:31
    grasp_guess = [];
    grasp_true = [];
    Gstr = strcat('LDA\sync_test1\GMotor',int2str(j),'.mat');
    Pstr = strcat('LDA\sync_test1\PMotor',int2str(j),'.mat');
    load(Gstr)
    grasp_data = data;
    load(Pstr)
    pinch_data = data;

    [grasp_guess,grasp_true] = LDA(grasp_data,pinch_data);

    if j == 1
        channel = j;
    else
        channel = j + 1;
    end
    figure(3)
    confusionchart(grasp_true,grasp_guess);
    figstr = strcat('Channel',int2str(channel));
    titlestr = strcat('Channel',int2str(channel));
    title(titlestr)
    fontsize(18,"points")
    save_folder = 'C:\Users\czhe0008\Documents\EEG\LDA\sync_test_fig\';
    saveas(gcf,strcat(save_folder,figstr,'.png'))
    close(3)

    all_grasp_guess = [all_grasp_guess;grasp_guess];
    all_grasp_true = [all_grasp_true;grasp_true];
end
confusionchart(all_grasp_true,all_grasp_guess,'Normalization','row-normalized');
figstr = strcat('sync_test_ChannelAll');
titlestr = strcat('sync test All Channels');
title(titlestr)
fontsize(18,"points")
save_folder = 'C:\Users\czhe0008\Documents\EEG\LDA\sync_test_fig\';
saveas(gcf,strcat(save_folder,figstr,'.png'))