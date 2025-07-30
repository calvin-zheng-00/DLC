close all; clear all; clc;
all_grasp_guess = [];
all_grasp_true = [];
for j = 1:31
    grasp_guess = [];
    grasp_true = [];
    Gstr = strcat('LDA\24_07_update\PPJun25\GMotor',int2str(j),'.mat');
    Pstr = strcat('LDA\24_07_update\PPJun25\PMotor',int2str(j),'.mat');
    load(Gstr)
    load(Pstr)
    % grasp_data = pre_data;
    % pinch_data = post_data;

    % motor = [grasp_data;pinch_data];
    % [coeff,score,latent,tsquared,explained,mu] = pca(motor, 'Rows', 'pairwise'); %pAngles
    % grasp_data = score(1:30,1:15);
    % pinch_data = score(31:60,1:15);
    
    for i = 1:30
        testdata = i;
        traindata = 1:30;
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
    if j == 1
        channel = j;
    else
        channel = j + 1;
    end
    % figure(1)
    % confusionchart(grasp_true,grasp_guess);
    % figstr = strcat('Channel',int2str(channel));
    % titlestr = strcat('Channel',int2str(channel));
    % title(titlestr)
    % fontsize(18,"points")
    % save_folder = 'C:\Users\czhe0008\Documents\EEG\LDA\GPJul2_fig\';
    % saveas(gcf,strcat(save_folder,figstr,'.png'))
    % close(1)

    all_grasp_guess = [all_grasp_guess;grasp_guess];
    all_grasp_true = [all_grasp_true;grasp_true];
end
confusionchart(all_grasp_true,all_grasp_guess,'Normalization','row-normalized');
figstr = strcat('GPJul2_ChannelAll');
titlestr = strcat('Jul 2 Grasp vs Pinch All Channels');
title(titlestr)
fontsize(18,"points")
save_folder = 'C:\Users\czhe0008\Documents\EEG\LDA\percentage_fig\';
saveas(gcf,strcat(save_folder,figstr,'.png'))