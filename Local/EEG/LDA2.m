% close all; clear all; clc;
all_grasp_guess = [];
all_grasp_true = [];
for j = 1:31
    grasp_guess = [];
    grasp_true = [];
    Gstr = strcat('LDA\new_window\GMotor',int2str(j),'.mat');
    Pstr = strcat('LDA\new_window\PMotor',int2str(j),'.mat');
    load(Gstr)
    load(Pstr)
    % grasp_data = grasp_std_data;
    % pinch_data = pinch_std_data;

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
        MdlLinear = fitcdiscr(motor,labels,'DiscrimType', 'linear', 'Gamma', 0.04, 'Delta', 0.8676);
        % [err, gamma_optimal, delta_optimal, numpred] = cvshrink(MdlLinear, ...
        % 'NumGamma', 25, 'NumDelta', 25, 'Verbose', 1);
        % plot(err,numpred,'k.')
        % xlabel('Error rate')
        % ylabel('Number of predictors')
    
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
    figure(1)
    subplot(1,2,1)
    confusionchart(grasp_true,grasp_guess,'Normalization','row-normalized');
    figstr = strcat('GP_2_Channel',int2str(channel));
    titlestr = strcat('Grasp vs Pinch 2 ch',int2str(channel));
    title(titlestr)
    subplot(1,2,2)
    confusionchart(grasp_true,grasp_guess);
    titlestr = strcat('Grasp vs Pinch 2 ch',int2str(channel));
    title(titlestr)
    fontsize(18,"points")
    set(gcf, 'Position', [100 100 800 350])
    save_folder = 'C:\Users\czhe0008\Documents\EEG\LDA\new_window_fig\';
    saveas(gcf,strcat(save_folder,figstr,'.png'))
    close(1)

    all_grasp_guess = [all_grasp_guess;grasp_guess];
    all_grasp_true = [all_grasp_true;grasp_true];
end
subplot(1,2,1)
confusionchart(all_grasp_true,all_grasp_guess,'Normalization','row-normalized');
figstr = strcat('GP_2_ChannelAll');
titlestr = strcat('Grasp vs Pinch 2 All Ch');
title(titlestr)
subplot(1,2,2)
confusionchart(all_grasp_true,all_grasp_guess);
titlestr = strcat('Grasp vs Pinch 2 All Ch');
title(titlestr)
fontsize(18,"points")
set(gcf, 'Position', [100 100 800 350])
save_folder = 'C:\Users\czhe0008\Documents\EEG\LDA\new_window_fig\';
saveas(gcf,strcat(save_folder,figstr,'.png'))