close all; clear all; clc;
all_grasp_guess = [];
all_grasp_true = [];
for new_ref = 1:31
    for j = 1:31
        if j == new_ref
            continue
        end
        Gstr = strcat('LDA\rereference\ref',int2str(new_ref),'\GMotor',int2str(j),'.mat');
        Pstr = strcat('LDA\rereference\ref',int2str(new_ref),'\PMotor',int2str(j),'.mat');
        load(Gstr)
        load(Pstr)
        grasp_data = audio_data;
    
        grasp_guess = [];
        grasp_true = [];
        for i = 1:size(grasp_data,1)
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
        if j == 1
            channel = j;
        else
            channel = j + 1;
        end
        % figure(1)
        % confusionchart(grasp_true,grasp_guess);
        % figstr = strcat('Reref',int2str(new_ref),'_GP_Channel',int2str(channel));
        % titlestr = strcat("Grasp vs Pinch reref",int2str(new_ref)," Channel",int2str(channel));
        % title(titlestr)
        % fontsize(18,"points")
        % save_folder = 'C:\Users\czhe0008\Documents\EEG\LDA\rereference_fig\';
        % saveas(gcf,strcat(save_folder,figstr,'.png'))
        % close(1)
    
        all_grasp_guess = [all_grasp_guess;grasp_guess];
        all_grasp_true = [all_grasp_true;grasp_true];
    end
    confusionchart(all_grasp_true,all_grasp_guess);%,'Normalization','row-normalized');
    figstr = strcat('Abs_Reref',int2str(new_ref),'_GP_ChannelAll');
    titlestr = strcat("Grasp vs Pinch reref",int2str(new_ref)," All Channels");
    title(titlestr)
    fontsize(18,"points")
    save_folder = 'C:\Users\czhe0008\Documents\EEG\LDA\rereference_fig\';
    saveas(gcf,strcat(save_folder,figstr,'.png'))
end