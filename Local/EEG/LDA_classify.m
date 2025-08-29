close all; clear all; clc;
% for i = 1:31
%     graspfile = strcat('LDA\test\GMotor',int2str(i),'.mat');
%     pinchfile = strcat('LDA\test\PMotor',int2str(i),'.mat');
%     load(graspfile)
%     load(pinchfile)
%     labels = [repmat("Grasp",size(grasp_data(1:25,:),1),1);repmat("Pinch",size(pinch_data(1:25,:),1),1)];
%     motor = [grasp_data(1:25,:);pinch_data(1:25,:)];
%     MdlLinear = fitcdiscr(motor,labels);
% 
%     % Predicting data outcomes
%     true_out = [repmat("Grasp",size(grasp_data(26:end,:),1),1);repmat("Pinch",size(pinch_data(26:end,:),1),1)];
%     testmotor = [grasp_data(26:end,:);pinch_data(26:end,:)];
%     guess_out = string(predict(MdlLinear,testmotor));
% 
%     figure(1)
%     confusionchart(true_out,guess_out)
%     if i == 1
%         figstr = strcat('LDAJuly',int2str(1));
%         titlestr = strcat('LDA July 02 ch',int2str(1));
%     else
%         figstr = strcat('LDAJuly',int2str(i+1));
%         titlestr = strcat('LDA July 02 ch',int2str(i+1));
%     end
%     title(titlestr)
%     fontsize(18,"points")
%     save_folder = 'C:\Users\czhe0008\Documents\EEG\LDA\15_07_update\July_new_format\';
%     if not(isfolder(save_folder))
%         mkdir(save_folder)
%     end
%     saveas(gcf,strcat(save_folder,figstr,'.png'))
%     close(1)
% end

i=8;
graspfile = strcat('LDA\test3\GMotor',int2str(i),'.mat');
pinchfile = strcat('LDA\test3\PMotor',int2str(i),'.mat');
load(graspfile)
load(pinchfile)
% grasp_data = grasp_data(:,[1,3]);
% pinch_data = pinch_data(:,[1,3]);
labels = [repmat("Grasp",size(grasp_data,1),1);repmat("Pinch",size(pinch_data,1),1)];
motor = [grasp_data;pinch_data];

rng('default') % For reproducibility
cv = cvpartition(labels,'HoldOut',0.10);
trainInds = training(cv);
sampleInds = test(cv);
trainingData = motor(trainInds,:);
sampleData = motor(sampleInds,:);

class = classify(sampleData,trainingData,labels(trainInds));
cm = confusionchart(labels(sampleInds),class);