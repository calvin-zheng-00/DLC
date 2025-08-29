close all; clear all; clc;
% Creating classifier
load GMotor.mat
load PMotor.mat
labels = [repmat("Grasp",size(grasp_data,1),1);repmat("Pinch",size(grasp_data,1),1)];
motor = [grasp_data;pinch_data];
MdlLinear = fitcdiscr(motor,labels,'DiscrimType','pseudolinear');
pinch_edit = pinch_data(:,3)>30;
pinch_data(pinch_edit,1)=pinch_data(pinch_edit,1)+2;
motor = [grasp_data;pinch_data];
MdlLinear2 = fitcdiscr(motor,labels,'DiscrimType','pseudolinear');
grasp_edit = grasp_data(:,3)>30;
grasp_data(grasp_edit,1)=grasp_data(grasp_edit,1)+2;
motor = [grasp_data;pinch_data];
MdlLinear3 = fitcdiscr(motor,labels,'DiscrimType','pseudolinear');

% Getting trial data. One set of grasp data and one set of pinch data
all_grasp_data = [];
all_pinch_data = [];
src_grasp = "data\grasp";
src_pinch = "data\pinch";
grasptemp = dir(fullfile(src_grasp,'*.mat'));
graspfolder = {grasptemp(~[grasptemp.isdir]).name};
pinchtemp = dir(fullfile(src_pinch,'*.mat'));
pinchfolder = {pinchtemp(~[pinchtemp.isdir]).name};
for k = 1:31
    % Gathering grasp data
    grasp_data = [];
    for i = 1:numel(graspfolder)
        T = fullfile(src_grasp,graspfolder{i});
        load(T)
        grasp_data = [grasp_data,mean(save_flag(1:end,k),2)];
    end
    grasp_data = grasp_data(1001:end,:);
    % Bad trial removal
    amp = max(grasp_data)-min(grasp_data);
    amp_std = std(amp);
    amp_mean = mean(amp);
    amp_min = amp_mean - 4*amp_std;
    amp_max = amp_mean + 4*amp_std;
    amp_thresh = (amp<amp_max)&(amp>amp_min);
    grasp_data = grasp_data(:,amp_thresh);

    % Gathering pinch data
    pinch_data = [];
    for i = 1:numel(pinchfolder)
        T = fullfile(src_pinch,pinchfolder{i});
        load(T)
        pinch_data = [pinch_data,mean(save_flag(1:end,k),2)];
    end
    pinch_data = pinch_data(1001:end,:);
    % Bad trial removal
    amp = max(pinch_data)-min(pinch_data);
    amp_std = std(amp);
    amp_mean = mean(amp);
    amp_min = amp_mean - 4*amp_std;
    amp_max = amp_mean + 4*amp_std;
    amp_thresh = (amp<amp_max)&(amp>amp_min);
    pinch_data = pinch_data(:,amp_thresh);

    % Making spectrograms
    grasp_mean = mean(grasp_data,2);
    [sg,fg,tg,psg] = spectrogram(grasp_mean,500,480,1000,1000,"psd",'yaxis');
    pinch_mean = mean(pinch_data,2);
    [sp,fp,tp,psp] = spectrogram(pinch_mean,500,480,1000,1000,"psd",'yaxis');

    % Mean Normalize
    pretrigger_mean_g = mean(psg(:,1:50),2);
    psg = psg./pretrigger_mean_g;
    pretrigger_mean_p = mean(psp(:,1:50),2);
    psp = psp./pretrigger_mean_p;

    % Reforming data into classes
    temp1 = [];
    for i = 50:151
        for j = 4:100
            temp2 = [psg(j,i),tg(i),fg(j),k];
            temp1 = [temp1;temp2];
        end
    end
    all_grasp_data = [all_grasp_data;temp1];
    temp1 = [];
    for i = 50:151
        for j = 4:100
            temp2 = [psp(j,i),tp(i),fp(j),k];
            temp1 = [temp1;temp2];
        end
    end
    all_pinch_data = [all_pinch_data;temp1];
end
% Predicting data outcomes
all_grasp_data2 = all_grasp_data;
all_grasp_edit = all_grasp_data2(:,3)>30;
all_grasp_data2(all_grasp_edit,1)=all_grasp_data2(all_grasp_edit,1)+2;
all_pinch_data2 = all_pinch_data;
all_pinch_edit = all_pinch_data2(:,3)>30;
all_pinch_data2(all_pinch_edit,1)=all_pinch_data2(all_pinch_edit,1)+2;

guess_grasp = predict(MdlLinear,all_grasp_data);
true_grasp = cellstr(repmat("Grasp",length(guess_grasp),1));
guess_pinch = predict(MdlLinear,all_pinch_data);
true_pinch = cellstr(repmat("Pinch",length(guess_pinch),1));
guess_grasp2 = predict(MdlLinear,all_grasp_data2);
true_grasp2 = cellstr(repmat("Grasp",length(guess_grasp2),1));
guess_pinch2 = predict(MdlLinear,all_pinch_data2);
true_pinch2 = cellstr(repmat("Pinch",length(guess_pinch2),1));

% figure
% confusionchart([true_grasp;true_pinch],[guess_grasp;guess_pinch])
% title('LDA with unedited training grasp, unedited test grasp')
% 
% figure
% confusionchart([true_grasp2;true_pinch],[guess_grasp2;guess_pinch])
% title('LDA with unedited training grasp, edited test grasp')
% 
% guess_grasp = predict(MdlLinear3,all_grasp_data);
% true_grasp = cellstr(repmat("Grasp",length(guess_grasp),1));
% guess_pinch = predict(MdlLinear3,all_pinch_data);
% true_pinch = cellstr(repmat("Pinch",length(guess_pinch),1));
% guess_grasp2 = predict(MdlLinear3,all_grasp_data2);
% true_grasp2 = cellstr(repmat("Grasp",length(guess_grasp2),1));
% 
% figure
% confusionchart([true_grasp;true_pinch],[guess_grasp;guess_pinch])
% title('LDA with edited training grasp, unedited test grasp')
% 
% figure
% confusionchart([true_grasp2;true_pinch],[guess_grasp2;guess_pinch])
% title('LDA with edited training grasp, edited test grasp')

figure
confusionchart(true_pinch,guess_pinch)
title('LDA with unedited training pinch, unedited test pinch')

figure
confusionchart(true_pinch2,guess_pinch2)
title('LDA with unedited training pinch, edited test pinch')

guess_pinch = predict(MdlLinear2,all_pinch_data);
true_pinch = cellstr(repmat("Pinch",length(guess_pinch),1));
guess_pinch2 = predict(MdlLinear2,all_pinch_data2);
true_pinch2 = cellstr(repmat("Pinch",length(guess_pinch2),1));

figure
confusionchart(true_pinch,guess_pinch)
title('LDA with edited training pinch, unedited test pinch')

figure
confusionchart(true_pinch2,guess_pinch2)
title('LDA with edited training pinch, edited test pinch')