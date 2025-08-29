close all; clear all; clc;
% Creating classifier
load GMotor.mat
load PMotor.mat
labels = [repmat("Grasp",size(grasp_data,1),1);repmat("Pinch",size(grasp_data,1),1)];
% pinch_edit = pinch_data(:,3)>30;
% pinch_data(pinch_edit,1)=pinch_data(pinch_edit,1)+2;
motor = [grasp_data;pinch_data];
MdlLinear = fitcdiscr(motor,labels,'DiscrimType','pseudolinear');

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

    if size(grasp_data,2)>size(pinch_data,2)
        grasp_data = grasp_data(:,1:size(pinch_data,2));
    elseif size(grasp_data,2)<size(pinch_data,2)
        pinch_data = pinch_data(:,1:size(grasp_data,2));
    end

    % Making spectrograms
    for l = 1:size(grasp_data,2)
        [sg,fg,tg,psg] = spectrogram(grasp_data(:,l),500,480,1000,1000,"psd",'yaxis');
        [sp,fp,tp,psp] = spectrogram(pinch_data(:,l),500,480,1000,1000,"psd",'yaxis');

        % Mean Normalize
        pretrigger_mean_g = mean(psg(:,1:50),2);
        psg = psg./pretrigger_mean_g;
        pretrigger_mean_p = mean(psp(:,1:50),2);
        psp = psp./pretrigger_mean_p;
    
        % Reforming data into classes
        temp1 = [];
        for i = 50:151
            for j = 4:100
                temp2 = [psg(j,i),tg(i),fg(j),k,l];
                temp1 = [temp1;temp2];
            end
        end
        all_grasp_data = [all_grasp_data;temp1];
        temp1 = [];
        for i = 50:151
            for j = 4:100
                temp2 = [psg(j,i),tg(i),fg(j),k,l];
                temp1 = [temp1;temp2];
            end
        end
        all_pinch_data = [all_pinch_data;temp1];
    end
end
% Predicting data outcomes
grasp_result = [];
pinch_result = [];
for l = 1:30
    test_grasp = all_grasp_data((all_grasp_data(:,5)==l),1:4);
    test_pinch = all_pinch_data((all_pinch_data(:,5)==l),1:4);
    % test_pinch2 = test_pinch;
    % test_pinch_edit = test_pinch2(:,3)>30;
    % test_pinch2(test_pinch_edit,1)=test_pinch2(test_pinch_edit,1)+2;
    
    guess_grasp = predict(MdlLinear,test_grasp);
    temp1 = sum(guess_grasp == "Grasp");
    temp2 = length(guess_grasp)/2;
    if sum(guess_grasp == "Grasp") > length(guess_grasp)/2
        grasp_result = [grasp_result;"Grasp"];
    else
        grasp_result = [grasp_result;"Pinch"];
    end
    guess_pinch = predict(MdlLinear,test_pinch);
    if sum(guess_pinch == "Pinch") > length(guess_pinch)/2
        pinch_result = [pinch_result;"Pinch"];
    else
        pinch_result = [pinch_result;"Grasp"];
    end
    % guess_pinch2 = predict(MdlLinear,test_pinch2);
    % true_pinch2 = cellstr(repmat("Pinch",length(guess_pinch2),1));
    % figure
    % confusionchart([true_grasp;true_pinch],[guess_grasp;guess_pinch])
    % title('LDA with unedited training pinch, unedited test pinch')
    % figure
    % confusionchart([true_grasp;true_pinch2],[guess_grasp;guess_pinch2])
    % title('LDA with unedited training pinch, edited test pinch')
end
true_grasp = cellstr(repmat("Grasp",length(grasp_result),1));
true_pinch = cellstr(repmat("Pinch",length(pinch_result),1));
figure
confusionchart([true_grasp;true_pinch],[cellstr(grasp_result);cellstr(pinch_result)])
title('30 True Grasp and 30 True Pinch LDA estimates')