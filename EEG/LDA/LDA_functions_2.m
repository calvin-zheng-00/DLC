function LDA_extract_adjust(src,session,start_bound,end_bound,out_loc)
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

    % Adjusting trials based on participant reaction
    temp_offset = [13,20,12,17,18,18,14,21,16,12,15,15,13,12,15,16,29,17,17,39,27,20,11,20,15,26,23,24,22,13]*0.025 + 0.5;
    % temp_offset = [12,12,9,18,11,10,23,13,15,17,11,24,10,13,13] * 0.025;

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
                    % filtered = bandstop(save_flag(:,k),[40 60],1000);
                    % filtered = bandstop(filtered,[90 110],1000);
                    % filtered = bandstop(filtered,[140 160],1000);
                    % filtered = bandstop(filtered,[190 210],1000);
                    % filtered = bandstop(filtered,[240 260],1000);
                    % filtered = bandstop(filtered,[290 310],1000);
                    filtered = save_flag(:,k);
                    [max_val,max_idx] = max(filtered(start_bound:end));
                    [min_val,min_idx] = min(filtered(start_bound:end));
                    data = [data;[max_val,max_idx,min_val,min_idx]];
                    % Create spectrogram
                    % [~,~,t,ps] = spectrogram(filtered,500,480,1000,1000,"power",'yaxis');
                    % % Normalize spectrogram
                    % cue = floor(size(ps,2)*(start_bound/(start_bound+end_bound)));
                    % pretrigger_mean = mean(ps(:,1:cue),2);
                    % ps = ps./pretrigger_mean;
                    % % [~, index] = min(abs(t-temp_offset(j)));
                    % [~, index] = min(abs(t-mean(temp_offset))); % This is the offset for audio only
                    % % Select features from spectrogram
                    % alpha = 8:12; % Alpha band
                    % beta = 13:30; % Beta band
                    % gamma = 31:80; % Gamma band
                    % data = [data;[mean(ps(alpha,index)),mean(ps(beta,index)),mean(ps(gamma,index)),max(filtered(start_bound:end))]];

                end
            end
        end
        dest = strcat(out_loc,int2str(k),'.mat');
        save(dest, "data")
    end
end

function LDA_extract(src,session,start_bound,end_bound,out_loc)
    maintemp = dir(fullfile(src,'*'));
    mainfolder = setdiff({maintemp([maintemp.isdir]).name},{'.','..'});
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
                signal = [];
                for k = 1:31%[1,4,14,15,19,20,25,27,28,30,31]
                    filtered = bandstop(save_flag(:,k),[40 60],1000);
                    % filtered = bandstop(filtered,[90 110],1000);
                    % filtered = bandstop(filtered,[140 160],1000);
                    % filtered = bandstop(filtered,[190 210],1000);
                    % filtered = bandstop(filtered,[240 260],1000);
                    % filtered = bandstop(filtered,[290 310],1000);

                    reach_mean_norm = filtered - mean(filtered);
                    plot(linspace(-1,2,length(reach_mean_norm)),reach_mean_norm,'LineWidth',2)
                    
                    % signal = [signal;filtered];
                    [max_val,max_idx] = max(filtered(start_bound:end));
                    [min_val,min_idx] = min(filtered(start_bound:end));
                    signal = [signal,max_val,max_idx,min_val,min_idx];
                end
                data = [data;signal];
                % mean_signal = mean(signal,2);
                % [max_val,max_idx] = max(mean_signal(start_bound:end));
                % [min_val,min_idx] = min(mean_signal(start_bound:end));
                % data = [data;[max_val,max_idx,min_val,min_idx]];
            end
        end
    end
    dest = strcat(out_loc,'.mat');
    save(dest, "data")
end

%Aug_19_grasp_40
%Aug_14_audio_clean
% LDA_extract("C:\Users\czhe0008\Documents\EEG\raw_data\reach\","12_06_grasp_ME",1000,2000,"C:\Users\czhe0008\Documents\EEG\LDA\All\GMotor_12_06_2")
close all; clear all; clc;
for j = 1:1
    Gstr = strcat('LDA\All\GMotor_25_09_1.mat');
    load(Gstr)
    grasp_data = data;
    Pstr = strcat('LDA\All\PMotor_25_09_1.mat');
    load(Pstr)
    pinch_data = data;
    class1 = "Grasp";
    class2 = "Pinch";
    labels = [repmat(class1,size(grasp_data,1),1);repmat(class2,size(pinch_data,1),1)];
    motor = [grasp_data;pinch_data];

    Mdl = fitcsvm(motor, labels);%, 'KernelFunction', 'rbf', 'Standardize', true);
    % Mdl = fitcdiscr(motor,labels);
    cvMdl = crossval(Mdl, 'kfold', 30); % 15-fold cross-validation
    cvError = kfoldLoss(cvMdl);
    [predictedLabels, scores] = kfoldPredict(cvMdl);
    predictedLabels = string(predictedLabels);

    figure(5)
    confusionchart(labels,predictedLabels,'Normalization','row-normalized');
    figstr = strcat("25_09_1_grasp_comp");
    titlestr = strcat("25 09 1 grasp comp");
    title(titlestr)
    fontsize(18,"points")
    save_folder = 'C:\Users\czhe0008\Documents\EEG\LDA\LDA_All_fig\';
    saveas(gcf,strcat(save_folder,figstr,'.png'))
    close(5)
end

for j = 1:1
    % LDA_filtered2\GMotor is for unfiltered data.
    % Gstr = strcat('LDA\LDA_filtered\G2Motor',int2str(j),'.mat');
    Gstr = strcat('LDA\All\GMotor2.mat');
    load(Gstr)
    grasp_data = data;
    % Astr = strcat('LDA\LDA_filtered\A2Motor',int2str(j),'.mat');
    Astr = strcat('LDA\All\PMotor2.mat');
    load(Astr)
    audio_data = data;
    class1 = "Grasp";
    class2 = "Pinch";
    % grasp_data(2,:) = [];
    % [predictedLabels,scores,cvError] = LDA(grasp_data,audio_data,class1,class2);
    labels = [repmat(class1,size(grasp_data,1),1);repmat(class2,size(audio_data,1),1)];
    motor = [grasp_data;audio_data];

    % for i = 1:124
    %     for m = 21:21
    %         if m == i
    %             continue
    %         end
    %         figure(6)
    %         hold on
    %         f1 = motor(:,m);
    %         f2 = motor(:,i);
    %         Mdl = fitcdiscr([f1,f2],labels);
    %         K = Mdl.Coeffs(1,2).Const;  
    %         L = Mdl.Coeffs(1,2).Linear;
    %         f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
    %         h2 = fimplicit(f, [min(f1) max(f1) min(f2) max(f2)]);
    %         h2.Color = 'k';
    %         h2.LineWidth = 2;
    %         h1 = gscatter(f1,f2,labels,'rb','ov',[],'off');
    %         titlestr = strcat("25 09 2 Feature ",int2str(m), " Feature ",int2str(i));
    %         xstr = strcat("Feature ",int2str(m));
    %         ystr = strcat("Feature ",int2str(i));
    %         title(titlestr)
    %         xlabel(xstr)
    %         ylabel(ystr)
    %         save_folder = 'C:\Users\czhe0008\Documents\EEG\LDA\LDA_filtered2_fig\';
    %         saveas(gcf,strcat(save_folder,titlestr,'.png'))
    %         close(6)
    %     end
    % end

    Mdl = fitcdiscr(motor,labels);
    cvMdl = crossval(Mdl, 'kfold', 15); % 15-fold cross-validation
    cvError = kfoldLoss(cvMdl);
    [predictedLabels, scores] = kfoldPredict(cvMdl);
    predictedLabels = string(predictedLabels);

    % [W, Lambda] = eig(Mdl.BetweenSigma, Mdl.Sigma);
    % Lambda = diag(Lambda);
    % [Lambda, SortOrder] = sort(Lambda, 'descend');
    % W = W(:, SortOrder);
    % Y = motor*W;
    % W is eigenvectors
    % Y is the equivalent of the score output in the pca function.
    % Y can visualize in the feature space by e.g. using the first 2-3
    % columns

    % figure(3)
    % plot(1:124,Mdl.Coeffs(1, 2).Linear);
    % titlestr = strcat("25 09 2 Coefficients");
    % xstr = strcat("Coefficients");
    % ystr = strcat("Weight");
    % title(titlestr)
    % xlabel(xstr)
    % ylabel(ystr)
    % xlim([1,124])
    % save_folder = 'C:\Users\czhe0008\Documents\EEG\LDA\LDA_filtered2_fig\';
    % saveas(gcf,strcat(save_folder,titlestr,'.png'))

    figure(5)
    confusionchart(labels,predictedLabels,'Normalization','row-normalized');
    % figstr = strcat("19_08_Avg_Channel_",int2str(channel));
    % titlestr = strcat("19 08 Avg Channel ",int2str(channel));
    figstr = strcat("25_09_2_grasp_comp");
    titlestr = strcat("25 09 2 grasp comp");
    title(titlestr)
    fontsize(18,"points")
    save_folder = 'C:\Users\czhe0008\Documents\EEG\LDA\LDA_All_fig\';
    saveas(gcf,strcat(save_folder,figstr,'.png'))
    close(5)
end


%% LDA code 2
% grasp_data = [];
% audio_data = [];
% grasp_labels = [];
% audio_labels = [];
% for j = 1:31
%     % LDA_filtered2\GMotor is for unfiltered data.
%     class1 = "Grasp";
%     class2 = "Audio";
%     Gstr = strcat('LDA\LDA_filtered\G2Motor',int2str(j),'.mat');
%     load(Gstr)
%     grasp_data = [grasp_data;data];
%     grasp_labels = [grasp_labels;repmat(class1,size(data,1),1)];
%     Astr = strcat('LDA\LDA_filtered\A2Motor',int2str(j),'.mat');
%     load(Astr)
%     audio_data = [audio_data;data];  
%     audio_labels = [audio_labels;repmat(class2,size(data,1),1)];
% end
% labels = [grasp_labels;audio_labels];
% motor = [grasp_data;audio_data];
% 
% Mdl = fitcdiscr(motor,labels);
% cvMdl = crossval(Mdl, 'kfold', 5); % 5-fold cross-validation
% cvError = kfoldLoss(cvMdl);
% [predictedLabels, scores] = kfoldPredict(cvMdl);
% predictedLabels = string(predictedLabels);
% 
% if j == 1
%     channel = j;
% else
%     channel = j + 1;
% end
% figure(5)
% confusionchart(labels,predictedLabels,'Normalization','row-normalized');
% figstr = strcat("19_08_All_Channels");%,int2str(channel));
% titlestr = strcat("19 08 All Channels");%,int2str(channel));
% title(titlestr)
% fontsize(18,"points")
% save_folder = 'C:\Users\czhe0008\Documents\EEG\LDA\LDA_filtered2_fig\';
% saveas(gcf,strcat(save_folder,figstr,'.png'))
% close(5)

%% LDA code
% % Selected features, all trials.
% for j = 30:30
%     % LDA_filtered2\GMotor is for unfiltered data.
%     Gstr = strcat('LDA\LDA_filtered\G2Motor',int2str(j),'.mat');
%     load(Gstr)
%     grasp_data = data;
%     Astr = strcat('LDA\LDA_filtered\A2Motor',int2str(j),'.mat');
%     load(Astr)
%     audio_data = data;
%     class1 = "Grasp";
%     class2 = "Audio";
%     % [predictedLabels,scores,cvError] = LDA(grasp_data,audio_data,class1,class2);
%     labels = [repmat(class1,size(grasp_data,1),1);repmat(class2,size(audio_data,1),1)];
%     motor = [grasp_data;audio_data];
% 
%     Mdl = fitcdiscr(motor,labels);
%     cvMdl = crossval(Mdl, 'kfold', 5); % 5-fold cross-validation
%     cvError = kfoldLoss(cvMdl);
%     [predictedLabels, scores] = kfoldPredict(cvMdl);
%     predictedLabels = string(predictedLabels);
% 
%     if j == 1
%         channel = j;
%     else
%         channel = j + 1;
%     end
%     figure(5)
%     confusionchart(labels,predictedLabels);
%     figstr = strcat("unfiltered_19_08_Channel_",int2str(channel));
%     titlestr = strcat("unfiltered 19 08 Channel ",int2str(channel));
%     title(titlestr)
%     fontsize(18,"points")
%     save_folder = 'C:\Users\czhe0008\Documents\EEG\LDA\LDA_filtered_fig\';
%     % saveas(gcf,strcat(save_folder,figstr,'.png'))
%     close(5)
% 
%     for i = 1:4 %19
%         for m = 1:4 %2
%             if m == i
%                 continue
%             end
%             figure(6)
%             hold on
%             f1 = motor(:,m);
%             f2 = motor(:,i);
%             Mdl = fitcdiscr([f1,f2],labels);
%             K = Mdl.Coeffs(1,2).Const;  
%             L = Mdl.Coeffs(1,2).Linear;
%             f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
%             h2 = fimplicit(f, [min(f1) max(f1) min(f2) max(f2)]);
%             h2.Color = 'k';
%             h2.LineWidth = 2;
%             h1 = gscatter(f1,f2,labels,'rb','ov',[],'off');
%             t = 13;
%             h3 = gscatter(f1(t),f2(t),labels(t),'g','*',[],'off');
%             titlestr = strcat("Unfiltered Ch30 Feature ",int2str(m), "_Feature ",int2str(i), ", highlight trial ",int2str(t));
%             xstr = strcat("Feature ",int2str(m));
%             ystr = strcat("Feature ",int2str(i));
%             title(titlestr)
%             xlabel(xstr)
%             ylabel(ystr)
%             save_folder = 'C:\Users\czhe0008\Documents\EEG\LDA\LDA_filtered_fig\';
%             saveas(gcf,strcat(save_folder,titlestr,'.png'))
%             close(6)
%         end
%     end
% end