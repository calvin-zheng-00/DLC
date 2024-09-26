close all; clear all; clc;

src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\data\high_thresh\angles';
pc_src = 'C:\Users\czhe0008\Documents\MATLAB\Openpose\PCA_coeffs\';
mean_src = 'C:\Users\czhe0008\Documents\MATLAB\Openpose\PCA_mean\';
pro_src = 'C:\Users\czhe0008\Documents\MATLAB\Openpose\Projection\';
time_src = 'C:\Users\czhe0008\Documents\MATLAB\Openpose\PCA_times\';
unpro_src = 'C:\Users\czhe0008\Documents\MATLAB\Openpose\Unprojected\';
sig_src = 'C:\Users\czhe0008\Documents\MATLAB\Openpose\Significance\';
allcoeff = 'C:\Users\czhe0008\Documents\MATLAB\Openpose\PCA_coeffs\pconcat_coeff.csv';
% minus9coeff = 'C:\Users\czhe0008\Documents\MATLAB\Openpose\PCA_coeffs\pminus9_coeff.csv';
bad_p = ['p1_high_angles','p2_high_angles'];
joints = ["Thumb carpormetacarpal flexion","Thumb metacarpophalangeal flexion","Thumb interphalangeal flexion","Index metacarpophalangeal flexion",...
    "Index proximal interphalangeal flexion","Index distal interphalangeal flexion","Middle metacarpophalangeal flexion","Middle proximal interphalangeal flexion",...
    "Middle distal interphalangeal flexion","Ring metacarpophalangeal flexion","Ring proximal interphalangeal flexion","Ring distal interphalangeal flexion",...
    "Little metacarpophalangeal flexion","Little proximal interphalangeal flexion","Little distal interphalangeal flexion","Elbow flexion",...
    "Index abduction","Middle abduction","Ring abduction","Little abduction","Thumb abduction","Thumb rotation", "Wrist flexion","Wrist Abduction",...
    "Wrist rotation","Shoulder flexion","Shoulder abduction","Shoulder rotation"];
colnames = categorical(joints);
colnames = reordercats(colnames,string(colnames));

allpcTable = readtable(allcoeff);
allpcarr = table2array(allpcTable(2:end,:));
% minus9pcTable = readtable(minus9coeff);
% minus9pcarr = table2array(minus9pcTable(2:end,:));
pctemp = dir(fullfile(pc_src,'*'));
pcfolder = setdiff({pctemp(~[pctemp.isdir]).name},{'.','..'});
meantemp = dir(fullfile(mean_src,'*'));
meanfolder = setdiff({meantemp(~[meantemp.isdir]).name},{'.','..'});
protemp = dir(fullfile(pro_src,'*'));
profolder = setdiff({protemp(~[protemp.isdir]).name},{'.','..'});
timetemp = dir(fullfile(time_src,'*'));
timefolder = setdiff({timetemp(~[timetemp.isdir]).name},{'.','..'});
unprotemp = dir(fullfile(unpro_src,'*'));
unprofolder = setdiff({unprotemp(~[unprotemp.isdir]).name},{'.','..'});
sigtemp = dir(fullfile(sig_src,'*'));
sigfolder = setdiff({sigtemp(~[sigtemp.isdir]).name},{'.','..'});

total_mean_corr = [];

% Looping through participants
for i = 1:numel(profolder)
    % Loading files
    participant = profolder{i}(1:3);
    file = fullfile(pro_src,profolder{i});
    fprintf(1, 'Now reading %s\n', profolder{i});
    proTable = readtable(file);
    proarr = table2array(proTable);
    file = fullfile(pc_src,pcfolder{i});
    pcTable = readtable(file);
    pcarr = table2array(pcTable(2:end,:));
    file = fullfile(mean_src,meanfolder{i});
    meanTable = readtable(file);
    meanarr = table2array(meanTable);
    file = fullfile(time_src,timefolder{i});
    timeTable = readtable(file);
    file = fullfile(unpro_src,unprofolder{i});
    unproTable = readtable(file);
    unproarr = table2array(unproTable);
    file = fullfile(sig_src,sigfolder{i});
    sig = matfile(file);
    sig = sig.curr_sig_80;
    sig_all = matfile('Significance\pAll_sig_80.mat');
    sig_all = sig_all.significant_80;
    
    %combined = [pcarr(:,1:sig), allpcarr(:,1:sig_all)];
    %[U,S,V] = svd(combined);
    %temp = proarr * U;

    %% New code
    unproarr(any(isnan(unproarr), 2), :) = [];
    temp_corr = [];
    for j = 8:8
        for k = 1:28
            reproarr = proarr(:,1:k)*pcarr(:,1:k)' + repmat(meanarr,height(proarr),1);
            reproarr(any(isnan(reproarr), 2), :) = [];
            figure 
            hold on
            plot((timeTable.start_time(1):timeTable.end_time(1))./40, reproarr(timeTable.start_time(1):timeTable.end_time(1),j).*(180/pi));
            plot((timeTable.start_time(1):timeTable.end_time(1))./40, unproarr(timeTable.start_time(1):timeTable.end_time(1),j).*(180/pi));
            legend('reproject','original')
            xlabel('time (seconds)')
            ylabel('angle (degrees)')
            %ylim([40 180]);
            title(strcat(joints(j)," PC", int2str(k)))
            saveas(gcf,strcat('PCA_images\middle\middle',int2str(k),'.png'))
            correlations = corrcoef(unproarr(timeTable.start_time(1):timeTable.end_time(1),j),reproarr(timeTable.start_time(1):timeTable.end_time(1),j));
            temp_corr = [temp_corr, correlations(1,2)];
        end
    end
    figure
    bar(1:28,temp_corr)
    title("Correlation between original and reprojected data")
    ylabel('Correlation')
    saveas(gcf,strcat('PCA_images\middle\middle_corr.png'))

    %%end of new code
    total_corr = [];
    mean_corr = [];
    % minus9_corr = [];
    num_pcs = width(unproarr);
    unproarr = normalize(unproarr,1);
    % Looping through PCs to reproject
    for j = 0:(num_pcs-1)
        reproarr = proarr(:,1:(num_pcs-j))*pcarr(:,1:(num_pcs-j))' + repmat(meanarr,height(proarr),1);
        reproarr(any(isnan(reproarr), 2), :) = [];
        unproarr(any(isnan(unproarr), 2), :) = [];
        correlations = corrcoef(unproarr,reproarr);
        total_corr = [total_corr, correlations(1,2)];

        meanreproarr = proarr(:,1:(num_pcs-j))*allpcarr(:,1:(num_pcs-j))' + repmat(meanarr,height(proarr),1);
        meanreproarr(any(isnan(meanreproarr), 2), :) = [];
        correlations = corrcoef(unproarr,meanreproarr);
        mean_corr = [mean_corr, correlations(1,2)];

        % if i == numel(profolder)
        %     min9reproarr = proarr(:,1:(num_pcs-j))*minus9pcarr(:,1:(num_pcs-j))' + repmat(meanarr,height(proarr),1);
        %     min9reproarr(any(isnan(min9reproarr), 2), :) = [];
        %     correlations = corrcoef(unproarr,min9reproarr);
        %     minus9_corr = [minus9_corr, correlations(1,2)];
        % end

    end

    total_mean_corr = [total_mean_corr;mean_corr];

    figure
    bar(flip(total_corr))
    titlestr = strcat("Correlation between original and reprojected data for ", participant);
    title(titlestr)
    xlabel('Number of principal components used in reprojection')
    ylabel('Correlation')

    figure
    bar(flip(mean_corr))
    titlestr = strcat("Correlation between original and generalised PC reprojected data for ", participant);
    title(titlestr)
    xlabel('Number of principal components used in reprojection')
    ylabel('Correlation')
    
    % if i == numel(profolder)
    %     figure
    %     bar(flip(minus9_corr.^2))
    %     title("Reprojection accuracy of generalised principal components")
    %     xlabel('Number of principal components used in reprojection')
    %     ylabel('Percentage of movement variance explained by reprojection')
    %     set(gca,'fontsize',14, 'TickDir', 'out')
    %     box off
    %     ylim([0 1])
    % end
end

figure
hold on
mean_var = flip(mean(total_mean_corr.^2, 1));
bar(mean_var, 'w')
% SEM = std(total_mean_corr.^2,0,1)/sqrt(size(total_mean_corr.^2,1));
% er = errorbar(mean_var, flip(SEM));
% er.Color = [0 0 0];
% er.LineStyle = 'none';
%title("Average correlation between original and generalised PC reprojected data")
%title("Reprojection accuracy of generalised principal components")
xlabel('Number of principal components used in reprojection')
ylabel('Percentage of variance explained in reprojection')
set(gca,'fontsize',14, 'TickDir', 'out')
box off
ylim([0 1])

% pc_path = strcat(pc_src,coeff);
% 
% pc = readtable(pc_path);
% pc = table2array(pc);
% pc = pc(2:end,:);

% maintemp = dir(fullfile(src,'*'));
% mainfolder = setdiff({maintemp([maintemp.isdir]).name},{'.','..'});
% % for i = 1:numel(mainfolder)
% for i = 5:5
%     if contains(bad_p,mainfolder{i})
%         continue
%     end
%     participantAngles = [];
%     subtemp = dir(fullfile(src,mainfolder{i},'*.csv'));
%     subfolder = {subtemp(~[subtemp.isdir]).name};
%     %for j = 1:numel(subfolder)
%     for j = 33:33
%         file = fullfile(src,mainfolder{i},subfolder{j});
%         fprintf(1, 'Now reading %s\n', subfolder{j});
%         thisTable = readtable(file);
%         participantAngles = [participantAngles;thisTable];
%     end
%     anglesarr = table2array(participantAngles);
%     repo_angles = anglesarr * pc;
% 
%     figure
%     hold on
%     plot(1:height(anglesarr), anglesarr(:,16));
%     legend('flex','abd','rot')
%     ylabel('Angle (radians)')
%     xlabel('timestep')
%     title(mainfolder{i})
%     set(gca,'fontsize',14, 'TickDir', 'out')
%     box off
% 
%     figure
%     hold on
%     plot(1:height(anglesarr), repo_angles(:,5));
%     ylabel('PC activation level')
%     xlabel('timestep')
%     title(mainfolder{i})
%     set(gca,'fontsize',14, 'TickDir', 'out')
%     box off
% end