close all; clear all; clc;

category = "All";
src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\categorized_angles_2';
src_3d = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\categorized';

allAngles = [];
grasp_angles = [];
p_sig_80 = [];
p_latent = [];
allExplained = [];
allCoeff = [];
maintemp = dir(fullfile(src,'*'));
mainfolder = setdiff({maintemp([maintemp.isdir]).name},{'.','..'});
% Looping through categories
for i = 1:numel(mainfolder)
    categoryAngles = [];
    subtemp = dir(fullfile(src,mainfolder{i},'*.csv'));
    subfolder = {subtemp(~[subtemp.isdir]).name};
    postemp = dir(fullfile(src_3d,mainfolder{i},'*.csv'));
    posfolder = {postemp(~[postemp.isdir]).name};
    activities = [];
    start_time = [];
    end_time = [];
    grasp_time = [];
    % Concatonating angles from multiple activities
    for j = 1:numel(subfolder)
        file = fullfile(src,mainfolder{i},subfolder{j});
        fprintf(1, 'Now reading %s\n', subfolder{j});
        thisTable = readtable(file);

        [filepath,name,ext] = fileparts(file);
        activities = [activities;convertCharsToStrings(name(1:end-7))];

        % Getting moment of grasp
        grasp_3d = fullfile(src_3d,mainfolder{i},posfolder{j});
        grasp_3d = table2array(readtable(grasp_3d));
        vel = grasp_3d(1:end-1,:)-grasp_3d(2:end,:);
        vel = sqrt((vel(:,1:3:end)).^2 + (vel(:,2:3:end)).^2 + (vel(:,3:3:end)).^2);
        [~, vel_rise] = max(vel(1:end-1, 26) > 10 & vel(2:end, 26) > 10);
        [~, vel_fall] = max(vel(vel_rise:end-1, 26) < 10 & vel(vel_rise+1:end, 26) < 10);
        reach = vel_fall + vel_rise + 6;
        grasp_time = [grasp_time;reach];

        % Finding start and end times of each activity
        if j == 1
            start_time = [1];
        else
            start_time = [start_time;end_time(end)+1];
        end
        %end_time = [end_time;start_time(end)+height(thisTable)-1];
        end_time = [end_time;start_time(end)+min([reach,40])-1];

        %temp = table2array(thisTable);
        %grasp_angles = cat(1,grasp_angles,temp(max([1,reach-40]):reach,:));
        categoryAngles = [categoryAngles;thisTable(max([1,reach-39]):reach,:)];

    end
    allAngles = [allAngles;categoryAngles];

    % Keeping tracking of the start and end times of each activity after
    % concatenation
    table_times = table(activities,start_time,end_time,grasp_time);
    writetable(table_times, strcat('PCA_times\',mainfolder{i},'_times.csv'))

    pAngles = table2array(categoryAngles);

    % Correlation
    % total_corr = [];
    % for j = 1:size(pAngles,2)
    %     mid_corr = [];
    %     for k = 1:size(pAngles,2)
    %         correlations = corrcoef(pAngles(:,j),pAngles(:,k),'rows','complete');
    %         mid_corr = [mid_corr, correlations(1,2)];
    %     end
    %     total_corr = [total_corr;mid_corr];
    % end

    %% PCA uncomment after variance
    %pAngles = normalize(pAngles,1);
    [coeff,score,latent,tsquared,explained,mu] = pca(pAngles, 'Rows', 'pairwise'); %pAngles
    p_latent = [p_latent,sum(explained(find(latent>1)))];
    p_sig_80 = [p_sig_80,find(cumsum(explained)>80,1)]; % Finds the amount of PCs needed to explain 80% of variance
    curr_sig_80 = find(cumsum(explained)>80,1);
    allExplained = [allExplained,explained];
    allCoeff(:,:,i) = coeff;
    save(strcat('Significance\',mainfolder{i},'_sig_80.mat'),'curr_sig_80');
    % example = matfile('saveA.mat');
    % C = example.A;

    % The score matrix contains the projected data
    table_projection = array2table(score);
    writetable(table_projection, strcat('Projection\',mainfolder{i},'_projected.csv'))
    % The mu matrix contains the estimated mean of each variable, needed
    % for reprojection
    table_mean = array2table(mu);
    writetable(table_mean, strcat('PCA_mean\',mainfolder{i},'_mean.csv'))
    % Saving the raw data
    writetable(categoryAngles, strcat('Unprojected\',mainfolder{i},'.csv'))

    % The PCA coefficients
    table_coeff = array2table(coeff,'VariableNames',string(1:length(coeff)));

    figure
    bar(explained, 'w');
    %title(category + " principal component percentage variance")
    ylabel('Percentage of variance explained')
    xlabel('Principal Component')
    titlestr = strcat("PCA of grasp type: ", mainfolder{i});
    title(titlestr, 'Interpreter', 'none')
    set(gca,'fontsize',14, 'TickDir', 'out')
    ylim([0 80]);
    box off
    writetable(table_coeff, strcat('PCA_coeffs\',mainfolder{i},'_coeff.csv'))
    saveas(gcf,strcat('PCA_images\',mainfolder{i},'_PCA.png'))

end

% theFiles = dir(fullfile(src, '*.csv'));
% baseFileName = theFiles(1).name;
% fullFileName = fullfile(theFiles(1).folder, baseFileName);
% fprintf(1, 'Now reading %s\n', fullFileName);
% allAngles = readtable(fullFileName);
% for k = 2 : length(theFiles)
%     baseFileName = theFiles(k).name;
%     fullFileName = fullfile(theFiles(k).folder, baseFileName);
%     fprintf(1, 'Now reading %s\n', fullFileName);
%     thisTable = readtable(fullFileName);
%     allAngles = [allAngles;thisTable];
% end

%% PCA just for grasp
% [coeff,score,latent,tsquared,explained,mu] = pca(grasp_angles, 'Rows', 'pairwise');
% 
% idx = find(cumsum(explained)>95,1);
% significant_80 = find(cumsum(explained)>80,1);
% idx2 = find(latent>1);
% var_ex = sum(explained(idx2));
% 
% % determining if two different task sets share PCAs?
% % newData is the new set of data I want to apply the existing PCA onto.
% %newDataPCA = (newData-mu)*coeff(:,1:idx);
% % Remove idx if I want all the components
% 
% save('All\pAll_sig_80.mat','significant_80');
% save("All\pconcat_coeff.mat","coeff")
% save("All\pconcat_project.mat","score")
% save("All\pconcat_mean.mat","mu")
% 
% figure
% bar(explained, 'w');
% %title(category + " principal component percentage variance")
% ylabel('Percentage of movement variance explained')
% xlabel('Principal Component')
% title('Concatenated PCs')
% set(gca,'fontsize',14, 'TickDir', 'out')
% ylim([0 80])
% box off
% %saveas(gcf,'PCA_images\concat_'+category+'_PCA.png')
% saveas(gcf,'reach_and_grasp\concat_'+category+'_PCA.png')

%% Original PCA for all activities
angles = table2array(allAngles);
%angles = normalize(angles,1);
[coeff,score,latent,tsquared,explained,mu] = pca(angles, 'Rows', 'pairwise');
% coeff are the principal component coefficients, also known as loadings
%           (eigenvectors)
% score are the principal component scores
% latent are the principal component variances (eigenvalues)
% explained are the percentage of the total variance explained by each principal component
% mu are the estimated mean of each variable in X
% Name-Value arguments:
% Example: 'Algorithm','eig','Centered',false,'Rows','all','NumComponents',3
% specifies that pca uses eigenvalue decomposition algorithm, not center the data,
% use all of the observations, and return only the first three principal components.

% Find the number of components required to explain at least 95%
% variability. (idx) (set cut off two 95% and 85%)
idx = find(cumsum(explained)>95,1);
significant_85 = find(cumsum(explained)>85,1);
significant_80 = find(cumsum(explained)>80,1);
idx2 = find(latent>1);
var_ex = sum(explained(idx2));

% determining if two different task sets share PCAs?
% newData is the new set of data I want to apply the existing PCA onto.
%newDataPCA = (newData-mu)*coeff(:,1:idx);
% Remove idx if I want all the components
save('All\pAll_sig_80.mat','significant_80');
save("All\pconcat_coeff.mat","coeff")
save("All\pconcat_project.mat","score")
save("All\pconcat_mean.mat","mu")

figure
bar(explained, 'w');
%title(category + " principal component percentage variance")
ylabel('Percentage of movement variance explained')
xlabel('Principal Component')
title('Concatenated PCs')
set(gca,'fontsize',14, 'TickDir', 'out')
ylim([0 80])
box off
saveas(gcf,'PCA_images\concat_'+category+'_PCA.png')

% Write (specifically explained) to outout
% PC cutoff
% How to determine if two different task sets share PCAs?
% Convert pixel locations into joint angles first, or is there a way to
% combine the x,y,z coords of every point?

%% Calculating Average PCs of every participant
% mean_exp = mean(allExplained, 2);
% SEM = std(allExplained,0,2)/sqrt(size(allExplained,2));
% 
% figure
% bar(mean_exp, 'w');
% hold on;
% er = errorbar(mean_exp, SEM);
% er.Color = [0 0 0];
% er.LineStyle = 'none';
% %title(category + " principal component percentage variance")
% ylabel('Percentage of movement variance explained')
% xlabel('Principal Component')
% title('Mean PCs')
% set(gca,'fontsize',14, 'TickDir', 'out')
% ylim([0 80])
% box off
% saveas(gcf,'PCA_images\average_'+category+'_PCA.png')

%% Calculating Correlations
% correlations = [];
% for i = 1:size(allCoeff,3)
%     if sum(sum(allCoeff(:,:,i))) == 0
%         continue
%     end
%     temp_corr = [];
%     for j = 1:size(allCoeff,2)
%         temp = corrcoef(coeff(:,j),allCoeff(:,j,i));
%         temp_corr = [temp_corr, temp(1,2)];
%     end
%     correlations = [correlations;temp_corr];
% end
% 
% figure
% bar(transpose(correlations(:, 1:6)));
% hold on;
% ylabel('Correlation Coefficient')
% xlabel('Principal Component')
% title('PC (concat) correlations')
% set(gca,'fontsize',14, 'TickDir', 'out')
% % legend('P10','P11','P12')
% box off
% saveas(gcf,'PCA_images\correlation_'+category+'_PCA.png')


%% Old code
% angles = table2array(minus9Angles);
% %angles = normalize(angles,1);
% [coeff,score,latent,tsquared,explained,mu] = pca(angles, 'Rows', 'pairwise');
% table_coeff = array2table(coeff,'VariableNames',string(1:length(coeff)));
% 
% figure
% bar(explained, 'w');
% %title(category + " principal component percentage variance")
% ylabel('Percentage of movement variance explained')
% xlabel('Principal Component')
% %title('Concatenated PCs')
% set(gca,'fontsize',14, 'TickDir', 'out')
% ylim([0 50])
% box off
% writetable(table_coeff, 'PCA_coeffs\pminus9_coeff.csv')

%% Calculating mean Correlations
% m_correlations = [];
% for i = 1:size(allCoeff,3)
%     if sum(sum(allCoeff(:,:,i))) == 0
%         continue
%     end
%     for k = (i+1):size(allCoeff,3)
%         temp_corr = [];
%         for j = 1:size(allCoeff,2)
%             temp = corrcoef(allCoeff(:,j,k),allCoeff(:,j,i));
%             temp_corr = [temp_corr, temp(1,2)];
%         end
%         figure
%         bar(transpose(temp_corr));
%         hold on;
%         ylabel('Correlation Coefficient')
%         xlabel('Principal Component')
%         %title('PC correlations ' + int2str(i) + '/' + int2str(k))
%         title('PC individual correlations')
%         set(gca,'fontsize',14, 'TickDir', 'out')
%         box off
%         %writetable(table_coeff, 'PCA_coeffs\correlation_'+category+'_coeff.csv')
%         %saveas(gcf,'PCA_images\correlation_'+category+'_PCA.png')
%     end
%     %m_correlations = [m_correlations;temp_corr];
% end