close all; clear all; clc;

category = "All";
%src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\data\high_thresh\angles';
src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\angles_participant';


allAngles = [];
p_sig_80 = [];
p_latent = [];
allExplained = [];
allCoeff = [];
maintemp = dir(fullfile(src,'*'));
participants = setdiff({maintemp([maintemp.isdir]).name},{'.','..'});
acttemp = dir(fullfile(src,'p23','*'));
actfolder = setdiff({acttemp([acttemp.isdir]).name},{'.','..'});
% Looping through activities
for i = 1:numel(actfolder)
    activityAngles = [];
    activities = [];
    start_time = [];
    end_time = [];
    % Looping through participants
    for j = 1:numel(participants)
        subsubtemp = dir(fullfile(src,participants{j},actfolder{i},'*.csv'));
        trials = {subsubtemp(~[subsubtemp.isdir]).name};
        current_file = strcat(participants{j},actfolder{i});
        fprintf(1, 'Now reading %s\n', current_file);
        % Concatonating angles from multiple trials
        for k = 1:numel(trials)
            file = fullfile(src,participants{i},activities{j},trials{k});
            thisTable = readtable(file);
            activityAngles = [activityAngles;thisTable];
        end
    end
    activityAngles = table2array(activityAngles);

    %% Current WIP
    % activityAngles = normalize(activityAngles,1);
    [coeff,score,latent,tsquared,explained,mu] = pca(activityAngles, 'Rows', 'pairwise');
    curr_sig_80 = find(cumsum(explained)>80,1); % Finds the amount of PCs needed to explain 80% of variance
    allExplained = [allExplained,explained];
    allCoeff(:,:,i) = coeff;
    save(strcat('Significance\',participants{i},'_sig_80.mat'),'curr_sig_80');

    % The score matrix contains the projected data
    table_projection = array2table(score);
    writetable(table_projection, strcat('Projection\',participants{i},'_projected.csv'))
    % The mu matrix contains the estimated mean of each variable, needed
    % for reprojection
    table_mean = array2table(mu);
    writetable(table_mean, strcat('PCA_mean\',participants{i},'_mean.csv'))
    % Saving the raw data
    writetable(activityAngles, strcat('Unprojected\',participants{i},'.csv'))

    % The PCA coefficients
    table_coeff = array2table(coeff,'VariableNames',string(1:length(coeff)));
    
    figure
    bar(explained, 'w');
    %title(category + " principal component percentage variance")
    ylabel('Percentage of variance explained')
    xlabel('Principal Component')
    titlestr = strcat("PCA of grasp type: ", participants{i});
    title(titlestr, 'Interpreter', 'none')
    set(gca,'fontsize',14, 'TickDir', 'out')
    ylim([0 80]);
    box off
    writetable(table_coeff, strcat('PCA_coeffs\',participants{i},'_coeff.csv'))
    saveas(gcf,strcat('PCA_images\',participants{i},'_PCA.png'))

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

angles = table2array(allAngles);
%angles = normalize(angles,1);
[coeff,score,latent,tsquared,explained,mu] = pca(angles, 'Rows', 'pairwise');
%[coeff,score,latent,tsquared,explained,mu] = pca(angles, 'Rows', 'complete');
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
save(strcat('Significance\p',category,'_sig_80.mat'),'significant_80');

table_coeff = array2table(coeff,'VariableNames',string(1:length(coeff)));

figure
bar(explained, 'w');
%title(category + " principal component percentage variance")
ylabel('Percentage of movement variance explained')
xlabel('Principal Component')
title('Concatenated PCs')
set(gca,'fontsize',14, 'TickDir', 'out')
ylim([0 80])
box off
writetable(table_coeff, 'PCA_coeffs\pconcat_coeff.csv')
saveas(gcf,'PCA_images\concat_'+category+'_PCA.png')

% Write (specifically explained) to outout
% PC cutoff
% How to determine if two different task sets share PCAs?
% Convert pixel locations into joint angles first, or is there a way to
% combine the x,y,z coords of every point?

%Calculating Average PCs of every participant
mean_exp = mean(allExplained, 2);
SEM = std(allExplained,0,2)/sqrt(size(allExplained,2));

figure
bar(mean_exp, 'w');
hold on;
er = errorbar(mean_exp, SEM);
er.Color = [0 0 0];
er.LineStyle = 'none';
%title(category + " principal component percentage variance")
ylabel('Percentage of movement variance explained')
xlabel('Principal Component')
title('Mean PCs')
set(gca,'fontsize',14, 'TickDir', 'out')
ylim([0 80])
box off
saveas(gcf,'PCA_images\average_'+category+'_PCA.png')

% Calculating Correlations
correlations = [];
for i = 1:size(allCoeff,3)
    if sum(sum(allCoeff(:,:,i))) == 0
        continue
    end
    temp_corr = [];
    for j = 1:size(allCoeff,2)
        temp = corrcoef(coeff(:,j),allCoeff(:,j,i));
        temp_corr = [temp_corr, temp(1,2)];
    end
    correlations = [correlations;temp_corr];
end

figure
bar(transpose(correlations(:, 1:6)));
hold on;
ylabel('Correlation Coefficient')
xlabel('Principal Component')
title('PC (concat) correlations')
set(gca,'fontsize',14, 'TickDir', 'out')
% legend('P10','P11','P12')
box off
saveas(gcf,'PCA_images\correlation_'+category+'_PCA.png')

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