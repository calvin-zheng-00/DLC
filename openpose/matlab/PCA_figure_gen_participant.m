% This file conducts PCA on the angle data

close all; clear all; clc;

category = "All";
% src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\categorized_angles';
% src3d = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\categorized';
src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\angles_participant';
src3d = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\3d';
grasp_type = readtable('C:\Users\czhe0008\Documents\MATLAB\Openpose\grasp_type_3.csv');
grasp_type = sortrows(grasp_type, "Var2");

error_counter1 = [];
error_counter2 = [];
allAngles = [];
allLatent = [];
allExplained = [];
allCoeff = [];
categoryAngles = [];
maintemp = dir(fullfile(src,'*'));
mainfolder = setdiff({maintemp([maintemp.isdir]).name},{'.','..'});
% Looping through participants
for i = 1:numel(mainfolder)
    acttemp = dir(fullfile(src,mainfolder{i},'*'));
    actfolder = setdiff({acttemp([acttemp.isdir]).name},{'.','..'});
    % Looping through activities
    for l = 1:numel(actfolder)
        subtemp = dir(fullfile(src,mainfolder{i},actfolder{l},'*.csv'));
        subfolder = {subtemp(~[subtemp.isdir]).name};
        sub3dtemp = dir(fullfile(src3d,mainfolder{i},actfolder{l},'*'));
        sub3dfolder = setdiff({sub3dtemp([sub3dtemp.isdir]).name},{'.','..'});
        for j = 1:numel(subfolder)
            % Concatonating angles from multiple activities
            file = fullfile(src,mainfolder{i},actfolder{l},subfolder{j});
            if ~exist(file)
                continue
            end
            fprintf(1, 'Now reading %s\n', subfolder{j});
            thisTable = readtable(file);
            thisTable(end,:) = [];
            thisTable(1,:) = [];
            % Getting moment of grasp
            error_3d = fullfile(src3d,mainfolder{i},actfolder{l},sub3dfolder{j},'df_err.csv');
            if ~exist(error_3d)
                continue
            end
            error_3d = table2array(readtable(error_3d));
            if mean(mean(error_3d(2:end-1,1:21))) > 40
                continue
            end
            % This section removes the arm
            % thisTable(:,22:28) = [];
            % categoryAngles = [categoryAngles;thisTable];
            grasp_3d = fullfile(src3d,mainfolder{i},actfolder{l},sub3dfolder{j},'df_3d.csv');
            grasp_3d = table2array(readtable(grasp_3d));
            vel = grasp_3d(1:end-1,:)-grasp_3d(2:end,:);
            vel = sqrt((vel(:,1:3:end)).^2 + (vel(:,2:3:end)).^2 + (vel(:,3:3:end)).^2);
            [~, vel_rise] = max(vel(1:end-1, 25) > 10 & vel(2:end, 25) > 10);
            [~, vel_fall] = max(vel(vel_rise:end-1, 25) < 10 & vel(vel_rise+1:end, 25) < 10);
            reach = vel_fall + vel_rise + 6;
            end_temp = min(reach+40,height(thisTable));
            thisTable(:,22:28) = [];
            categoryAngles = [categoryAngles;thisTable(max([1,reach-39]):end_temp,:)];



            % grasp_3d = fullfile(src3d,mainfolder{i},sub3dfolder{j});
            % grasp_3d = fullfile('C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\categorized\Instructions_COO4\20241101T101232-101239_filtered.csv');
            % grasp_3d = table2array(readtable(grasp_3d));
            % grasp_3d(end,:) = [];
            % grasp_3d(1,:) = [];
            % reach_3d = grasp_3d(:,75);
            % resting = mean(reach_3d(end-10:end))-30;
            % over_threshold = find(reach_3d<resting)';
            % over_threshold2 = double(diff(over_threshold)==1);
            % reach = strfind(over_threshold2,[1,1,1,1,1,1,1,1,1,1]);
            % if isempty(reach)
            %     error_counter1 = [error_counter1,sub3dfolder{j}];
            %     continue;
            % end
            % reach = over_threshold(reach(1))-5;
            % if (reach + 80) > height(thisTable) || reach < 10
            %     error_counter2 = [error_counter2,sub3dfolder{j}];
            %     continue;
            % end
            % thisTable(:,22:28) = [];
            % categoryAngles = [categoryAngles;thisTable(reach:reach+80,:)];
        end
    end
    pAngles = table2array(categoryAngles);
    allAngles = [allAngles;categoryAngles];
    [coeff,score,latent,tsquared,explained,mu] = pca(pAngles, 'Rows', 'pairwise'); %pAngles
    allLatent = [allLatent, latent];
    allExplained = [allExplained,explained];
    allCoeff = cat(3,allCoeff,coeff);
    categoryAngles = [];
end

save("All\allExplained.mat","allExplained")
save("All\allLatent.mat","allLatent")
save("All\allCoeff.mat","allCoeff")

outall = cumsum(allExplained);

list95 = [];
list80 = [];
for i = 1:width(outall)
    list95 = [list95,find(outall(:,i)>95, 1)];
    list80 = [list80,find(outall(:,i)>80, 1)];
end

figure
plot(outall)
xlim([1 28])
ylim([0 100])
xlabel('Number of Principal Components')
ylabel('Cumulative % of Variance Explained')
title('Cumulative Variance Explained')
% legend({'Power Palm 1','Power Pad 1','Power Pad 2','Power Pad 3','Precision Pad 1','Precision Pad 2','Precision Pad 3','Precision Side 1','Power Palm 2','Power Palm 3','Intermediate Side 1','Finger Press','Palm press'}, 'Location','southeast');
legend({'Bathing','Cleaning','Cooking','Digital','Dining','Dressing','General','Mobility','Labor','Office'}, 'Location','southeast');
box off
set(gcf, 'Position', [100 100 300 400])

%% Original PCA for all activities
angles = table2array(allAngles);
[coeff,score,latent,tsquared,explained,mu] = pca(angles, 'Rows', 'pairwise');

% save("All\pconcat_coeff.mat","coeff")
% save("All\pconcat_project.mat","score")
% save("All\pconcat_mean.mat","mu")
% save("All\pconcat_explained.mat","explained")
% save("All\pconcat_latent.mat","latent")
% save("All\pconcat_tsquared.mat","tsquared")

figure
bar(explained, 'w');
ylabel('Percentage of movement variance explained')
xlabel('Principal Component')
title('Principal Components Of The Hand')
set(gca,'fontsize',14, 'TickDir', 'out')
ylim([0 80])
set(gcf, 'Position', [100 100 500 500])
box off