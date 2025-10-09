function PCA_filt(src,dest,threshold,increment)
    % Filter input data to remove bad trials.
    % src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\3d';
    % dest = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\filter';
    % threshold = 30;   % When to start flagging jumps
    % increment = 20;  % How fast to expand acceptable range
    
    parttemp = dir(fullfile(src,'*'));
    partfolder = setdiff({parttemp([parttemp.isdir]).name},{'.','..'});

    for participant_i = 1:numel(partfolder)
        % Extra loop incase of categorised 3d data
        acttemp = dir(fullfile(src,partfolder{participant_i},'*'));
        actfolder = setdiff({acttemp([acttemp.isdir]).name},{'.','..'});
        for activity_k = 1:numel(actfolder)
            trialtemp = dir(fullfile(src,partfolder{participant_i},actfolder{activity_k},'*'));
            trialfolder = setdiff({trialtemp([trialtemp.isdir]).name},{'.','..'});
            for trial_j = 1:numel(trialfolder)
                if isfile(dir(fullfile(src,partfolder{participant_i},actfolder{activity_k},trialfolder{trial_j},'df_3d.csv'))) && isfile(dir(fullfile(src,partfolder{participant_i},actfolder{activity_k},trialfolder{trial_j},'df_err.csv')))
                    % Getting 3d position values
                    data = fullfile(src,partfolder{participant_i},actfolder{activity_k},trialfolder{trial_j},'df_3d.csv');
                    df = readtable(data);
                    % Getting error data
                    data_err = fullfile(src,partfolder{participant_i},actfolder{activity_k},trialfolder{trial_j},'df_err.csv');
                    df_err = readtable(data_err);
                    % Removing the timestep ends since those sections have
                    % bad tracking.
                    df([1:5,end-4:end],:) = [];
                    df_err([1:5,end-4:end],:) = [];
                    % Removing unused markers
                    df(:,[1:3,79:end]) = [];
                    df_err = df_err(:,2:26);
                    df2_check = 1;
                    for joint_idx = 1:width(df)
                        col = df.(joint_idx);
                        diff = threshold;
                        check = 1;
                        col = fillmissing(col,'linear',1,'EndValues','nearest');
                        % Checking for large jumps
                        for i = 2:length(col)
                            if abs(col(i) - col(check)) > diff
                                diff = diff + increment;
                            else
                                grad = (col(i)-col(check))/(i-check);
                                for j = check+1:i
                                    col(j) = col(j-1) + grad;
                                end
                                diff = threshold;
                                check = i;
                            end
                        end
                        % Butterworth
                        fc = 5;
                        fs = 40;    
                        [b, a] = butter(2,fc/(fs/2));
                        filtered = filtfilt(b,a,col);

                        colname = df.Properties.VariableNames{joint_idx};
                        if df2_check == 1
                            df2 = table(filtered);
                            df2.Properties.VariableNames = convertCharsToStrings(colname);
                            df2_check = 0;
                        else
                            df2.(colname) = filtered;
                        end
                    end
                    writetable(df2,fullfile(dest,strcat(actfolder{activity_k},"_filtered.csv")))
                end
            end
        end
    end
end

function PCA_participant()
    % Finding PCs from individuals.
end

function cosine_cluster(X)
%% Example: Clustering Hand Synergies Across Subjects
    % Assume each column of 'synergies' is a synergy vector
    % Rows = joints, Columns = different synergies from all subjects
    % Example: 15 joints, 12 synergies total (4 subjects × 3 synergies each)
    synergies = X;  % replace with your actual data
    
    %% Step 1: Compute pairwise cosine distance
    % cosine distance = 1 - cosine similarity
    cosDist = pdist(synergies', 'cosine');  % transpose to make synergies as rows
    
    %% Step 2: Cluster synergies using hierarchical clustering
    Z = linkage(cosDist, 'average');  % hierarchical linkage method
    
    % Plot dendrogram to visualize clustering
    figure;
    dendrogram(Z);
    title('Hierarchical Clustering of Synergies (Cosine Distance)');
    
    %% Step 3: Assign clusters
    numClusters = 3;  % choose number of general synergies you expect
    clusterIdx = cluster(Z, 'maxclust', numClusters);
    
    % Display which synergy belongs to which cluster
    disp('Synergy cluster assignments:');
    disp(clusterIdx);
    
    %% Step 4: Compute cluster centroids (general synergies)
    generalSynergies = zeros(size(synergies,1), numClusters);
    for k = 1:numClusters
        generalSynergies(:,k) = mean(synergies(:, clusterIdx == k), 2);
    end
    
    % Optional: normalize general synergies to unit length
    generalSynergies = generalSynergies ./ vecnorm(generalSynergies);
    
    disp('General synergies (cluster centroids):');
    disp(generalSynergies);
end

function PCA_cluster(X)
    % clustering PCs.
    %X = abs(X);
    Y = pdist(X,"cosine");
    Z = linkage(Y,"average");
    dendrogram(Z,72)
    xlabel('Synergies')
    ylabel('Cosine Distance (Pi Rad)')
    title('Synergy clustering')
    C = cophenet(Z,Y);
    I = inconsistent(Z);
    T = cluster(Z,"cutoff",1.1);
    % T = cluster(Z,"maxclust",2)
    % Or
    T = clusterdata(X,Cutoff=cutoff);
    T = clusterdata(X,MaxClust=maxclust);
end

close all; clear all; clc;
src = 'C:\Users\czhe0008\Documents\MATLAB\Openpose\per_par\PCA_coeffs';
allPCs = [];
parttemp = dir(fullfile(src,'*.csv'));
partfolder = {parttemp(~[parttemp.isdir]).name};
for participant_i = 1:numel(partfolder)
    data = fullfile(src,partfolder{participant_i});
    df = readtable(data);
    df = table2array(df);
    df = df(2:end,1:4);
    allPCs = [allPCs,df];
end
PCA_cluster(allPCs');

% % Filter input data to remove bad trials.
% src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\3d';
% dest = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\filter';
% threshold = 30;   % When to start flagging jumps
% increment = 20;  % How fast to expand acceptable range
% 
% parttemp = dir(fullfile(src,'*'));
% partfolder = setdiff({parttemp([parttemp.isdir]).name},{'.','..'});
% 
% total_count = [];
% for participant_i = 1:numel(partfolder)
%     fprintf("participant %d\n",participant_i)
%     % Extra loop incase of categorised 3d data
%     par_count = [];
%     acttemp = dir(fullfile(src,partfolder{participant_i},'*'));
%     actfolder = setdiff({acttemp([acttemp.isdir]).name},{'.','..'});
%     for activity_k = 1:numel(actfolder)
%         trialtemp = dir(fullfile(src,partfolder{participant_i},actfolder{activity_k},'*'));
%         trialfolder = setdiff({trialtemp([trialtemp.isdir]).name},{'.','..'});
%         trial_count = 0;
%         for trial_j = 1:numel(trialfolder)
%             if isfile(fullfile(src,partfolder{participant_i},actfolder{activity_k},trialfolder{trial_j},'df_3d.csv')) && isfile(fullfile(src,partfolder{participant_i},actfolder{activity_k},trialfolder{trial_j},'df_err.csv'))
%                 data_err = fullfile(src,partfolder{participant_i},actfolder{activity_k},trialfolder{trial_j},'df_err.csv');
%                 df_err = readtable(data_err);
%                 df_err([1:5,end-4:end],:) = [];
%                 df_err = df_err(:,2:26);
%                 df_err = table2array(df_err);
%                 df_err_comp = df_err;
%                 df_err_comp_size = sum(sum(~isnan(df_err)));
%                 df_err(df_err>20) = NaN;
%                 df_err_size = sum(sum(~isnan(df_err)));
%                 if df_err_size > df_err_comp_size*0.95
%                     trial_count = trial_count + 1;
%                 end
%             end
%         end
%         par_count = [par_count;trial_count];
%     end
%     total_count = [total_count,par_count];
% end
