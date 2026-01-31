function PCA_testing(src,src_3d,dest,ALLEEG)
    %% Finding PCs from participants.
    % src - folder location of all the angle data
    % src_3d - folder location of all the positional data
    % dest - output file locations
    data_set = 1;
    loops_event = size(ALLEEG(data_set).event);
    loops_event = loops_event(2);
    event_time = [];
    
    for i = 2:loops_event
        if strcmp(ALLEEG(data_set).event(i).type, 'S  3')
            temp = ALLEEG(data_set).event(i).latency;
            temp = round((temp*40)/1000);
            temp2 = ALLEEG(data_set).event(i-1).latency;
            temp2 = round((temp2*40)/1000);
            event_time = [event_time,(temp-temp2)];
        end
    end
    % allAngles = [];
    categoryAngles = [];
    angletemp = dir(fullfile(src,'*.csv'));
    anglefolder = {angletemp(~[angletemp.isdir]).name};
    postemp = dir(fullfile(src_3d,'*.csv'));
    posfolder = {postemp(~[postemp.isdir]).name};
    % Looping through trials
    for i = 1:numel(anglefolder)
        
        activities = [];
        start_time = 1;
        reach_time = [];
        end_time = [];
        grasp_time = [];

        file = fullfile(src,anglefolder{i});
        fprintf(1, 'Now reading %s\n', anglefolder{i});
        thisTable = readtable(file);

        [filepath,name,ext] = fileparts(file);
        activities = [activities;convertCharsToStrings(name(1:end-7))];

        % Getting moment of grasp
        % grasp_3d = fullfile(src_3d,posfolder{i});
        % grasp_3d = table2array(readtable(grasp_3d));
        % vel = grasp_3d(1:end-1,:)-grasp_3d(2:end,:);
        % vel = sqrt((vel(:,1:3:end)).^2 + (vel(:,2:3:end)).^2 + (vel(:,3:3:end)).^2);
        % [~, vel_rise] = max(vel(1:end-1, 25) > 10 & vel(2:end, 25) > 10);
        % [~, vel_fall] = max(vel(vel_rise:end-1, 25) < 10 & vel(vel_rise+1:end, 25) < 10);
        % reach = vel_fall + vel_rise + 6;
        % grasp_time = [grasp_time;reach];
        % 
        % % Finding start and end times of each activity
        % if ~(i == 1)
        %     start_time = [start_time;end_time(end)+1];
        % end
        % reach_time = [reach_time;start_time(end)+min([reach,40])-1];
        % end_temp = min(reach+10,height(thisTable));
        % end_time = [end_time;reach_time(end)+min(10,height(thisTable)-reach)];

        % categoryAngles = [categoryAngles;thisTable(max([1,reach-39]):end_temp,:)];
        categoryAngles = [categoryAngles;thisTable((event_time(i)-39):(event_time(i)+80),:)];
    end
    % allAngles = [allAngles;categoryAngles];
    %Normalizing angles
    [zAngles,mean_raw,std_raw] = zscore(table2array(categoryAngles));

    %% PCA uncomment after variance
    [coeff,score,latent,tsquared,explained,mu] = pca(zAngles, 'Rows', 'pairwise');

    % Keeping tracking of the start and end times of each activity after
    % concatenation
    % if ~isfolder(fullfile(dest,'PCA_times\'))
    %     mkdir(fullfile(dest,'PCA_times\'));
    % end
    % writetable(table(activities,start_time,reach_time,end_time,grasp_time), fullfile(dest,'PCA_times\',strcat(anglefolder{i},'.csv')))

    % The score matrix contains the projected data
    if ~isfolder(fullfile(dest,'Projection\'))
        mkdir(fullfile(dest,'Projection\'));
    end
    writetable(array2table(score), fullfile(dest,'Projection\',strcat('Projection','.csv')))
    % The mu matrix contains the estimated mean of each variable, needed
    % for reprojection
    if ~isfolder(fullfile(dest,'PCA_mean\'))
        mkdir(fullfile(dest,'PCA_mean\'));
    end
    writetable(array2table(mu), fullfile(dest,'PCA_mean\',strcat('mean','.csv')))
    % Saving the raw data
    if ~isfolder(fullfile(dest,'Unprojected\'))
        mkdir(fullfile(dest,'Unprojected\'));
    end
    writetable(categoryAngles, fullfile(dest,'Unprojected\',strcat('Unprojected','.csv')))
    % Latent
    if ~isfolder(fullfile(dest,'Latent\'))
        mkdir(fullfile(dest,'Latent\'));
    end
    writetable(array2table(latent), fullfile(dest,'Latent\',strcat('Latent','.csv')))
    % The PCA coefficients
    if ~isfolder(fullfile(dest,'PCA_coeffs\'))
        mkdir(fullfile(dest,'PCA_coeffs\'));
    end
    writetable(array2table(coeff,'VariableNames',string(1:length(coeff))), fullfile(dest,'PCA_coeffs\',strcat('coeffs','.csv')))

    figure
    bar(explained, 'w');
    ylabel('Percentage of variance explained')
    xlabel('Principal Component')
    titlestr = strcat("PCA of grasp type: ", '');
    title(titlestr, 'Interpreter', 'none')
    set(gca,'fontsize',14, 'TickDir', 'out')
    ylim([0 80]);
    box off
    if ~isfolder(fullfile(dest,'PCA_explained\'))
        mkdir(fullfile(dest,'PCA_explained\'));
    end
    saveas(gcf,fullfile(dest,'PCA_explained\',strcat('ex','_PCA.png')))
    writetable(array2table(explained), fullfile(dest,'PCA_explained\',strcat('ex','.csv')))
    
    % [~,mean_raw,std_raw] = zscore(table2array(allAngles));
    if ~isfolder(fullfile(dest,'All\'))
        mkdir(fullfile(dest,'All\'));
    end
    writetable(array2table(mean_raw), fullfile(dest,'All\mu_global.csv'))
    % Standard deviation to reverse the normalization later
    if ~isfolder(fullfile(dest,'All\'))
        mkdir(fullfile(dest,'All\'));
    end
    writetable(array2table(std_raw), fullfile(dest,'All\sigma.csv'))
end

function PCA_participant(src,src_3d,dest)
    %% Finding PCs from participants.
    % src - folder location of all the angle data
    % src_3d - folder location of all the positinal data
    % dest - output file locations
    allAngles = [];
    maintemp = dir(fullfile(src,'*'));
    mainfolder = setdiff({maintemp([maintemp.isdir]).name},{'.','..'});
    % Looping through participants
    for i = 1:numel(mainfolder)
        categoryAngles = [];
        activities = [];
        start_time = 1;
        reach_time = [];
        end_time = [];
        grasp_time = [];
        acttemp = dir(fullfile(src,mainfolder{i},'*'));
        actfolder = setdiff({acttemp([acttemp.isdir]).name},{'.','..'});
        % Looping through activities
        for l = 1:numel(actfolder)
            subtemp = dir(fullfile(src,mainfolder{i},actfolder{l},'*.csv'));
            subfolder = {subtemp(~[subtemp.isdir]).name};
            postemp = dir(fullfile(src_3d,mainfolder{i},actfolder{l},'*.csv'));
            posfolder = {postemp(~[postemp.isdir]).name};
            
            % Concatonating angles from multiple activities
            for j = 1:numel(subfolder)
                file = fullfile(src,mainfolder{i},actfolder{l},subfolder{j});
                fprintf(1, 'Now reading %s\n', subfolder{j});
                thisTable = readtable(file);
        
                [filepath,name,ext] = fileparts(file);
                activities = [activities;convertCharsToStrings(name(1:end-7))];
        
                % Getting moment of grasp
                grasp_3d = fullfile(src_3d,mainfolder{i},actfolder{l},posfolder{j});
                grasp_3d = table2array(readtable(grasp_3d));
                vel = grasp_3d(1:end-1,:)-grasp_3d(2:end,:);
                vel = sqrt((vel(:,1:3:end)).^2 + (vel(:,2:3:end)).^2 + (vel(:,3:3:end)).^2);
                [~, vel_rise] = max(vel(1:end-1, 25) > 10 & vel(2:end, 25) > 10);
                [~, vel_fall] = max(vel(vel_rise:end-1, 25) < 10 & vel(vel_rise+1:end, 25) < 10);
                reach = vel_fall + vel_rise + 6;
                grasp_time = [grasp_time;reach];
        
                % Finding start and end times of each activity
                if ~((l == 1) && (j == 1))
                    start_time = [start_time;end_time(end)+1];
                end
                reach_time = [reach_time;start_time(end)+min([reach,40])-1];
                end_temp = min(reach+10,height(thisTable));
                end_time = [end_time;reach_time(end)+min(10,height(thisTable)-reach)];
        
                categoryAngles = [categoryAngles;thisTable(max([1,reach-39]):end_temp,:)];
            end
        end
        allAngles = [allAngles;categoryAngles];
        %Normalizing angles
        [zAngles,mean_raw,std_raw] = zscore(table2array(categoryAngles));
    
        %% PCA uncomment after variance
        [coeff,score,latent,tsquared,explained,mu] = pca(zAngles, 'Rows', 'pairwise');

        % Keeping tracking of the start and end times of each activity after
        % concatenation
        if ~isfolder(fullfile(dest,'PCA_times\'))
            mkdir(fullfile(dest,'PCA_times\'));
        end
        writetable(table(activities,start_time,reach_time,end_time,grasp_time), fullfile(dest,'PCA_times\',strcat(mainfolder{i},'.csv')))
    
        % The score matrix contains the projected data
        if ~isfolder(fullfile(dest,'Projection\'))
            mkdir(fullfile(dest,'Projection\'));
        end
        writetable(array2table(score), fullfile(dest,'Projection\',strcat(mainfolder{i},'.csv')))
        % The mu matrix contains the estimated mean of each variable, needed
        % for reprojection
        if ~isfolder(fullfile(dest,'PCA_mean\'))
            mkdir(fullfile(dest,'PCA_mean\'));
        end
        writetable(array2table(mean_raw), fullfile(dest,'PCA_mean\',strcat(mainfolder{i},'.csv')))
        % Standard deviation to reverse the normalization later
        if ~isfolder(fullfile(dest,'sigma\'))
            mkdir(fullfile(dest,'sigma\'));
        end
        writetable(array2table(std_raw), fullfile(dest,'sigma\',strcat(mainfolder{i},'.csv')))
        % Saving the raw data
        if ~isfolder(fullfile(dest,'Unprojected\'))
            mkdir(fullfile(dest,'Unprojected\'));
        end
        writetable(categoryAngles, fullfile(dest,'Unprojected\',strcat(mainfolder{i},'.csv')))
        % Latent
        if ~isfolder(fullfile(dest,'Latent\'))
            mkdir(fullfile(dest,'Latent\'));
        end
        writetable(array2table(latent), fullfile(dest,'Latent\',strcat(mainfolder{i},'.csv')))
        % The PCA coefficients
        if ~isfolder(fullfile(dest,'PCA_coeffs\'))
            mkdir(fullfile(dest,'PCA_coeffs\'));
        end
        writetable(array2table(coeff,'VariableNames',string(1:length(coeff))), fullfile(dest,'PCA_coeffs\',strcat(mainfolder{i},'.csv')))
    
        figure
        bar(explained, 'w');
        ylabel('Percentage of variance explained')
        xlabel('Principal Component')
        titlestr = strcat("PCA of grasp type: ", mainfolder{i});
        title(titlestr, 'Interpreter', 'none')
        set(gca,'fontsize',14, 'TickDir', 'out')
        ylim([0 80]);
        box off
        if ~isfolder(fullfile(dest,'PCA_explained\'))
            mkdir(fullfile(dest,'PCA_explained\'));
        end
        saveas(gcf,fullfile(dest,'PCA_explained\',strcat(mainfolder{i},'_PCA.png')))
        writetable(array2table(explained), fullfile(dest,'PCA_explained\',strcat(mainfolder{i},'.csv')))
    
    end
    [~,mean_raw,std_raw] = zscore(table2array(allAngles));
    if ~isfolder(fullfile(dest,'All\'))
        mkdir(fullfile(dest,'All\'));
    end
    writetable(array2table(mean_raw), fullfile(dest,'All\mu_global.csv'))
    % Standard deviation to reverse the normalization later
    if ~isfolder(fullfile(dest,'All\'))
        mkdir(fullfile(dest,'All\'));
    end
    writetable(array2table(std_raw), fullfile(dest,'All\sigma.csv'))
end

function PCA_act(src,src_3d,dest)
    %% Finding PCs from activities.
    % src - folder location of all the angle data
    % src_3d - folder location of all the positinal data
    % dest - output file locations
    
    maintemp = dir(fullfile(src,'*'));
    mainfolder = setdiff({maintemp([maintemp.isdir]).name},{'.','..'});
    % Looping through participants
    
    acttemp = dir(fullfile(src,mainfolder{1},'*'));
    actfolder = setdiff({acttemp([acttemp.isdir]).name},{'.','..'});
    % Looping through activities
    for l = 1:numel(actfolder)
        activities = [];
        start_time = [];
        reach_time = [];
        end_time = [];
        grasp_time = [];
        categoryAngles = [];
        for i = 1:numel(mainfolder)
            subtemp = dir(fullfile(src,mainfolder{i},actfolder{l},'*.csv'));
            subfolder = {subtemp(~[subtemp.isdir]).name};
            postemp = dir(fullfile(src_3d,mainfolder{i},actfolder{l},'*.csv'));
            posfolder = {postemp(~[postemp.isdir]).name};
            
            % Concatonating angles from multiple activities
            for j = 1:numel(subfolder)
                file = fullfile(src,mainfolder{i},actfolder{l},subfolder{j});
                fprintf(1, 'Now reading %s\n', subfolder{j});
                thisTable = readtable(file);
        
                [filepath,name,ext] = fileparts(file);
                activities = [activities;convertCharsToStrings(name(1:end-7))];
        
                % Getting moment of grasp
                grasp_3d = fullfile(src_3d,mainfolder{i},actfolder{l},posfolder{j});
                grasp_3d = table2array(readtable(grasp_3d));
                vel = grasp_3d(1:end-1,:)-grasp_3d(2:end,:);
                vel = sqrt((vel(:,1:3:end)).^2 + (vel(:,2:3:end)).^2 + (vel(:,3:3:end)).^2);
                [~, vel_rise] = max(vel(1:end-1, 25) > 10 & vel(2:end, 25) > 10);
                [~, vel_fall] = max(vel(vel_rise:end-1, 25) < 10 & vel(vel_rise+1:end, 25) < 10);
                reach = vel_fall + vel_rise + 6;
                grasp_time = [grasp_time;reach];
        
                % Finding start and end times of each activity
                if isempty(start_time)
                    start_time = 1;
                else
                    start_time = [start_time;end_time(end)+1];
                end
                reach_time = [reach_time;start_time(end)+min([reach,40])-1];
                end_temp = min(reach+10,height(thisTable));
                end_time = [end_time;reach_time(end)+min(10,height(thisTable)-reach)];

                categoryAngles = [categoryAngles;thisTable(max([1,reach-39]):end_temp,:)];
            end
        end
        if ~isfolder(fullfile(dest,'PCA_times\'))
            mkdir(fullfile(dest,'PCA_times\'));
        end
        writetable(table(activities,start_time,reach_time,end_time,grasp_time), fullfile(dest,'PCA_times\',strcat(actfolder{l},'.csv')))

        %Normalizing angles
        [zAngles,mean_raw,std_raw] = zscore(table2array(categoryAngles));

        [coeff,score,latent,tsquared,explained,mu] = pca(zAngles, 'Rows', 'pairwise');
        % The score matrix contains the projected data
        if ~isfolder(fullfile(dest,'Projection\'))
            mkdir(fullfile(dest,'Projection\'));
        end
        writetable(array2table(score), fullfile(dest,'Projection\',strcat(actfolder{l},'.csv')))
        % The mu matrix contains the estimated mean of each variable, needed
        % for reprojection
        if ~isfolder(fullfile(dest,'PCA_mean\'))
            mkdir(fullfile(dest,'PCA_mean\'));
        end
        writetable(array2table(mean_raw), fullfile(dest,'PCA_mean\',strcat(actfolder{l},'.csv')))
        % Standard deviation to reverse the normalization later
        if ~isfolder(fullfile(dest,'sigma\'))
            mkdir(fullfile(dest,'sigma\'));
        end
        writetable(array2table(std_raw), fullfile(dest,'sigma\',strcat(actfolder{l},'.csv')))
        % Latent
        if ~isfolder(fullfile(dest,'Latent\'))
            mkdir(fullfile(dest,'Latent\'));
        end
        writetable(array2table(latent), fullfile(dest,'Latent\',strcat(actfolder{l},'.csv')))
        % Saving the raw data
        if ~isfolder(fullfile(dest,'Unprojected\'))
            mkdir(fullfile(dest,'Unprojected\'));
        end
        writetable(categoryAngles, fullfile(dest,'Unprojected\',strcat(actfolder{l},'.csv')))
        % The PCA coefficients
        if ~isfolder(fullfile(dest,'PCA_coeffs\'))
            mkdir(fullfile(dest,'PCA_coeffs\'));
        end
        writetable(array2table(coeff,'VariableNames',string(1:length(coeff))), fullfile(dest,'PCA_coeffs\',strcat(actfolder{l},'.csv')))
    
        figure
        bar(explained, 'w');
        ylabel('Percentage of variance explained')
        xlabel('Principal Component')
        titlestr = strcat("PCA of grasp type: ", actfolder{l});
        title(titlestr, 'Interpreter', 'none')
        set(gca,'fontsize',14, 'TickDir', 'out')
        ylim([0 80]);
        box off
        if ~isfolder(fullfile(dest,'PCA_explained\'))
            mkdir(fullfile(dest,'PCA_explained\'));
        end
        saveas(gcf,fullfile(dest,'PCA_explained\',strcat(actfolder{l},'_PCA.png')))
        writetable(array2table(explained), fullfile(dest,'PCA_explained\',strcat(actfolder{l},'.csv')))
    end
end

function leaves = getLeaves(Z, node, nLeaves)
    %% Needed for cosine_cluster function to colour the leaves correctly
    children = Z(node,1:2);
    leaves = [];
    for j = 1:2
        if children(j) <= nLeaves
            leaves(end+1) = children(j);
        else
            leaves = [leaves, getLeaves(Z, children(j)-nLeaves, nLeaves)];
        end
    end
end

function generalSynergies = cosine_cluster(sub_PC, task_PC)
    %% Cluster function
    % sub_PC - All subject PCs
    % task_PC - All task PCs.

    %% Step 1: Transpose
    X = task_PC'; % transpose to make synergies as rows
    Y = sub_PC'; % transpose to make synergies as rows

    %% Step 2: Cluster task synergies using hierarchical clustering
    cosDist = pdist(X, 'cosine');
    numClusters = size(Y,1);
    Z = linkage(cosDist, 'complete');  % hierarchical linkage method
    clusterIdx = cluster(Z, 'maxclust', numClusters);

    %% Step 3: Reducing task pcs
    medoids = zeros(numClusters, size(X,2));
    medoidIdx    = zeros(numClusters,1);

    for k = 1:numClusters
        members = find(clusterIdx == k);
        if numel(members) == 1
            % Only one member in cluster
            medoidIdx(k) = members;
            medoids(k,:) = X(members,:);
            continue;
        end

        % Compute pairwise distances *within normalized cluster*
        D_sub = pdist2(X(members,:), X(members,:), 'cosine');

        % Medoid = member with smallest total distance to others
        [~, bestIdx] = min(sum(D_sub, 2));
        medoidIdx(k) = members(bestIdx);

        % Store corresponding rows from both normalized and raw data
        medoids(k,:) = X(medoidIdx(k), :);
    end

    allPCs = [medoids;Y];
    
    %% Step 1: Compute pairwise cosine distance
    % cosine distance = 1 - cosine similarity
    cosDist = pdist(allPCs, 'cosine');  % transpose to make synergies as rows
    
    %% Step 2: Cluster synergies using hierarchical clustering
    Z = linkage(cosDist, 'complete');  % hierarchical linkage method
    
    % 3️⃣ Find the largest distance jump
    distances = Z(:,3);
    [~, idxMaxJump] = max(diff(distances));
    
    % Number of clusters = total merges - index of biggest jump + 1
    numClusters = size(Z,1) - idxMaxJump + 1;
    
    fprintf('Largest linkage distance jump at step %d → %d clusters\n', ...
            idxMaxJump, numClusters);
    
    % 4️⃣ Assign clusters based on that threshold
    clusterIdx = cluster(Z, 'maxclust', numClusters);
    
    % Get MATLAB's default dendrogram leaf order (perm)
    [~,~,perm] = dendrogram(Z, 0);
    close;   % we only needed perm, not the plot
    
    % Build reorderIdx so that:
    %   - cluster 1 leaves come first (in their original order),
    %   - then cluster 2 leaves, etc.
    reorderIdx = zeros(size(perm));
    pos = 1;
    for c = 1:numClusters
        % leaves that belong to cluster c, in the original dendrogram order
        idxThisCluster = perm(clusterIdx(perm) == c);
        nThis          = numel(idxThisCluster);
    
        reorderIdx(pos:pos+nThis-1) = idxThisCluster;
        pos = pos + nThis;
    end

    figure;
    [H,~,leafOrder] = dendrogram(Z, 0, 'Reorder', reorderIdx);
    xlabel('Observations');
    ylabel('Distance');
    title('Dendrogram with clusters ordered 1..K');

    clusterColors = [
        0.00 0.45 0.74;   % 1 — Strong Blue
        0.85 0.33 0.10;   % 2 — Strong Orange
        0.93 0.69 0.13;   % 3 — Warm Yellow
        0.49 0.18 0.56;   % 4 — Deep Purple
        0.47 0.67 0.19    % 5 — Earthy Green
    ];
    ax = gca;

    % Cluster of each leaf in the *plotted* order
    leafClusters = clusterIdx(leafOrder);
    nLeaves = length(leafClusters);

    % Map each original observation index → its leaf position in the rendered dendrogram
    obsToLeafPos = zeros(nLeaves,1);
    for k = 1:nLeaves
        obsToLeafPos( leafOrder(k) ) = k;
    end

    for i = 1:size(Z,1)
        % Get leaves under this branch: ORIGINAL observation indices
        leaves_original = getLeaves(Z, i, nLeaves);
    
        % Convert original indices → leaf positions
        leafPositions = obsToLeafPos(leaves_original);
    
        % Which clusters appear under this branch (correctly indexed!)
        cl = unique( leafClusters(leafPositions) );
    
        if numel(cl) == 1
            % Pure branch
            set(H(i), 'Color', clusterColors(cl,:));
        else
            % Mixed branch
            set(H(i), 'Color', [0.4 0.4 0.4]);
        end
    end

    % Plot cutoff line
    hold on;
    yline(distances(idxMaxJump), 'r--', 'LineWidth', 1.5);
    title(sprintf('Hierarchical Clustering (nClusters = %d)', numClusters));
    xlabel('Observations'); ylabel('Linkage Distance');


    %% Step 4: Compute cluster centroids (general synergies)
    generalSynergies = zeros(numClusters, size(allPCs,2));
    
    for k = 1:numClusters
        members = (clusterIdx == k);
        if sum(members) == 0
            continue;
        end
        % Average the original (unnormalized, unaligned) PCs in this cluster
        Ck = mean(allPCs(members,:), 1);
        Ck_std = std(allPCs(members,:), 1)./sqrt(sum(members));
        generalSynergies(k,:) = Ck;
    end
end

function [coeff_general, latent_general, explained_general] = general_2PCA(src_par, src_act, src_la_par, src_la_act)
    %% Function to find general synergies using two stage PCA
    % src_par - File location of the participant PCA_coeffs
    % src_act = File location of the activity PCA_coeffs
    % src_la_par - File location of the participant Latents
    % src_la_act = File location of the activity Latents

    % Stack the first k PCs from each participant
    sub_PC = [];
    task_PC = [];
    latemp = dir(fullfile(src_la_par,'*.csv'));
    lafolder = {latemp(~[latemp.isdir]).name};
    maintemp = dir(fullfile(src_par,'*.csv'));
    mainfolder = {maintemp(~[maintemp.isdir]).name};
    for i = 1:numel(mainfolder)
        la_file = fullfile(src_la_par,lafolder{i});
        laArray = table2array(readtable(la_file));
        laArray = laArray > 1;
        k_new = sum(laArray);
        pc_file = fullfile(src_par,mainfolder{i});
        thisTable = readtable(pc_file);
        sub_PC = [sub_PC,table2array(thisTable(2:end,1:k_new))];
    end
    latemp2 = dir(fullfile(src_la_act,'*.csv'));
    lafolder2 = {latemp2(~[latemp2.isdir]).name};
    maintemp2 = dir(fullfile(src_act,'*.csv'));
    mainfolder2 = {maintemp2(~[maintemp2.isdir]).name};
    for i = 1:numel(mainfolder2)
        la_file = fullfile(src_la_act,lafolder2{i});
        laArray = table2array(readtable(la_file));
        laArray = laArray > 1;
        k_new = sum(laArray);
        pc_file = fullfile(src_act,mainfolder2{i});
        thisTable = readtable(pc_file);
        task_PC = [task_PC,table2array(thisTable(2:end,1:k_new))];
    end

     %% Step 1: Transpose
    X = task_PC'; % transpose to make synergies as rows
    Y = sub_PC'; % transpose to make synergies as rows

    %% Step 2: Cluster task synergies using hierarchical clustering
    cosDist = pdist(X, 'cosine');
    numClusters = size(Y,1);
    Z = linkage(cosDist, 'complete');  % hierarchical linkage method
    clusterIdx = cluster(Z, 'maxclust', numClusters);

    %% Step 3: Reducing task pcs
    medoids = zeros(numClusters, size(X,2));
    medoidIdx    = zeros(numClusters,1);

    for k = 1:numClusters
        members = find(clusterIdx == k);
        if numel(members) == 1
            % Only one member in cluster
            medoidIdx(k) = members;
            medoids(k,:) = X(members,:);
            continue;
        end

        % Compute pairwise distances *within normalized cluster*
        D_sub = pdist2(X(members,:), X(members,:), 'cosine');

        % Medoid = member with smallest total distance to others
        [~, bestIdx] = min(sum(D_sub, 2));
        medoidIdx(k) = members(bestIdx);

        % Store corresponding rows from both normalized and raw data
        medoids(k,:) = X(medoidIdx(k), :);
    end

    allPCs = [medoids;Y];

    [coeff_general, ~, latent_general , ~, explained_general, ~] = pca(allPCs, 'Rows', 'pairwise');
    latent_sig = latent_general > 1;
    k_new = sum(latent_sig);
    coeff_general = coeff_general(:,1:k_new);
end

function coeff_general = general_PCA_cluster(src_par, src_act, src_la_par, src_la_act)
    %% Function to find general synergies using hierarchical clustering
    % src_par - File location of the participant PCA_coeffs
    % src_act = File location of the activity PCA_coeffs
    % src_la_par - File location of the participant Latents
    % src_la_act = File location of the activity Latents

    % Stack the first k PCs from each participant
    PCs_subject = [];
    PCs_task = [];
    latemp = dir(fullfile(src_la_par,'*.csv'));
    lafolder = {latemp(~[latemp.isdir]).name};
    maintemp = dir(fullfile(src_par,'*.csv'));
    mainfolder = {maintemp(~[maintemp.isdir]).name};
    for i = 1:numel(mainfolder)
        la_file = fullfile(src_la_par,lafolder{i});
        laArray = table2array(readtable(la_file));
        laArray = laArray > 1;
        k_new = sum(laArray);
        pc_file = fullfile(src_par,mainfolder{i});
        thisTable = readtable(pc_file);
        PCs_subject = [PCs_subject,table2array(thisTable(2:end,1:k_new))];
    end
    latemp2 = dir(fullfile(src_la_act,'*.csv'));
    lafolder2 = {latemp2(~[latemp2.isdir]).name};
    maintemp2 = dir(fullfile(src_act,'*.csv'));
    mainfolder2 = {maintemp2(~[maintemp2.isdir]).name};
    for i = 1:numel(mainfolder2)
        la_file = fullfile(src_la_act,lafolder2{i});
        laArray = table2array(readtable(la_file));
        laArray = laArray > 1;
        k_new = sum(laArray);
        pc_file = fullfile(src_act,mainfolder2{i});
        thisTable = readtable(pc_file);
        PCs_task = [PCs_task,table2array(thisTable(2:end,1:k_new))];
    end
    coeff_general = cosine_cluster2(PCs_subject,PCs_task);
    coeff_general = coeff_general';
end

function output = angle2pos(reproarr,file,time,shoulder_focus)
    %% Converting joint angles into 3d positions
    % reproarr - PCs that have been reprojected back into angle space
    % file - 3d position file location to get limb lengths and palm
    % starting position
    % time - time stamp within the 3d joint position file to calculate the
    % limb lengths
    % shoulder_focus - whether to include the arm within the
    % reconstruction. 1 for yes, 0 for no.

    grasp_time = 1;
    grasp_time_3d = time;
    T = readtable(file);
    original = table2array(T);
    original(:,61:63) = [];
    lengths = pdist([original(grasp_time_3d,1),original(grasp_time_3d,2),original(grasp_time_3d,3);original(grasp_time_3d,70),original(grasp_time_3d,71),original(grasp_time_3d,72)]);
    % Finding the vectors to describe [TCMC, IMCP, MMCP, RMCP, LMCP]
    knuckles = [original(grasp_time_3d,1),original(grasp_time_3d,2),original(grasp_time_3d,3)]-[original(grasp_time_3d,70),original(grasp_time_3d,71),original(grasp_time_3d,72)];
    % Realigning original data
    for k = 0:22
        if (mod(k+1,4) == 0) && (k < 19)
            limb = [original(grasp_time_3d,(k+1)*3+1),original(grasp_time_3d,(k+1)*3+2),original(grasp_time_3d,(k+1)*3+3);original(grasp_time_3d,70),original(grasp_time_3d,71),original(grasp_time_3d,72)];
            knuckles = cat(1,knuckles,[original(grasp_time_3d,(k+1)*3+1),original(grasp_time_3d,(k+1)*3+2),original(grasp_time_3d,(k+1)*3+3)]-[original(grasp_time_3d,70),original(grasp_time_3d,71),original(grasp_time_3d,72)]);
        elseif (mod(k+1,4) == 0) && (k == 19)
                continue
        else
            limb = [original(grasp_time_3d,k*3+1),original(grasp_time_3d,k*3+2),original(grasp_time_3d,k*3+3);original(grasp_time_3d,k*3+4),original(grasp_time_3d,k*3+5),original(grasp_time_3d,k*3+6)];
        end
        lengths = cat(1,lengths,pdist(limb));
    end
    
    %% Adjusting angles
    palm_plane = -cross(knuckles(2,:),knuckles(5,:));
    z_axis = [0 0 1];
    % Determine the angle between the vector and the z-axis
    theta = acos(dot(palm_plane, z_axis)/(norm(palm_plane)*norm(z_axis)));
    % Determine the axis of rotation
    axis = cross(palm_plane, z_axis)/norm(cross(palm_plane, z_axis));
    % Construct the rotation matrix using Rodrigues' formula
    K = [0 -axis(3) axis(2); axis(3) 0 -axis(1); -axis(2) axis(1) 0];
    R1 = eye(3) + sin(theta)*K + (1-cos(theta))*K*K;
    % We are rotating the palm_plane by angle theta around vector K
    knuckles_t = transpose(R1*knuckles');
    
    y_axis = [0 1 0];
    % Determine the angle between the vector and the z-axis
    theta = acos(dot(knuckles_t(2,:), y_axis)/(norm(knuckles_t(2,:))*norm(y_axis)));
    % Determine the axis of rotation
    axis = cross(knuckles_t(2,:), y_axis)/norm(cross(knuckles_t(2,:), y_axis));
    % Construct the rotation matrix using Rodrigues' formula
    K = [0 -axis(3) axis(2); axis(3) 0 -axis(1); -axis(2) axis(1) 0];
    R2 = eye(3) + sin(theta)*K + (1-cos(theta))*K*K;
    
    original_R = R2*R1;
    % We are rotating the palm_plane by angle theta around vector K
    palm_plane = (original_R*palm_plane')';
    palm_plane = palm_plane./norm(palm_plane);
    knuckles = transpose(original_R*knuckles');
    
    PIP = [];
    DIP = [];
    Tip = [];
    
    for k = 2:5
        m = (k-1)*4+2;
        v2_np = palm_plane.*((knuckles(k,1).*palm_plane(1) + knuckles(k,2).*palm_plane(2) + knuckles(k,3).*palm_plane(3)));   % Vector component normal to the palm plane
        v2_p = [knuckles(k,1)-v2_np(1),knuckles(k,2)-v2_np(2),knuckles(k,3)-v2_np(3)];                                       % Vector component on the palm plane
        v2mag = sqrt(v2_p(1).^2 + v2_p(2).^2 + v2_p(3).^2);
        v2norm = v2_p./v2mag;
        
        x_axis = [1 0 0];
        z_axis = [0 0 1];
        % Determine the angle between the vector and the x-axis
        theta = acos(dot(x_axis, v2norm)/(norm(x_axis)*norm(v2norm)));
        % Determine the axis of rotation
        axis = cross(x_axis, v2norm)/norm(cross(x_axis, v2norm));
        % Construct the rotation matrix using Rodrigues' formula
        K = [0 -axis(3) axis(2); axis(3) 0 -axis(1); -axis(2) axis(1) 0];
        R1 = eye(3) + sin(theta)*K + (1-cos(theta))*K*K;
        new_z_axis = (R1*z_axis')';
        
        % Determine the angle between the vector and the y-axis
        theta = acos(dot(new_z_axis, palm_plane)/(norm(new_z_axis)*norm(palm_plane)));
        % Determine the axis of rotation
        axis = cross(new_z_axis, palm_plane)/norm(cross(new_z_axis, palm_plane));
        % Construct the rotation matrix using Rodrigues' formula
        K = [0 -axis(3) axis(2); axis(3) 0 -axis(1); -axis(2) axis(1) 0];
        R2 = eye(3) + sin(theta)*K + (1-cos(theta))*K*K;
        R_knuckle = R2*R1;
        [x,y,z] = sph2cart(reproarr(grasp_time,k+16),reproarr(grasp_time,(k-1)*3+1),lengths(m));
        PIP_temp = [x,y,z];
        PIP_temp = (R_knuckle*PIP_temp')';
    
        finger_axis_np = palm_plane.*(PIP_temp(1).*palm_plane(1) + PIP_temp(2).*palm_plane(2) + PIP_temp(3).*palm_plane(3));  % Vector component normal to the thumb plane
        finger_axis_p = [PIP_temp(1)-finger_axis_np(1),PIP_temp(2)-finger_axis_np(2),PIP_temp(3)-finger_axis_np(3)];                                    % Vector component on the thumb plane
        finger_axismag = sqrt(finger_axis_p(1).^2 + finger_axis_p(2).^2 + finger_axis_p(3).^2);
        finger_axis = finger_axis_p./finger_axismag;
        K = [0 -palm_plane(3) palm_plane(2); palm_plane(3) 0 -palm_plane(1); -palm_plane(2) palm_plane(1) 0];
        R = eye(3) + sin(-pi/2)*K + (1-cos(-pi/2))*K*K;
        finger_axis = (R*finger_axis')';
        PIP = [PIP;PIP_temp + knuckles(k,:)];
        K = [0 -finger_axis(3) finger_axis(2); finger_axis(3) 0 -finger_axis(1); -finger_axis(2) finger_axis(1) 0];
        R = eye(3) + sin(reproarr(grasp_time,(k-1)*3+2))*K + (1-cos(reproarr(grasp_time,(k-1)*3+2)))*K*K;
        DIP = [DIP;(R*(PIP(k-1,:)-knuckles(k,:))')'.*(lengths(m+1)/lengths(m)) + PIP(k-1,:)];
        R = eye(3) + sin(reproarr(grasp_time,(k-1)*3+3))*K + (1-cos(reproarr(grasp_time,(k-1)*3+3)))*K*K;
        Tip = [Tip;(R*(DIP(k-1,:)-PIP(k-1,:))')'.*(lengths(m+2)/lengths(m+1)) + DIP(k-1,:)];
    end
    
    %% Fix definition
    thumb_rot_axis = knuckles(2,:) - knuckles(1,:);
    thumb_rot_axis = thumb_rot_axis./norm(thumb_rot_axis);
    thumb_plane = cross(knuckles(1,:),knuckles(2,:));
    thumb_plane = thumb_plane./norm(thumb_plane);
    thumb_abd_axis = cross(thumb_plane,thumb_rot_axis);
    thumb_abd_axis = (thumb_abd_axis./norm(thumb_abd_axis));
    
    x_axis = [1 0 0];
    y_axis = [0 1 0];
    % Determine the angle between the vector and the z-axis
    theta = acos(dot(x_axis, thumb_rot_axis)/(norm(x_axis)*norm(thumb_rot_axis)));
    % Determine the axis of rotation
    axis = cross(x_axis, thumb_rot_axis)/norm(cross(x_axis, thumb_rot_axis));
    % Construct the rotation matrix using Rodrigues' formula
    K = [0 -axis(3) axis(2); axis(3) 0 -axis(1); -axis(2) axis(1) 0];
    R1 = eye(3) + sin(theta)*K + (1-cos(theta))*K*K;
    % We are rotating the palm_plane by angle theta around vector K
    new_y_axis = (R1*y_axis')';
    
    % Determine the angle between the vector and the z-axis
    theta = acos(dot(new_y_axis, thumb_abd_axis)/(norm(new_y_axis)*norm(thumb_abd_axis)));
    % Determine the axis of rotation
    axis = cross(new_y_axis, thumb_abd_axis)/norm(cross(new_y_axis, thumb_abd_axis));
    % Construct the rotation matrix using Rodrigues' formula
    K = [0 -axis(3) axis(2); axis(3) 0 -axis(1); -axis(2) axis(1) 0];
    R2 = eye(3) + sin(theta)*K + (1-cos(theta))*K*K;
    % We are rotating the palm_plane by angle theta around vector K
    R_thumb = R2*R1;
    
    [x,y,z] = sph2cart(reproarr(grasp_time,1),reproarr(grasp_time,16),lengths(2));
    TMCP = [x,y,z];
    TMCP = (R_thumb*TMCP')';
    TMCP_mag = sqrt(TMCP(:,1).^2 + TMCP(:,2).^2 + TMCP(:,3).^2);
    
    thumb_plane_np = TMCP.*((thumb_plane(:,1).*TMCP(:,1) + thumb_plane(:,2).*TMCP(:,2) + thumb_plane(:,3).*TMCP(:,3))./TMCP_mag.^2);  % thumb plane component along v1
    thumb_plane_p = thumb_plane - thumb_plane_np;                               % thumb plane component orthogonal to v1
    thumb_plane2 = thumb_plane_p./norm(thumb_plane_p);
    TMCP_norm = TMCP./TMCP_mag;
    K = [0 -TMCP_norm(3) TMCP_norm(2); TMCP_norm(3) 0 -TMCP_norm(1); -TMCP_norm(2) TMCP_norm(1) 0];
    R = eye(3) + sin(reproarr(grasp_time,17))*K + (1-cos(reproarr(grasp_time,17)))*K*K;
    thumb_plane2 = (R*thumb_plane2')';
    thumb_plane2 = -thumb_plane2./norm(thumb_plane2);
    TMCP = TMCP + knuckles(1,:);
    
    K = [0 -thumb_plane2(3) thumb_plane2(2); thumb_plane2(3) 0 -thumb_plane2(1); -thumb_plane2(2) thumb_plane2(1) 0];
    R = eye(3) + sin(reproarr(grasp_time,2))*K + (1-cos(reproarr(grasp_time,2)))*K*K;
    TIP = (R*(TMCP-knuckles(1,:))')'.*(lengths(3)/lengths(2)) + TMCP;
    R = eye(3) + sin(reproarr(grasp_time,3))*K + (1-cos(reproarr(grasp_time,3)))*K*K;
    TT = (R*(TIP-TMCP)')'.*(lengths(4)/lengths(3)) + TIP;
    
    % Elbow starts in the direction of MMCP to W
    x_axis = [1 0 0];
    z_axis = [0 0 1];
    wrist_forward_np = palm_plane.*((knuckles(3,1).*palm_plane(1) + knuckles(3,2).*palm_plane(2) + knuckles(3,3).*palm_plane(3)));   % Vector component normal to the palm plane
    wrist_forward_p = [knuckles(3,1)-wrist_forward_np(1),knuckles(3,2)-wrist_forward_np(2),knuckles(3,3)-wrist_forward_np(3)];                                       % Vector component on the palm plane
    wrist_forwardmag = sqrt(wrist_forward_p(1).^2 + wrist_forward_p(2).^2 + wrist_forward_p(3).^2);
    wrist_forwardnorm = wrist_forward_p./wrist_forwardmag;
    wrist_fe_axis = cross(palm_plane, wrist_forwardnorm,2);
    % Determine the angle between the vector and the x-axis
    theta = acos(dot(x_axis, wrist_fe_axis)/(norm(x_axis)*norm(wrist_fe_axis)));
    % Determine the axis of rotation
    axis = cross(x_axis, wrist_fe_axis)/norm(cross(x_axis, wrist_fe_axis));
    % Construct the rotation matrix using Rodrigues' formula
    K = [0 -axis(3) axis(2); axis(3) 0 -axis(1); -axis(2) axis(1) 0];
    R1 = eye(3) + sin(theta)*K + (1-cos(theta))*K*K;
    new_z_axis = (R1*z_axis')';
    
    % Determine the angle between the vector and the y-axis
    theta = acos(dot(new_z_axis, palm_plane)/(norm(new_z_axis)*norm(palm_plane)));
    % Determine the axis of rotation
    axis = cross(new_z_axis, palm_plane)/norm(cross(new_z_axis, palm_plane));
    % Construct the rotation matrix using Rodrigues' formula
    K = [0 -axis(3) axis(2); axis(3) 0 -axis(1); -axis(2) axis(1) 0];
    R2 = eye(3) + sin(theta)*K + (1-cos(theta))*K*K;
    R_Elbow = R2*R1;
    [x,y,z] = sph2cart(reproarr(grasp_time,22),reproarr(grasp_time,23),lengths(23));
    E = [x,y,z];
    E = (R_Elbow*E')';
    
    v1_np = palm_plane.*(knuckles(3,1).*palm_plane(:,1) + knuckles(3,2).*palm_plane(:,2) + knuckles(3,3).*palm_plane(:,3));   % Vector component normal to the palm plane
    v1_p = [knuckles(3,1)-v1_np(:,1),knuckles(3,2)-v1_np(:,2),knuckles(3,3)-v1_np(:,3)];                                       % Vector component on the palm plane
    v1mag = sqrt(v1_p(:,1).^2 + v1_p(:,2).^2 + v1_p(:,3).^2);
    wrist_forward = v1_p./v1mag;
    wrist_fe_axis = cross(palm_plane,wrist_forward,2);  % Wrist flexion axis
    wrist_fe_axis_mag = sqrt(wrist_fe_axis(:,1).^2 + wrist_fe_axis(:,2).^2 + wrist_fe_axis(:,3).^2);
    wrist_fe_axis = wrist_fe_axis./wrist_fe_axis_mag;
    forearm_np = E.*((wrist_fe_axis(:,1).*E(:,1) + wrist_fe_axis(:,2).*E(:,2) + wrist_fe_axis(:,3).*E(:,3))./(lengths(23).^2));   % Vector component normal to the palm plane
    forearm_p = [wrist_fe_axis(:,1)-forearm_np(:,1),wrist_fe_axis(:,2)-forearm_np(:,2),wrist_fe_axis(:,3)-forearm_np(:,3)];        % Vector component on the palm plane
    forearmmag = sqrt(forearm_p(:,1).^2 + forearm_p(:,2).^2 + forearm_p(:,3).^2);
    forearm_ref = forearm_p./forearmmag;
    
    E_rot = E/lengths(23);
    K = [0 -E_rot(3) E_rot(2); E_rot(3) 0 -E_rot(1); -E_rot(2) E_rot(1) 0];
    R = eye(3) + sin(-reproarr(grasp_time,24))*K + (1-cos(-reproarr(grasp_time,24)))*K*K;
    E_axis = (R*forearm_ref')';
    K = [0 -E_axis(3) E_axis(2); E_axis(3) 0 -E_axis(1); -E_axis(2) E_axis(1) 0];
    R = eye(3) + sin(reproarr(grasp_time,25))*K + (1-cos(reproarr(grasp_time,25)))*K*K;
    S = (R*E')'.*lengths(22)./lengths(23) + E;
    
    s_third_axis = cross(E_axis,(S-E),2);
    x_axis = [1 0 0];
    z_axis = [0 0 1];
    % Determine the angle between the vector and the z-axis
    theta = acos(dot(x_axis, S-E)/(norm(x_axis)*norm(S-E)));
    % Determine the axis of rotation
    axis = cross(x_axis, S-E)/norm(cross(x_axis, S-E));
    % Construct the rotation matrix using Rodrigues' formula
    K = [0 -axis(3) axis(2); axis(3) 0 -axis(1); -axis(2) axis(1) 0];
    R1 = eye(3) + sin(theta)*K + (1-cos(theta))*K*K;
    % We are rotating the palm_plane by angle theta around vector K
    new_z_axis = (R1*z_axis')';
    
    % Determine the angle between the vector and the z-axis
    theta = acos(dot(new_z_axis, E_axis)/(norm(new_z_axis)*norm(E_axis)));
    % Determine the axis of rotation
    axis = cross(new_z_axis, E_axis)/norm(cross(new_z_axis, E_axis));
    % Construct the rotation matrix using Rodrigues' formula
    K = [0 -axis(3) axis(2); axis(3) 0 -axis(1); -axis(2) axis(1) 0];
    R2 = eye(3) + sin(theta)*K + (1-cos(theta))*K*K;
    % We are rotating the palm_plane by angle theta around vector K
    R_Shoulder = R2*R1;
    
    [x,y,z] = sph2cart(reproarr(grasp_time,26),reproarr(grasp_time,27),lengths(21));
    C = [x,y,z];
    C = (R_Shoulder*C')' + S;

    W = [0,0,0];

    %% Chest centred
    if shoulder_focus == 1
        shoulder = (C-S)./norm(C-S);
        x_axis = [-1 0 0];
        % Determine the angle between the vector and the x-axis
        theta = acos(dot(shoulder, x_axis)/(norm(shoulder)*norm(x_axis)));
        % Determine the axis of rotation
        axis = cross(shoulder, x_axis)/norm(cross(shoulder, x_axis));
        % Construct the rotation matrix using Rodrigues' formula
        K = [0 -axis(3) axis(2); axis(3) 0 -axis(1); -axis(2) axis(1) 0];
        R1 = eye(3) + sin(theta)*K + (1-cos(theta))*K*K;
        % We are rotating the arm_plane by angle theta around vector K
        % upper_arm = transpose(R1*norm(S-E)');
        % z_axis = [0 0 1];
        % % Determine the angle between the vector and the z-axis
        % theta = acos(dot(upper_arm(2,:), z_axis)/(norm(upper_arm(2,:))*norm(z_axis)));
        % % Determine the axis of rotation
        % axis = cross(upper_arm(2,:), z_axis)/norm(cross(upper_arm(2,:), z_axis));
        % % Construct the rotation matrix using Rodrigues' formula
        % K = [0 -axis(3) axis(2); axis(3) 0 -axis(1); -axis(2) axis(1) 0];
        % R2 = eye(3) + sin(theta)*K + (1-cos(theta))*K*K;
    
        original_R = R1;%*R2;
        % Chest centered plotting
        origin = transpose(original_R*S');
        knuckles = transpose(original_R*knuckles')-origin;
        TMCP = transpose(original_R*TMCP')-origin;
        TIP = transpose(original_R*TIP')-origin;
        TT = transpose(original_R*TT')-origin;
        PIP = transpose(original_R*PIP')-origin;
        DIP = transpose(original_R*DIP')-origin;
        Tip = transpose(original_R*Tip')-origin;
        E = transpose(original_R*E')-origin;
        S = transpose(original_R*S')-origin;
        C = transpose(original_R*C')-origin;
        W = W-origin;
    end

    output = [TMCP;TIP;TT;knuckles;PIP;DIP;Tip;C;S;E;W];
end

function draw_skeleton(reproarr, colour, time, file, s_focus)
    %% Drawing 3d position of upper limb from joint angles
    % reproarr - PCs that have been reprojected back into angle space
    % colour - the colour of the skeleton
    % time - time stamp within the 3d joint position file to calculate the
    % limb lengths
    % file - 3d position file location to get limb lengths and palm
    % starting position
    % s_focus - whether to include the arm within the reconstruction. 1 for
    % yes, 0 for no.

    reproarr = reproarr';
    %% Drawing reprojection
    % [TMCP;TIP;TT;TMCP;IMCP;MMCP;RMCP;LMCP;IPIP;MPIP;RPIP;LPIP;IDIP;MDIP;RDIP;LDIP;IT;MT;RT;LT;C;S;E;W];
    pos = angle2pos(reproarr,file,time,s_focus);
    hold on
    ThumbPlotre = plot3([pos(24,1),pos(4,1),pos(1,1),pos(2,1),pos(3,1)], ...
        [pos(24,2),pos(4,2),pos(1,2),pos(2,2),pos(3,2)], ...
        [pos(24,3),pos(4,3),pos(1,3),pos(2,3),pos(3,3)], ...
        '-o', 'MarkerSize',3,'MarkerFaceColor',	colour, 'Color',colour);
    IndexPlotre = plot3([pos(24,1),pos(5,1),pos(9,1),pos(13,1),pos(17,1)], ...
        [pos(24,2),pos(5,2),pos(9,2),pos(13,2),pos(17,2)], ...
        [pos(24,3),pos(5,3),pos(9,3),pos(13,3),pos(17,3)], ...
        '-o', 'MarkerSize',3,'MarkerFaceColor',	colour, 'Color',colour);
    MiddlePlotre = plot3([pos(24,1),pos(6,1),pos(10,1),pos(14,1),pos(18,1)], ...
        [pos(24,2),pos(6,2),pos(10,2),pos(14,2),pos(18,2)], ...
        [pos(24,3),pos(6,3),pos(10,3),pos(14,3),pos(18,3)], ...
        '-o', 'MarkerSize',3,'MarkerFaceColor',	colour, 'Color',colour);
    RingPlotre = plot3([pos(24,1),pos(7,1),pos(11,1),pos(15,1),pos(19,1)], ...
        [pos(24,2),pos(7,2),pos(11,2),pos(15,2),pos(19,2)], ...
        [pos(24,3),pos(7,3),pos(11,3),pos(15,3),pos(19,3)], ...
        '-o', 'MarkerSize',3,'MarkerFaceColor',	colour, 'Color',colour);
    LittlePlotre = plot3([pos(24,1),pos(8,1),pos(12,1),pos(16,1),pos(20,1)], ...
        [pos(24,2),pos(8,2),pos(12,2),pos(16,2),pos(20,2)], ...
        [pos(24,3),pos(8,3),pos(12,3),pos(16,3),pos(20,3)], ...
        '-o', 'MarkerSize',3,'MarkerFaceColor',	colour, 'Color',colour);
    if s_focus
        RightArmPlotre = plot3([pos(24,1),pos(23,1),pos(22,1),pos(21,1)], ...
            [pos(24,2),pos(23,2),pos(22,2),pos(21,2)], ...
            [pos(24,3),pos(23,3),pos(22,3),pos(21,3)], ...
            '-o', 'MarkerSize',3,'MarkerFaceColor',colour, 'Color',colour);
    end

    set(gca,'DataAspectRatio',[1 1 1])
    xlabel('x (mm)')
    ylabel('y (mm)')
    zlabel('z (mm)')
    titlestr = ["Mean Hand Shape","Hand Shape with Synergy","Minimum"];
    title(titlestr(1))
    view([240 20])
end

function general_synergy(co_gen, mu_gen, k, amp, colour, time, file)
    %% Drawing synergy
    % co_gen - matrix of synergies
    % mu_gen - the mean posture
    % k - position of the PC to visualize within the matrix of PCs
    % amp - amplitude of the PC
    % colour - the colour of the skeleton
    % time - time stamp within the 3d joint position file to calculate the
    % limb lengths
    % file - 3d position file location to get limb lengths and palm
    % starting position

    mean_posture = mu_gen;
    
    % Visualize postures
    posture_plus  = mean_posture + amp * co_gen(:,k);
    
    figure
    draw_skeleton(posture_plus,colour, time, file, 1);
    draw_skeleton(mean_posture,'k', time, file, 1);
    title(['Posture variation along PC ' num2str(k)]);

    joint_names = ["TCMC f","TMCP f","TIP f","IMCP f","IPIP f","IDIP f","MMCP f",...
    "MPIP f","MDIP f","RMCP f","RPIP f","RDIP f","LMCP f","LPIP f","LDIP f",...
    "TMCP a","TCMC r","IMCP a","MMCP a","RMCP a","LMCP a","W a","W f","W r",...
    "RE f","RS f","RS a","RS r"];
    figure
    co_gen(26,k) = co_gen(26,k)*-1; % Flipping shoulder around since it was defined incorrectly in the wrong direction.
    bar(joint_names,co_gen(:,k),'w');
    box off
    xlabel('Joint');
    ylabel('Loading');
    title(['PC ' num2str(k) ' Coefficients'])
    set(gca,'fontsize',14, 'TickDir', 'out')
end

function distance = reco_err(co_gen, mu_gen, sigma, k, time, idx, score_file, file_true)
    %% Calculating reconstruction error
    % co_gen - matrix of synergies
    % mu_gen - the mean posture
    % sigma - standard deviation of the mean posture
    % k - position of the PC to visualize within the matrix of PCs
    % time - time stamp within the 3d joint position file to calculate the
    % limb lengths
    % idx - moment of reconstruction within the score file
    % score_file - file of the projected data output from the initial PCA
    % file_true - file of the angles data to compare with the projected
    % data to find the error

    try
        score = table2array(readtable(score_file));
        score = score(idx,:);
        recon = (score(:,1:k) * co_gen(:,1:k)')';
        reproarr = recon .* repmat(sigma, 1, size(recon,2)) + repmat(mu_gen, 1, size(recon,2));
        reproarr = reproarr';
        angles_true = table2array(readtable(file_true));
        angles_true = angles_true(time,:);
        distance = angles_true - reproarr;
        %% Drawing reprojection
        % pos_guess = angle2pos(reproarr,file,time,1);
        % pos_true = angle2pos(angles_true,file,time,1);
        % error = pos_true-pos_guess;
        % distance = [];
        % for k = 1:size(error,1)
        %     distance = [distance,sqrt(error(k,1)^2+error(k,2)^2+error(k,3)^2)];
        % end
        % figure
        % hold on
        % draw_skeleton(reproarr','r', time, 0);
        % draw_skeleton(angles_true','k', time, 0);
    catch
        distance = zeros(1,28);
    end
end

function recon_err_all(co_gen, mu_gen, sigma, pro_src, time_src, angle_src)
    %% Reconstruction error across all participants
    % sigma - array containing the standard deviation of the mean posture data
    % mu_gen - array containing the mean posture data
    % pro_src - folder location containing all the projected data from
    % participants
    % time_src - folder location containing all time data from participants
    % angle_src - folder location containing all raw angle data

    protemp = dir(fullfile(pro_src,'*.csv'));
    profolder = {protemp(~[protemp.isdir]).name};
    timetemp = dir(fullfile(time_src,'*.csv'));
    timefolder = {timetemp(~[timetemp.isdir]).name};
    maintemp = dir(fullfile(angle_src,'*'));
    mainfolder = setdiff({maintemp([maintemp.isdir]).name},{'.','..'});
    total_dis = [];
    for i = 1:numel(mainfolder)
        %% For each participant
        timeTable = readtable(fullfile(time_src,timefolder{i}));
        acttemp = dir(fullfile(angle_src,mainfolder{i},'*'));
        actfolder = setdiff({acttemp([acttemp.isdir]).name},{'.','..'});
        for l = 1:numel(actfolder)
            %% For each activity
            subtemp = dir(fullfile(angle_src,mainfolder{i},actfolder{l},'*.csv'));
            subfolder = {subtemp(~[subtemp.isdir]).name};
            for j = 1:numel(subfolder)
                file = fullfile(angle_src,mainfolder{i},actfolder{l},subfolder{j});
                activity = subfolder{j};
                activity = activity(1:end-11);
                reach_time = timeTable.reach_time(strcmp(timeTable.activities, activity));
                grasp_time = timeTable.grasp_time(strcmp(timeTable.activities, activity));
                distance = reco_err(co_gen, mu_gen, sigma, 4, grasp_time, reach_time, fullfile(pro_src,profolder{i}), file);
                total_dis = [total_dis;distance];
            end
        end
    end
    %% Graphing error distance
    joint_names = ["TCMC_f","TMCP_f","TIP_f","IMCP_f","IPIP_f","IDIP_f","MMCP_f",...
        "MPIP_f","MDIP_f","RMCP_f","RPIP_f","RDIP_f","LMCP_f","LPIP_f","LDIP_f",...
        "TMCP_a","TCMC_r","IMCP_a","MMCP_a","RMCP_a","LMCP_a","W_a","W_f","W_r",...
        "RE_f","RS_f","RS_a","RS_r"];
    meanPC = mean(rad2deg(total_dis));
    meanPC(26) = meanPC(26)*-1;
    semPC = std(rad2deg(total_dis))./sqrt(size(total_dis,1));
    meanPC = flip(meanPC);
    semPC = flip(semPC);
    figure
    hold on
    barh(flip(joint_names),meanPC, 'w');
    errorbar(meanPC, 1:28, semPC, 'k.', 'horizontal', 'LineWidth', 1.2, 'CapSize', 10);
    ylabel("Joints")
    xlabel("Angle error (degrees)")
    title("Task Focused Joint Reconstruction Error")
    set(gca,'fontsize',14, 'TickDir', 'out')
end


%% Characterizing top 4 synergies
PCA_testing('C:\Users\czhe0008\Documents\EEG\3d\11_12_angles\','C:\Users\czhe0008\Documents\EEG\3d\11_12_filt\','C:\Users\czhe0008\Documents\EEG\PCA\11_12\', ALLEEG)
% [co_gen,count1,count2] = general_PCA_cluster('C:\Users\czhe0008\Documents\EEG\PCA\11_12\PCA_coeffs\', 'C:\Users\czhe0008\Documents\EEG\PCA\11_12\PCA_coeffs\', 'C:\Users\czhe0008\Documents\EEG\PCA\11_12\Latent\', 'C:\Users\czhe0008\Documents\EEG\PCA\11_12\Latent\'); % cluster arm
% mu_gen = table2array(readtable('file location of mu_global.csv'))';
% file = table2array(readtable('file location of any trial to use as the base'))';
% time = moment of time within the chosen file to use as the base
clusterColors = [
    0.00 0.45 0.74;   % 1 — Strong Blue
    0.85 0.33 0.10;   % 2 — Strong Orange
    0.93 0.69 0.13;   % 3 — Warm Yellow
    0.49 0.18 0.56;   % 4 — Deep Purple
    0.47 0.67 0.19    % 5 — Earthy Green
];
for i = 1:4
    general_synergy(co_gen, mu_gen, i, 1, clusterColors(i,:), time, file);
    % general_synergy(co_gen2, mu_gen, i, 1, clusterColors(i,:), time, file);
end
