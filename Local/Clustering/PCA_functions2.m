function generalSynergies = cosine_cluster(synergies)
%% Example: Clustering Hand Synergies Across Subjects
    % Assume each column of 'synergies' is a synergy vector
    % Rows = joints, Columns = different synergies from all subjects
    % Example: 15 joints, 12 synergies total (4 subjects × 3 synergies each)
    
    %% Step 1: Compute pairwise cosine distance
    % cosine distance = 1 - cosine similarity
    cosDist = pdist(synergies', 'cosine');  % transpose to make synergies as rows
    
    %% Step 2: Cluster synergies using hierarchical clustering
    Z = linkage(cosDist, 'average');  % hierarchical linkage method

    % Plot dendrogram to visualize clustering
    figure;
    dendrogram(Z,0,'ColorThreshold',0.9);
    title('Hierarchical Clustering of Synergies (Cosine Distance)');
    ylabel('Cosine Distance')
    set(gca,'fontsize',14, 'TickDir', 'out')
    
    %% Step 3: Assign clusters
    % numClusters = 4;  % choose number of general synergies you expect
    % clusterIdx = cluster(Z, 'maxclust', numClusters);
    clusterIdx = cluster(Z,"cutoff",0.9,'criterion','distance');
    numClusters = numel(unique(clusterIdx));

    % Display which synergy belongs to which cluster
    % disp('Synergy cluster assignments:');
    % disp(clusterIdx);
    
    %% Step 4: Compute cluster centroids (general synergies)
    generalSynergies = zeros(size(synergies,1), numClusters);
    % generalLatent = zeros(size(latent,1), numClusters);
    % generalMu = zeros(size(mu,1), numClusters);
    for k = 1:numClusters
        generalSynergies(:,k) = mean(synergies(:, clusterIdx == k), 2);
        % generalMu(:,k) = mean(mu(:, clusterIdx == k), 2);
        % generalLatent(k) = mean(latent(:, clusterIdx == k), 2);
    end
    
    % Optional: normalize general synergies to unit length
    % generalSynergies = generalSynergies ./ vecnorm(generalSynergies);
    
    % disp('General synergies (cluster centroids):');
    % disp(generalSynergies);
end

function generalSynergies = cosine_cluster2(sub_PC, task_PC)
%% Example: Clustering Hand Synergies Across Subjects
    % Assume each column of 'synergies' is a synergy vector
    % Rows = joints, Columns = different synergies from all subjects
    % Example: 15 joints, 12 synergies total (4 subjects × 3 synergies each)

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

    figure;
    dendrogram(Z,0,'ColorThreshold',distances(idxMaxJump+1));
    hold on;
    yline(distances(idxMaxJump), 'r--', 'LineWidth', 1.5);
    title(sprintf('Hierarchical Clustering (nClusters = %d)', numClusters));
    xlabel('Observations'); ylabel('Linkage Distance');
    % legend('Linkage','Cutoff threshold');
    
    %% Step 4: Compute cluster centroids (general synergies)
    generalSynergies = zeros(numClusters, size(allPCs,2));
    
    for k = 1:numClusters
        members = (clusterIdx == k);
        if sum(members) == 0
            continue;
        end
        % Average the original (unnormalized, unaligned) PCs in this cluster
        Ck = mean(allPCs(members,:), 1);
        generalSynergies(k,:) = Ck;
    end
end

function [coeff_general, latent_general, explained_general] = general_2PCA2()
    src = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Participant_PC\PCA_coeffs\';
    src2 = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Activity_PC\PCA_coeffs\';
    src_la = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Participant_PC\Latent\';
    src_la2 = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Activity_PC\Latent\';
    % Stack the first k PCs from each participant
    sub_PC = [];
    task_PC = [];
    latemp = dir(fullfile(src_la,'*.csv'));
    lafolder = {latemp(~[latemp.isdir]).name};
    maintemp = dir(fullfile(src,'*.csv'));
    mainfolder = {maintemp(~[maintemp.isdir]).name};
    for i = 1:numel(mainfolder)
        la_file = fullfile(src_la,lafolder{i});
        laArray = table2array(readtable(la_file));
        laArray = laArray > 1;
        k_new = sum(laArray);
        pc_file = fullfile(src,mainfolder{i});
        thisTable = readtable(pc_file);
        sub_PC = [sub_PC,table2array(thisTable(2:end,1:k_new))];
    end
    latemp2 = dir(fullfile(src_la2,'*.csv'));
    lafolder2 = {latemp2(~[latemp2.isdir]).name};
    maintemp2 = dir(fullfile(src2,'*.csv'));
    mainfolder2 = {maintemp2(~[maintemp2.isdir]).name};
    for i = 1:numel(mainfolder2)
        la_file = fullfile(src_la2,lafolder2{i});
        laArray = table2array(readtable(la_file));
        laArray = laArray > 1;
        k_new = sum(laArray);
        pc_file = fullfile(src2,mainfolder2{i});
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
    % mu_general = mean(Mus,1)';
end

function [coeff_general, latent_general, explained_general] = general_2PCA()
    src = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Participant_PC\PCA_coeffs\';
    src2 = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Activity_PC\PCA_coeffs\';
    % Stack the first 3 PCs from each participant
    k = 3;
    PCs = [];
    maintemp = dir(fullfile(src,'*.csv'));
    mainfolder = {maintemp(~[maintemp.isdir]).name};
    for i = 1:numel(mainfolder)
        file = fullfile(src,mainfolder{i});
        thisTable = readtable(file);
        PCs = [PCs,table2array(thisTable(2:end,1:k))];
    end
    maintemp2 = dir(fullfile(src2,'*.csv'));
    mainfolder2 = {maintemp2(~[maintemp2.isdir]).name};
    for i = 1:numel(mainfolder2)
        file = fullfile(src2,mainfolder2{i});
        thisTable = readtable(file);
        PCs = [PCs,table2array(thisTable(2:end,1:k))];
    end
    [coeff_general, ~, latent_general , ~, explained_general, ~] = pca(PCs', 'Centered', true);
end

function coeff_general = general_PCA_cluster(k)
    src = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\18_par\Participant_PC\PCA_coeffs\';
    src2 = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\18_par\Activity_PC\PCA_coeffs\';
    % src_mu = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\18_par\Participant_PC\PCA_mean\';
    % src_mu2 = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\18_par\Activity_PC\PCA_mean\';
    % src_la = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\18_par\Participant_PC\Latent\';
    % src_la2 = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\18_par\Activity_PC\Latent\';
    % Stack the first k PCs from each participant
    PCs = [];
    % Mus = [];
    % Las = [];
    % mutemp = dir(fullfile(src_mu,'*.csv'));
    % mufolder = {mutemp(~[mutemp.isdir]).name};
    % latemp = dir(fullfile(src_la,'*.csv'));
    % lafolder = {latemp(~[latemp.isdir]).name};
    maintemp = dir(fullfile(src,'*.csv'));
    mainfolder = {maintemp(~[maintemp.isdir]).name};
    for i = 1:numel(mainfolder)
        file = fullfile(src,mainfolder{i});
        thisTable = readtable(file);
        PCs = [PCs,table2array(thisTable(2:end,1:k))];
        % file = fullfile(src_mu,mufolder{i});
        % thisArray = table2array(readtable(file));
        % for j = 1:k
        %     Mus = [Mus,thisArray'];
        % end
        % file = fullfile(src_la,lafolder{i});
        % thisArray = table2array(readtable(file));
        % Las = [Las,thisArray(1:k)'];
    end
    % mutemp2 = dir(fullfile(src_mu2,'*.csv'));
    % mufolder2 = {mutemp2(~[mutemp2.isdir]).name};
    % latemp2 = dir(fullfile(src_la2,'*.csv'));
    % lafolder2 = {latemp2(~[latemp2.isdir]).name};
    maintemp2 = dir(fullfile(src2,'*.csv'));
    mainfolder2 = {maintemp2(~[maintemp2.isdir]).name};
    for i = 1:numel(mainfolder2)
        file = fullfile(src2,mainfolder2{i});
        thisTable = readtable(file);
        PCs = [PCs,table2array(thisTable(2:end,1:k))];
        % file = fullfile(src_mu2,mufolder2{i});
        % thisArray = table2array(readtable(file));
        % for j = 1:k
        %     Mus = [Mus,thisArray'];
        % end
        % file = fullfile(src_la2,lafolder2{i});
        % thisArray = table2array(readtable(file));
        % Las = [Las,thisArray(1:k)'];
    end
    [coeff_general] = cosine_cluster(PCs);
end

function [coeff_general,count1,count2] = general_PCA_cluster2()
    src = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Participant_PC\PCA_coeffs\';
    src2 = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Activity_PC\PCA_coeffs\';
    src_la = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Participant_PC\Latent\';
    src_la2 = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Activity_PC\Latent\';
    % Stack the first k PCs from each participant
    PCs_subject = [];
    PCs_task = [];
    count1 = [];
    count2 = [];
    latemp = dir(fullfile(src_la,'*.csv'));
    lafolder = {latemp(~[latemp.isdir]).name};
    maintemp = dir(fullfile(src,'*.csv'));
    mainfolder = {maintemp(~[maintemp.isdir]).name};
    for i = 1:numel(mainfolder)
        la_file = fullfile(src_la,lafolder{i});
        laArray = table2array(readtable(la_file));
        laArray = laArray > 1;
        k_new = sum(laArray);
        pc_file = fullfile(src,mainfolder{i});
        thisTable = readtable(pc_file);
        PCs_subject = [PCs_subject,table2array(thisTable(2:end,1:k_new))];
        count1 = [count1,k_new];
    end
    latemp2 = dir(fullfile(src_la2,'*.csv'));
    lafolder2 = {latemp2(~[latemp2.isdir]).name};
    maintemp2 = dir(fullfile(src2,'*.csv'));
    mainfolder2 = {maintemp2(~[maintemp2.isdir]).name};
    for i = 1:numel(mainfolder2)
        la_file = fullfile(src_la2,lafolder2{i});
        laArray = table2array(readtable(la_file));
        laArray = laArray > 1;
        k_new = sum(laArray);
        pc_file = fullfile(src2,mainfolder2{i});
        thisTable = readtable(pc_file);
        PCs_task = [PCs_task,table2array(thisTable(2:end,1:k_new))];
        count2 = [count2,k_new];
    end
    coeff_general = cosine_cluster2(PCs_subject,PCs_task);
    coeff_general = coeff_general';
end

% function [coeff_general, latent_general, explained_general,mu_general] = participant_PCA()
%     src = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Participant_hand_PC\PCA_coeffs\';
%     % Stack the first 3 PCs from each participant
%     PCs = [];
%     maintemp = dir(fullfile(src,'*.csv'));
%     mainfolder = {maintemp(~[maintemp.isdir]).name};
%     for i = 1:numel(mainfolder)
%         file = fullfile(src,mainfolder{i});
%         thisTable = readtable(file);
%         PCs = [PCs,table2array(thisTable(2:end,1:3))];
%     end
%     % coeff_general = cosine_cluster(PCs);
%     joint_names = ["TCMC_f","TMCP_f","TIP_f","IMCP_f","IPIP_f","IDIP_f","MMCP_f",...
%         "MPIP_f","MDIP_f","RMCP_f","RPIP_f","RDIP_f","LMCP_f","LPIP_f","LDIP_f",...
%         "TMCP_a","TCMC_r","IMCP_a","MMCP_a","RMCP_a","LMCP_a","W_a","W_f","W_r",...
%         "RE_f","RS_f","RS_a","RS_r"];
%     joint_names2 = ["TCMC_f","TMCP_f","TIP_f","IMCP_f","IPIP_f","IDIP_f","MMCP_f",...
%         "MPIP_f","MDIP_f","RMCP_f","RPIP_f","RDIP_f","LMCP_f","LPIP_f","LDIP_f",...
%         "TMCP_a","TCMC_r","IMCP_a","MMCP_a","RMCP_a","LMCP_a"];
%     for i = 1:4
%         figure
%         bar(joint_names2,coeff_general(:,i));
%         % xticklabels(joint_names);
%         xlabel('Joint');
%         ylabel('Loading');
%         title(['PC ' num2str(i) ' joint contributions']);
%     end
%     [coeff_general, ~, latent_general , ~, explained_general, mu_general] = pca(PCs', 'Centered', true);
% end

% function [coeff_general, latent_general, explained_general, mu_general] = task_PCA()
%     src = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Activity_hand_PC\PCA_coeffs\';
%     % Stack the first 3 PCs from each participant
%     W = [];
%     maintemp = dir(fullfile(src,'*.csv'));
%     mainfolder = {maintemp(~[maintemp.isdir]).name};
%     for i = 1:numel(mainfolder)
%         file = fullfile(src,mainfolder{i});
%         thisTable = readtable(file);
%         W = [W,table2array(thisTable(2:end,1:3))];
%     end
%     % coeff_general = cosine_cluster(W);
%     joint_names = ["TCMC_f","TMCP_f","TIP_f","IMCP_f","IPIP_f","IDIP_f","MMCP_f",...
%         "MPIP_f","MDIP_f","RMCP_f","RPIP_f","RDIP_f","LMCP_f","LPIP_f","LDIP_f",...
%         "TMCP_a","TCMC_r","IMCP_a","MMCP_a","RMCP_a","LMCP_a","W_a","W_f","W_r",...
%         "RE_f","RS_f","RS_a","RS_r"];
%     joint_names2 = ["TCMC_f","TMCP_f","TIP_f","IMCP_f","IPIP_f","IDIP_f","MMCP_f",...
%         "MPIP_f","MDIP_f","RMCP_f","RPIP_f","RDIP_f","LMCP_f","LPIP_f","LDIP_f",...
%         "TMCP_a","TCMC_r","IMCP_a","MMCP_a","RMCP_a","LMCP_a"];
%     for i = 1:3
%         figure
%         bar(joint_names2,coeff_general(:,i));
%         % xticklabels(joint_names);
%         xlabel('Joint');
%         ylabel('Loading');
%         title(['PC ' num2str(i) ' joint contributions']);
%     end
%     [coeff_general, ~, latent_general , ~, explained_general, mu_general] = pca(W', 'Centered', true);
% end

function output = angle2pos(reproarr,file,time,shoulder_focus)
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

function draw_skeleton(reproarr, colour, time, s_focus)
    reproarr = reproarr';
    file = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Bottle\filter\20241108T112102-112115_filtered.csv';
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

function general_synergy(co_gen, mu_gen, k, amp)
    mu_gen = mu_gen(1:21);
    % K is the PC to visualize
    t = 98;
    mean_posture = mu_gen;
    
    % Visualize postures
    posture_plus  = mean_posture + amp * co_gen(:,k);
    % posture_minus = mean_posture - alpha * co_gen(:,k);
    posture_plus = [posture_plus;0;0;0;0;0;0;0];
    mean_posture = [mean_posture;0;0;0;0;0;0;0];
    
    figure
    draw_skeleton(posture_plus,'r', t, 0);
    % draw_skeleton(posture_minus','b');
    draw_skeleton(mean_posture,'k', t, 0);
    title(['Posture variation along PC ' num2str(k)]);
    % legend('','','','','Synergy', 'Mean')

    joint_names = ["TCMC_f","TMCP_f","TIP_f","IMCP_f","IPIP_f","IDIP_f","MMCP_f",...
    "MPIP_f","MDIP_f","RMCP_f","RPIP_f","RDIP_f","LMCP_f","LPIP_f","LDIP_f",...
    "TMCP_a","TCMC_r","IMCP_a","MMCP_a","RMCP_a","LMCP_a","W_a","W_f","W_r",...
    "RE_f","RS_f","RS_a","RS_r"];
    % joint_names2 = ["TCMC_f","TMCP_f","TIP_f","IMCP_f","IPIP_f","IDIP_f","MMCP_f",...
    % "MPIP_f","MDIP_f","RMCP_f","RPIP_f","RDIP_f","LMCP_f","LPIP_f","LDIP_f",...
    % "TMCP_a","TCMC_r","IMCP_a","MMCP_a","RMCP_a","LMCP_a"];
    % figure
    % bar(joint_names,co_gen(:,k),'w');
    % xlabel('Joint');
    % ylabel('Loading');
    % title(['PC ' num2str(k) ' Coefficients'])
    % set(gca,'fontsize',14, 'TickDir', 'out')
end

function synergy(co_gen, la_gen, mu_gen, sigma, k, dev)
    % mu_gen = mu_gen(1:21);
    % sigma = sigma(1:21);
    % K is the PC to visualize
    t = 98;
    amp = dev * sqrt(la_gen(k));  % 2 std deviations (can multiply to increase amplitude)
    mean_posture = mu_gen;
    
    % Visualize postures
    posture_plus  = mean_posture + amp * (co_gen(:,k).*sigma);
    % posture_minus = mean_posture - alpha * co_gen(:,k);
    % posture_plus = [posture_plus;0;0;0;0;0;0;0];
    % mean_posture = [mean_posture;0;0;0;0;0;0;0];
    
    figure
    draw_skeleton(posture_plus,'r', t, 0);
    % draw_skeleton(posture_minus','b');
    draw_skeleton(mean_posture,'k', t, 0);
    title(['Posture variation along PC ' num2str(k)]);
    % legend('','','','','Synergy', 'Mean')

    joint_names = ["TCMC_f","TMCP_f","TIP_f","IMCP_f","IPIP_f","IDIP_f","MMCP_f",...
    "MPIP_f","MDIP_f","RMCP_f","RPIP_f","RDIP_f","LMCP_f","LPIP_f","LDIP_f",...
    "TMCP_a","TCMC_r","IMCP_a","MMCP_a","RMCP_a","LMCP_a","W_a","W_f","W_r",...
    "RE_f","RS_f","RS_a","RS_r"];
    figure
    bar(joint_names,co_gen(:,k));
    xlabel('Joint');
    ylabel('Loading');
    title(['PC ' num2str(k) ' Coefficients'])
    set(gca,'fontsize',14, 'TickDir', 'out')
end

function distance = reco_err(co_gen, mu_gen, sigma, k, time, idx, score_file, file_true)
    try
        % score_file = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Participant_PC\Projection\p27.csv';
        score = table2array(readtable(score_file));
        % idx = 10315;
        score = score(idx,:); %10315 for bottle, 16406 for screw
        recon = (score(:,1:k) * co_gen(:,1:k)')';
        reproarr = recon .* repmat(sigma, 1, size(recon,2)) + repmat(mu_gen, 1, size(recon,2));
        reproarr = reproarr';
        % file = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Bottle\filter\20241108T112102-112115_filtered.csv'; % 20241108T112102-112115 for bottle
        % file_true = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Bottle\angles\20241108T112102-112115_angles.csv'; % 20241106T104644-104656 for screw
        angles_true = table2array(readtable(file_true));
        angles_true = angles_true(time,:);
        %% Drawing reprojection
        % [TMCP;TIP;TT;TMCP;IMCP;MMCP;RMCP;LMCP;IPIP;MPIP;RPIP;LPIP;IDIP;MDIP;RDIP;LDIP;IT;MT;RT;LT;C;S;E;W];
        % pos_guess = angle2pos(reproarr,file,time,1);
        % pos_true = angle2pos(angles_true,file,time,1);
        % error = pos_true-pos_guess;
        % distance = [];
        % for k = 1:size(error,1)
        %     distance = [distance,sqrt(error(k,1)^2+error(k,2)^2+error(k,3)^2)];
        % end
        distance = angles_true - reproarr;
        % figure
        % hold on
        % draw_skeleton(reproarr','r', time, 0);
        % draw_skeleton(angles_true','k', time, 0);
    catch
        distance = zeros(1,28);
    end
end

% function distance = avg_reco_err(co_gen, mu_gen, sigma, k, time)
%     score_file = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Participant_PC\Projection\p27.csv';
%     score = table2array(readtable(score_file));
%     score = score(10315,:);
%     recon = (score(:,1:k) * co_gen(:,1:k)')';
%     reproarr = recon .* repmat(sigma, 1, size(recon,2)) + repmat(mu_gen, 1, size(recon,2));
%     reproarr = reproarr';
%     file = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Bottle\filter\20241108T112102-112115_filtered.csv';
%     file_true = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Bottle\angles\20241108T112102-112115_angles.csv';
%     angles_true = table2array(readtable(file_true));
%     angles_true = angles_true(time,:);
%     %% Drawing reprojection
%     % [TMCP;TIP;TT;TMCP;IMCP;MMCP;RMCP;LMCP;IPIP;MPIP;RPIP;LPIP;IDIP;MDIP;RDIP;LDIP;IT;MT;RT;LT;C;S;E;W];
%     pos_guess = angle2pos(reproarr,file,time,1);
%     pos_true = angle2pos(angles_true,file,time,1);
%     error = pos_true-pos_guess;
%     distance = [];
%     for k = 1:size(error,1)
%         distance = [distance,sqrt(error(k,1)^2+error(k,2)^2+error(k,3)^2)];
%     end
%     figure
%     hold on
%     draw_skeleton(reproarr','r', time, 1);
%     draw_skeleton(angles_true','k', time, 1);
% end

% function grasp_time = PCA_test()
%     src = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Bottle\angles';
%     src_3d = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Bottle\filter';
% 
%     maintemp = dir(fullfile(src,'*.csv'));
%     mainfolder = {maintemp(~[maintemp.isdir]).name};
%     subtemp = dir(fullfile(src_3d,'*.csv'));
%     subfolder = {subtemp(~[subtemp.isdir]).name};
%     categoryAngles = [];
%     start_time = 1;
%     reach_time = [];
%     end_time = [];
%     grasp_time = [];
%     % Looping through participants
%     for i = 1:numel(mainfolder)
% 
%         file = fullfile(src,mainfolder{i});
%         fprintf(1, 'Now reading %s\n', mainfolder{i});
%         thisTable = readtable(file);
% 
%         % Getting moment of grasp
%         grasp_3d = fullfile(src_3d,subfolder{i});
%         grasp_3d = table2array(readtable(grasp_3d));
%         vel = grasp_3d(1:end-1,:)-grasp_3d(2:end,:);
%         vel = sqrt((vel(:,1:3:end)).^2 + (vel(:,2:3:end)).^2 + (vel(:,3:3:end)).^2);
%         [~, vel_rise] = max(vel(1:end-1, 25) > 10 & vel(2:end, 25) > 10);
%         [~, vel_fall] = max(vel(vel_rise:end-1, 25) < 10 & vel(vel_rise+1:end, 25) < 10);
%         reach = vel_fall + vel_rise + 6;
%         grasp_time = [grasp_time;reach];
% 
%         % Finding start and end times of each activity
%         if ~(i == 1)
%             start_time = [start_time;end_time(end)+1];
%         end
%         reach_time = [reach_time;start_time(end)+min([reach,40])-1];
%         end_temp = min(reach+10,height(thisTable));
%         end_time = [end_time;reach_time(end)+min(10,height(thisTable)-reach)];
% 
%         categoryAngles = [categoryAngles;thisTable(max([1,reach-39]):end_temp,:)];
% 
%     end
% end

% p27 bottle reach: [96,67,98,117,107];

% reproarr = [-0.644487314646993,-0.564798001517081,-0.483042678048414,-0.341762097913448,...
%     -0.574457965884443,-0.407124410322963,-0.375727314721728,-0.702834756306882,...
%     -0.470555642570272,-0.378033045829321,-0.737599930294856,-0.507852856439589,...
%     -0.384770519180077,-0.702187811271578,-0.618198219075565,0.0480621283487461,...
%     1.02140100238645,-0.317864752221818,-0.163472360579374,0.00271304425727295,...
%     0.155539752147326,1.58007649686837,0.0338090199801796,2.06019120911894,...
%     -1.04507672128763,-0.532834086106159,-1.05940403170959,2.78456656182536];
% draw_skeleton(reproarr','r');
% 

%% Graphing error distance
close all; clear all; clc;
% writematrix(total_dis,'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\temp\cluster_recon_err.csv')
% writematrix(total_dis,'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\temp\task_recon_err.csv')

src = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\error\p27_tasks';
maintemp = dir(fullfile(src,'*.csv'));
mainfolder = {maintemp(~[maintemp.isdir]).name};
total_dis = [];
for i = 1:(numel(mainfolder)-1)
    file = table2array(readtable(fullfile(src,mainfolder{i})));
    total_dis = [total_dis;file];
end

% total_dis = readmatrix('C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\error\p27_pca2_err.csv');
joint_names = ["TCMC_f","TMCP_f","TIP_f","IMCP_f","IPIP_f","IDIP_f","MMCP_f",...
    "MPIP_f","MDIP_f","RMCP_f","RPIP_f","RDIP_f","LMCP_f","LPIP_f","LDIP_f",...
    "TMCP_a","TCMC_r","IMCP_a","MMCP_a","RMCP_a","LMCP_a","W_a","W_f","W_r",...
    "RE_f","RS_f","RS_a","RS_r"];
meanPC = mean(rad2deg(total_dis));
tot_mean = mean(abs(meanPC));
semPC = std(rad2deg(total_dis))./sqrt(size(total_dis,1));
meanPC = flip(meanPC);
semPC = flip(semPC);
f = figure;
hold on
barh(flip(joint_names),meanPC, 'w');
errorbar(meanPC, 1:28, semPC, 'k.', 'horizontal', 'LineWidth', 1.2, 'CapSize', 10);
ylabel("Joints")
xlabel("Angle error (degrees)")
title("Task Focused Joint Reconstruction Error")
set(gca,'fontsize',14, 'TickDir', 'out')
f.Position = [100 100 400 600];

hand = mean(abs(meanPC(8:end)));
arm = mean(abs(meanPC(1:7)));

%% Characterizing top 4 synergies
[co_gen,count1,count2] = general_PCA_cluster2();
% [co_gen,la_gen,ex_gen] = general_2PCA();
% co_gen = co_gen(:,1:4);
% sigma = table2array(readtable('C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Participant_PC\All\sigma.csv'))';
% mu_gen = table2array(readtable('C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Participant_PC\All\mu_global.csv'))';
% for i = 1:size(co_gen,2)
%     general_synergy(co_gen, mu_gen, i, 2); % 2.5 for clustering, 8 for 2 stage
% end


%% two stage PCA variance explained
% joint_names = ["TCMC_f","TMCP_f","TIP_f","IMCP_f","IPIP_f","IDIP_f","MMCP_f",...
%     "MPIP_f","MDIP_f","RMCP_f","RPIP_f","RDIP_f","LMCP_f","LPIP_f","LDIP_f",...
%     "TMCP_a","TCMC_r","IMCP_a","MMCP_a","RMCP_a","LMCP_a","W_a","W_f","W_r",...
%     "RE_f","RS_f","RS_a","RS_r"];
% for i = 1:5
%     figure
%     bar(joint_names,co_gen(:,i));
%     % xticklabels(joint_names);
%     xlabel('Joint');
%     ylabel('Loading');
%     title(['PC ' num2str(i) ' joint contributions']);
% end
% ex_gen = cumsum(ex_gen);
% figure
% bar(ex_gen, 'w')
% ylim([0 100])
% xlabel('PC');
% ylabel('Variance Explained');

%% Reconstructing example task
% subject_coeff = table2array(readtable('C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Participant_PC\PCA_coeffs\p26.csv'));
% subject_coeff = subject_coeff(2:end,:);
% object_coeff = table2array(readtable('C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Activity_PC\PCA_coeffs\Instructions_WRI4.csv')); % DIN8 and PHY8
% object_coeff = object_coeff(2:end,:);
% distance = reco_err(object_coeff, mu_gen, sigma, 4, 98); % 98 for bottle, 108 for screw

%% Reconstructing all tasks
sigma = table2array(readtable('C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Participant_PC\All\sigma.csv'))';
mu_gen = table2array(readtable('C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Participant_PC\All\mu_global.csv'))';
pro_src = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Participant_PC\Projection\';
time_src = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Participant_PC\PCA_times\';
pc_src = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Participant_PC\PCA_coeffs\';
angle_src = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\angle\';
pctemp = dir(fullfile(pc_src,'*.csv'));
pcfolder = {pctemp(~[pctemp.isdir]).name};
protemp = dir(fullfile(pro_src,'*.csv'));
profolder = {protemp(~[protemp.isdir]).name};
timetemp = dir(fullfile(time_src,'*.csv'));
timefolder = {timetemp(~[timetemp.isdir]).name};
maintemp = dir(fullfile(angle_src,'*'));
mainfolder = setdiff({maintemp([maintemp.isdir]).name},{'.','..'});
total_dis = [];
skip = 0;
count = 0;
for i = 1:numel(mainfolder)
    %% For each participant
    timeTable = readtable(fullfile(time_src,timefolder{i}));
    % pcTable = table2array(readtable(fullfile(pc_src,pcfolder{i})));
    % co_gen = pcTable(2:end,:);
    acttemp = dir(fullfile(angle_src,mainfolder{i},'*'));
    actfolder = setdiff({acttemp([acttemp.isdir]).name},{'.','..'});
    for l = 1:numel(actfolder)
        %% For each activity
        % timeTable = readtable(fullfile(time_src,timefolder{l}));
        % pcTable = table2array(readtable(fullfile(pc_src,pcfolder{l})));
        % co_gen = pcTable(2:end,:);
        subtemp = dir(fullfile(angle_src,mainfolder{i},actfolder{l},'*.csv'));
        subfolder = {subtemp(~[subtemp.isdir]).name};
        for j = 1:numel(subfolder)
            % skip = skip + 1;
            % if skip < 2482
            %     continue
            % end
            file = fullfile(angle_src,mainfolder{i},actfolder{l},subfolder{j});
            % angle_arr = table2array(readtable(file));
            activity = subfolder{j};
            activity = activity(1:end-11);
            reach_time = timeTable.reach_time(strcmp(timeTable.activities, activity));
            grasp_time = timeTable.grasp_time(strcmp(timeTable.activities, activity));
            distance = reco_err(co_gen, mu_gen, sigma, 4, grasp_time, reach_time, fullfile(pro_src,profolder{i}), file);
            total_dis = [total_dis;distance];
        end
    end
end
writematrix(total_dis,'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\error\cluster_recon_err2.csv')

%% Reconstructing unseen task
pro_src = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Participant_PC\Projection\p27.csv';
time_src = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Participant_PC\PCA_times\p27.csv';
pc_src = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Activity_PC\PCA_coeffs\';
angle_src = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\angle\p27';
pctemp = dir(fullfile(pc_src,'*.csv'));
pcfolder = {pctemp(~[pctemp.isdir]).name};
maintemp = dir(fullfile(angle_src,'*'));
mainfolder = setdiff({maintemp([maintemp.isdir]).name},{'.','..'});
timeTable = readtable(time_src);
sigma = table2array(readtable('C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Participant_PC\sigma\p27.csv'))';
mu_gen = table2array(readtable('C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Participant_PC\PCA_mean\p27.csv'))';
% for h = 1:numel(pcfolder)
%     pcTable = table2array(readtable(fullfile(pc_src,pcfolder{h})));
%     co_gen = pcTable(2:end,:);
total_dis = [];
for i = 1:numel(mainfolder)
    subtemp = dir(fullfile(angle_src,mainfolder{i},'*.csv'));
    subfolder = {subtemp(~[subtemp.isdir]).name};
    for j = 1:numel(subfolder)
        file = fullfile(angle_src,mainfolder{i},subfolder{j});
        % angle_arr = table2array(readtable(file));
        activity = subfolder{j};
        activity = activity(1:end-11);
        reach_time = timeTable.reach_time(strcmp(timeTable.activities, activity));
        grasp_time = timeTable.grasp_time(strcmp(timeTable.activities, activity));
        distance = reco_err(co_gen, mu_gen, sigma, 4, grasp_time, reach_time, pro_src, file);
        total_dis = [total_dis;distance];
    end
end
writematrix(total_dis,strcat('C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\error\cluster_recon_err2_p27.csv'));
% end