%% inter_subject_PC_similarity.m
% Compare PCA subspaces between participants using principal angles

close all; clear all; clc;
% src = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Activity_PC\Latent';
% maintemp = dir(fullfile(src,'*.csv'));
% mainfolder = {maintemp(~[maintemp.isdir]).name};
% latent = [];
% for i = 1:numel(mainfolder)
%     latent = [latent,table2array(readtable(fullfile(src,mainfolder{i})))];
% end
% latent = mean(latent,2);
% latent2 = latent > 1;

% src = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Activity_PC\PCA_explained';
% maintemp = dir(fullfile(src,'*.csv'));
% mainfolder = {maintemp(~[maintemp.isdir]).name};
% Ex = [];
% nsub = numel(mainfolder);
% for i = 1:nsub
%     Ex = [Ex,table2array(readtable(fullfile(src,mainfolder{i})))];
% end
% Ex_mean = mean(Ex,2);
% Ex_std = std(Ex,0,2);
% 
% figure
% bar(Ex_mean, 'w'); hold on;
% errorbar(1:size(Ex_mean), Ex_mean, Ex_std, 'k.', 'LineWidth', 1.2, 'CapSize', 10);
% xlabel('PCs');
% ylabel('Cumulative Variance Explained');
% title('Activities');
% % grid on;
% box off;
% set(gca,'fontsize',14, 'TickDir', 'out')

% ---- USER INPUTS ----
k = 7;  % number of PCs to use per subject (adjust as needed)
% nSubjects = 10; % total number of participants
% Example: each cell contains [nJoints x nPCs] loadings
% PCs{1}, PCs{2}, ..., PCs{nSubjects}
% Replace this with your real PCA outputs
% load('subject_PCs.mat', 'PCs');  

% ---- COMPUTE PAIRWISE SUBSPACE ANGLES ----
src = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Participant_PC\PCA_coeffs';
maintemp = dir(fullfile(src,'*.csv'));
mainfolder = {maintemp(~[maintemp.isdir]).name};
nSubjects = numel(mainfolder);
subspaceAngles = zeros(nSubjects,nSubjects);
for i = 1:nSubjects
    for j = 1:nSubjects
        if i ~= j
            coeffA = table2array(readtable(fullfile(src,mainfolder{i})));
            coeffA = coeffA(2:end,1:k);
            coeffB = table2array(readtable(fullfile(src,mainfolder{j})));
            coeffB = coeffB(2:end,1:k);

            % orthonormalize (important for numerical stability)
            [QA, ~] = qr(coeffA, 0);
            [QB, ~] = qr(coeffB, 0);

            % singular values of QA' * QB give cosines of principal angles
            s = svd(QA' * QB);
            s = min(max(s, -1), 1); % clip numerical errors
            angles = real(acos(s));

            % average principal angle between subspaces
            subspaceAngles(i,j) = mean(angles);
        else
            subspaceAngles(i,j) = 0;
        end
    end
end

% ---- CONVERT TO SIMILARITY MATRIX ----
% Similarity = 1 - normalized angle
maxAngle = pi/2; % maximum possible subspace angle
similarityMat = 1 - (subspaceAngles / maxAngle);

%% ==========================================================
%  SIGNIFICANCE TESTS ON INTER-SUBJECT PCA SIMILARITY
%  Using t-test and ANOVA to detect significantly different subjects
%  ==========================================================
% 
nSubs = size(similarityMat, 1);
D = 1 - similarityMat;          % convert to dissimilarity if desired
meanSim = mean(similarityMat - diag(diag(similarityMat)), 2);

% --- 1️⃣ One-sample t-test: each subject vs others
fprintf('\nOne-sample t-test for individual subject similarity:\n');
pVals_t = nan(nSubs,1);
tStats = nan(nSubs,1);

for i = 1:nSubs
    others = setdiff(1:nSubs, i);
    subjSim = similarityMat(i, others);
    groupSim = mean(similarityMat(others, others), 'all'); % grand mean of others

    [~, pVals_t(i), ~, stats] = ttest(subjSim, groupSim);
    tStats(i) = stats.tstat;
end

% Bonferroni correction
alpha = 0.05 / nSubs;
outliers_t = find(pVals_t < alpha);

fprintf('Significant outliers (t-test, Bonferroni corrected): %s\n', mat2str(outliers_t));
fprintf('Mean similarity per subject (M±SD): %.3f ± %.3f\n', mean(meanSim), std(meanSim));

% --- Visualization
figure; 
bar(meanSim,'FaceColor',[0.6 0.6 0.9]);
hold on;
xline(outliers_t,'r--','LineWidth',1.5);
xlabel('Subject'); ylabel('Mean similarity to others');
title('Subject-wise Mean Similarity (t-test outliers in red)');


% ---- VISUALIZE HEATMAP ----
figure;
imagesc(similarityMat);
axis equal tight;
colorbar;
clim([0 1]);
title(sprintf('Inter-Task Subspace Similarity (Top %d PCs)', k));
xlabel('Task'); ylabel('Task');
set(gca,'fontsize',14, 'TickDir', 'out')

D = 1 - similarityMat;  % distance matrix
D = (D + D')/2;
Y = mdscale(D, 2, 'Criterion', 'stress');

figure;
scatter(Y(:,1), Y(:,2), 80, 'MarkerEdgeColor',[0 0 0]);
text(Y(:,1)+0.02, Y(:,2), arrayfun(@(x) sprintf('S%d',x), 1:nSubjects, 'UniformOutput', false));
title(sprintf('MDS Projection of Inter-Task Synergy Similarity (Top %d PCs)', k));
xlabel('Dimension 1'); ylabel('Dimension 2');
% grid on;
axis equal;
set(gca,'fontsize',14, 'TickDir', 'out')
% [10:(nSubjects+7),8,9]

% ---- SUMMARY ----
meanSim = mean(similarityMat(triu(true(size(similarityMat)),1)), 'omitnan');
fprintf('Mean inter-task subspace similarity: %.2f\n', meanSim);

% ------------------------------
% 4. VARIANCE OF PC LOADINGS ACROSS SUBJECTS (BAR GRAPH)
% ------------------------------
nJoints = 28;
for pcIdx = 1:k
    allPCs = zeros(nJoints, nSubjects);
    for s = 1:nSubjects
        coeff = table2array(readtable(fullfile(src,mainfolder{s})));
        allPCs(:,s) = coeff(2:end,pcIdx);
    end
    
    meanPC = mean(allPCs, 2);
    stdPC  = std(allPCs, 0, 2)./sqrt(size(allPCs,2)); % Getting SEM rather than std

    % --- BAR GRAPH WITH ERROR BARS ---
    figure('Name', sprintf('PC%d Loading Variability', pcIdx), 'Color', 'w');
    bar(meanPC, 'w'); hold on;
    errorbar(1:nJoints, meanPC, stdPC, 'k.', 'LineWidth', 1.2, 'CapSize', 10);
    xticklabels({"TCMC f","TMCP f","TIP f","IMCP f","IPIP f","IDIP f","MMCP f",...
    "MPIP f","MDIP f","RMCP f","RPIP f","RDIP f","LMCP f","LPIP f","LDIP f",...
    "TMCP a","TCMC r","IMCP a","MMCP a","RMCP a","LMCP a","W a","W f","W r",...
    "RE f","RS f","RS a","RS r"})
    
    xlabel('Joint Index');
    ylabel('Mean Loading Weight');
    title(sprintf('Variance of PC%d Loadings Across Tasks', pcIdx));
    % grid on;
    box off;
    set(gca,'fontsize',14, 'TickDir', 'out')

    % improve x-axis readability
    xticks(1:nJoints);
    xlim([0 nJoints+1]);
end

% ------------------------------
% 5. SUMMARY OUTPUT
% ------------------------------
fprintf('\n=== Summary ===\n');
fprintf('Tasks: %d\n', nSubjects);
fprintf('Joints per PC: %d\n', nJoints);
fprintf('Mean Inter-Task Similarity: %.2f\n', meanSim);
fprintf('Displayed Variability for Top %d PCs.\n', k);