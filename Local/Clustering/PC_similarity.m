%% inter_subject_PC_similarity.m
% Compare PCA subspaces between participants using principal angles

close all; clear all; clc;

% ---- USER INPUTS ----
k = 4;  % number of PCs to use per subject (adjust as needed)
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

% ---- VISUALIZE HEATMAP ----
figure;
imagesc(similarityMat);
axis equal tight;
colorbar;
clim([0 1]);
title(sprintf('Inter-Subject Subspace Similarity (Top %d PCs)', k));
xlabel('Participant'); ylabel('Participant');

D = 1 - similarityMat;  % distance matrix
D = (D + D')/2;
Y = mdscale(D, 2, 'Criterion', 'stress');

figure;
scatter(Y(:,1), Y(:,2), 80, 'filled');
text(Y(:,1)+0.02, Y(:,2), arrayfun(@(x) sprintf('S%d',x), 1:nSubjects, 'UniformOutput', false));
title(sprintf('MDS Projection of Inter-Subject Synergy Similarity (Top %d PCs)', k));
xlabel('Dimension 1'); ylabel('Dimension 2');
grid on; axis equal;

% ---- SUMMARY ----
meanSim = mean(similarityMat(triu(true(size(similarityMat)),1)), 'omitnan');
fprintf('Mean inter-subject subspace similarity: %.2f\n', meanSim);

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
    stdPC  = std(allPCs, 0, 2);

    % --- BAR GRAPH WITH ERROR BARS ---
    figure('Name', sprintf('PC%d Loading Variability', pcIdx), 'Color', 'w');
    bar(meanPC, 'FaceColor', [0.2 0.45 0.8], 'EdgeColor', 'none'); hold on;
    errorbar(1:nJoints, meanPC, stdPC, 'k.', 'LineWidth', 1.2, 'CapSize', 10);
    xticklabels({"TCMC f","TMCP f","TIP f","IMCP f","IPIP f","IDIP f","MMCP f",...
    "MPIP f","MDIP f","RMCP f","RPIP f","RDIP f","LMCP f","LPIP f","LDIP f",...
    "TMCP a","TCMC r","IMCP a","MMCP a","RMCP a","LMCP a","W a","W f","W r",...
    "RE f","RS f","RS a","RS r"})
    
    xlabel('Joint Index');
    ylabel('Mean Loading Weight');
    title(sprintf('Variance of PC%d Loadings Across Participants', pcIdx));
    grid on;
    box off;

    % improve x-axis readability
    xticks(1:nJoints);
    xlim([0 nJoints+1]);
end

% ------------------------------
% 5. SUMMARY OUTPUT
% ------------------------------
fprintf('\n=== Summary ===\n');
fprintf('Participants: %d\n', nSubjects);
fprintf('Joints per PC: %d\n', nJoints);
fprintf('Mean Inter-Subject Similarity: %.2f\n', meanSim);
fprintf('Displayed Variability for Top %d PCs.\n', k);