close all; clear all; clc;

function mean_PC_sig()
    %% ==========================================================
    %  Test if mean PC loadings across subjects differ from 0
    %  ==========================================================
    %  Input:
    %     L : [nSubjects x nPC x nJoints]
    %         Each L(s,p,j) = loading for subject s, PC p, joint j
    %  Output:
    %     pVals_rank : p-values for each PC
    %  ==========================================================
    src = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Activity_PC\PCA_coeffs';
    maintemp = dir(fullfile(src,'*.csv'));
    mainfolder = {maintemp(~[maintemp.isdir]).name};
    nSubjects = numel(mainfolder);
    nJoints = 28;
    nPC = 4;
    for pcIdx = 1:nPC
        allPCs = zeros(nJoints, nSubjects);
        for s = 1:nSubjects
            coeff = table2array(readtable(fullfile(src,mainfolder{s})));
            allPCs(:,s) = coeff(2:end,pcIdx);
        end
        pVals_joint = zeros(size(allPCs,1),1);
        meanLoad  = zeros(size(allPCs,1),1);
        
        for j = 1:size(allPCs,1)
            jointVals = allPCs(j,:);             % loadings for this joint across subjects
            meanLoad(j) = mean(jointVals);
            pVals_joint(j) = ranksum(jointVals, zeros(size(jointVals)));
        end
        
        % --- Optional multiple-comparison correction (Bonferroni)
        pVals_adj = pVals_joint * size(allPCs,1);
        pVals_adj(pVals_adj > 1) = 1;
        
        % --- Display summary for this PC
        fprintf('\nSignificance test for PC%d joint loadings:\n', pcIdx);
        fprintf('-------------------------------------------\n');
        sigCount = sum(pVals_adj < 0.05);
        fprintf('%d / %d joints significant (p < 0.05, Bonferroni-corrected)\n', ...
                sigCount, numel(pVals_adj));
        
        % --- (Optional) visualize which joints are significant
        figure('Color','w');
        bar(meanLoad,'FaceColor',[0.7 0.7 0.9]); hold on;
        sigIdx = find(pVals_adj < 0.05);
        plot(sigIdx, meanLoad(sigIdx), 'r*', 'MarkerSize',8);
        xlabel('Joint Index'); ylabel('Median Loading');
        title(sprintf('PC%d: Joint Loadings vs 0 (Wilcoxon rank-sum)', pcIdx));
        grid on;
    end
end


function Mean_PC_permutation_sig()
    % ------------------------------
    % 4. VARIANCE OF PC LOADINGS ACROSS SUBJECTS (BAR GRAPH)
    % ------------------------------
    src = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Participant_PC\PCA_coeffs';
    maintemp = dir(fullfile(src,'*.csv'));
    mainfolder = {maintemp(~[maintemp.isdir]).name};
    nSubjects = numel(mainfolder);
    nJoints = 28;
    nPC = 4;
    meanCorr = zeros(nPC,1);
    pVals_perm = zeros(nPC,1);
    useCosine = true;   % set false to use Pearson correlation
    alignSign = true;   % align sign direction before comparison
    nBoot = 1000;   % number of permutations
    for pcIdx = 1:nPC
        allPCs = zeros(nJoints, nSubjects);
        for s = 1:nSubjects
            coeff = table2array(readtable(fullfile(src,mainfolder{s})));
            allPCs(:,s) = coeff(2:end,pcIdx);
        end
        allPCs = allPCs';
        % --- Align signs so that all PCs have positive mean or same reference
        if alignSign
            ref = mean(allPCs,1);
            for s = 1:nSubjects
                if dot(allPCs(s,:), ref) < 0
                    allPCs(s,:) = -allPCs(s,:);
                end
            end
        end
    
        % ------------------------------------------------------
        % 2️⃣ Compute pairwise similarity (cosine or correlation)
        % ------------------------------------------------------
        if useCosine
            % normalize each vector
            mat_norm = allPCs ./ vecnorm(allPCs,2,2);
            simMat = mat_norm * mat_norm';  % cosine similarity
        else
            simMat = corr(allPCs');
        end
    
        meanSim(pcIdx) = mean(simMat(tril(true(nSubjects),-1)), 'all');
    
        % ------------------------------------------------------
        % 3️⃣ Permutation test: shuffle joint order
        % ------------------------------------------------------
        bootSim = zeros(nBoot,1);
        for b = 1:nBoot
            shuff = allPCs(:, randperm(nJoints));
            if useCosine
                shuff_norm = shuff ./ vecnorm(shuff,2,2);
                s = shuff_norm * shuff_norm';
            else
                s = corr(shuff');
            end
            bootSim(b) = mean(s(tril(true(nSubjects),-1)), 'all');
        end
    
        % small p → observed similarity higher than random
        pVals_perm(pcIdx) = mean(bootSim >= meanSim(pcIdx));
    end
    
    %% ----------------------------------------------------------
    % 4️⃣ Visualization
    % ----------------------------------------------------------
    figure('Name','Inter-subject PC Shape Similarity (Aligned)','Color','w');
    yyaxis left
    bar(1:nPC, meanSim, 'FaceColor',[0.7 0.7 0.9]);
    ylabel('Mean Inter-subject Similarity');
    yyaxis right
    plot(1:nPC, pVals_perm, '-ok','MarkerFaceColor','k');
    ylabel('Permutation p-value');
    xlabel('Principal Component #');
    title('Shape Similarity of PC Loadings Across Subjects (Aligned)');
    grid on; set(gca, 'XTick', 1:nPC);
    
    % Add significance markers
    hold on;
    for pcIdx = 1:nPC
        if pVals_perm(pcIdx) < 0.001
            stars = '***';
        elseif pVals_perm(pcIdx) < 0.01
            stars = '**';
        elseif pVals_perm(pcIdx) < 0.05
            stars = '*';
        else
            stars = '';
        end
        text(pcIdx, meanSim(pcIdx)+0.02, stars, 'HorizontalAlignment','center', ...
             'FontSize',12,'FontWeight','bold','Color','r');
    end
    hold off;
    
    %% ----------------------------------------------------------
    % 5️⃣ Output summary
    % ----------------------------------------------------------
    % fprintf('\nPC Shape Similarity Results (Aligned + %s)\n', ...
    %         useCosine*'Cosine' + ~useCosine*'Correlation');
    fprintf('----------------------------------------------------------\n');
    for pcIdx = 1:nPC
        fprintf('PC%-2d: meanSim = %.3f, p = %.4f\n', pcIdx, meanSim(pcIdx), pVals_perm(pcIdx));
    end
    
    fprintf('\nInterpretation:\n');
    fprintf('  PCs with high mean similarity and low p (<0.05)\n');
    fprintf('  have consistent shape patterns across subjects (general synergies).\n');
    fprintf('  High similarity but high p means shapes look similar but not significantly\n');
    fprintf('  more similar than random — possible shared structure but not confirmed.\n');
end



function recon_err_sig()
    %% ==========================================================
    %  SIGNIFICANCE TESTS: Comparing reconstruction error across 4 PC sets
    %  ==========================================================
    %  Inputs:
    %   err_general1  [nTrial x nJoint]  - general PCs (e.g., clustering)
    %   err_general2  [nTrial x nJoint]  - general PCs (e.g., 2-stage PCA)
    %   err_specific1 [nTrial x nJoint]  - specific PCs (method 1)
    %   err_specific2 [nTrial x nJoint]  - specific PCs (method 2)
    %
    %  Smaller error = better reconstruction
    %  ==========================================================
    
    par_src = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\error\p27';
    maintemp = dir(fullfile(par_src,'*.csv'));
    mainfolder = {maintemp(~[maintemp.isdir]).name};
    err_specific1 = [];
    for i = 1:numel(mainfolder)
        file = rad2deg(readmatrix(fullfile(par_src,mainfolder{i})));
        % err_specific1 = [err_specific1;file];
        err_specific1 = cat(3, err_specific1, file);
    end
    task_src = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\error\p27_tasks';
    maintemp = dir(fullfile(task_src,'*.csv'));
    mainfolder = {maintemp(~[maintemp.isdir]).name};
    err_specific2 = [];
    for i = 1:numel(mainfolder)
        file = rad2deg(readmatrix(fullfile(task_src,mainfolder{i})));
        % err_specific2 = [err_specific2;file];
        err_specific2 = cat(3, err_specific2, file);
    end
    err_specific1 = squeeze(mean(err_specific1, 3)); % average over subjects
    err_specific2 = squeeze(mean(err_specific2, 3)); % average over tasks
    err_general1 = rad2deg(readmatrix('C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\error\p27_pca2_err.csv'));
    err_general2 = rad2deg(readmatrix('C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\error\cluster_recon_err2_p27.csv'));
  
    
    % err_general1 = mean(mean(abs(err_general1)));
    % err_general2 = mean(mean(abs(err_general2)));
    % err_specific1 = mean(mean(abs(err_specific1)));
    % err_specific2 = mean(mean(abs(err_specific2)));
    
    % joint_names = ["TCMC_f","TMCP_f","TIP_f","IMCP_f","IPIP_f","IDIP_f","MMCP_f",...
    %     "MPIP_f","MDIP_f","RMCP_f","RPIP_f","RDIP_f","LMCP_f","LPIP_f","LDIP_f",...
    %     "TMCP_a","TCMC_r","IMCP_a","MMCP_a","RMCP_a","LMCP_a","W_a","W_f","W_r",...
    %     "RE_f","RS_f","RS_a","RS_r"];
    % meanPC = mean(err_specific2);
    % tot_mean = mean(abs(meanPC));
    % semPC = std(err_specific2)./sqrt(size(err_specific2,1));
    % meanPC = flip(meanPC);
    % semPC = flip(semPC);
    % f = figure;
    % hold on
    % barh(flip(joint_names),meanPC);
    % errorbar(meanPC, 1:28, semPC, 'k.', 'horizontal', 'LineWidth', 1.2, 'CapSize', 10);
    % xlabel("Joints")
    % ylabel("Angle error (degrees)")
    % title("Task Focused Joint Reconstruction Error")
    % set(gca,'fontsize',14, 'TickDir', 'out')
    % f.Position = [100 100 400 600];
    
    % err_general1 = rad2deg(readmatrix('C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\error\pca2_recon_err.csv'));
    % err_general2 = rad2deg(readmatrix('C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\error\cluster_recon_err2.csv'));
    % err_specific1 = rad2deg(readmatrix('C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\error\sub_recon_err.csv'));
    % err_specific2 = rad2deg(readmatrix('C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\error\task_recon_err.csv'));
    % 
    % err_general1 = err_general1(1:size(err_specific2,1),:);
    % err_general2 = err_general2(1:size(err_specific2,1),:);
    % err_specific1 = err_specific1(1:size(err_specific2,1),:);
    
    % ------------------------------------------------------------
    % 1️⃣ Compute mean reconstruction error per trial
    % ------------------------------------------------------------
    meanErr.general1  = mean(err_general1, 2);  % average over joints
    meanErr.general2  = mean(err_general2, 2);
    meanErr.specific1 = mean(err_specific1, 2);
    meanErr.specific2 = mean(err_specific2, 2);
    
    % Collect for convenience
    allErr = [meanErr.general1, meanErr.general2, meanErr.specific1, meanErr.specific2];
    labels = {'General_1','General_2','Specific_1','Specific_2'};
    nSets = numel(labels);
    nTrials = size(allErr,1);
    
    % ------------------------------------------------------------
    % 2️⃣ Check normality (per error set)
    % ------------------------------------------------------------
    isNormal = false(1, nSets);
    for i = 1:nSets
        thisData = allErr(:,i);
        if all(isfinite(thisData)) && numel(unique(thisData)) > 1
            [h,p] = lillietest(thisData);
            isNormal(i) = (p > 0.05);
            if isNormal(i)
                fprintf('%s: Lilliefors p = %.4f -> Normal\n', labels{i}, p);
            else
                fprintf('%s: Lilliefors p = %.4f -> Non-normal\n', labels{i}, p);
            end
        else
            fprintf('%s: Skipped normality test (insufficient variation)\n', labels{i});
        end
    end
    
    useParametric = all(isNormal);
    if useParametric
        disp('All groups approximately normal -> using paired t-tests.');
    else
        disp('Non-normal data detected -> using Wilcoxon signed-rank tests.');
    end
    
    %% --- Pairwise comparisons (unchanged logic) ---
    pairs = [1 3; 1 4; 2 3; 2 4];
    pVals = zeros(size(pairs,1),1);
    
    fprintf('\nPairwise comparisons between general and specific PC sets:\n');
    for ii = 1:size(pairs,1)
        A = allErr(:, pairs(ii,1));
        B = allErr(:, pairs(ii,2));
        if useParametric
            [~, pVals(ii)] = ttest(A, B, 'Tail','right'); % general worse?
            testType = 'paired t-test';
        else
            pVals(ii) = signrank(A, B, 'Tail','right');
            testType = 'Wilcoxon signed-rank';
        end
        fprintf('  %s vs %s -> p = %.5f (%s)\n', ...
            labels{pairs(ii,1)}, labels{pairs(ii,2)}, pVals(ii), testType);
    end
    
    % Bonferroni correction
    pAdj = min(pVals * size(pairs,1), 1);
    fprintf('\nBonferroni-corrected p-values:\n');
    for ii = 1:numel(pAdj)
        fprintf('  %s vs %s -> p = %.5f\n', labels{pairs(ii,1)}, labels{pairs(ii,2)}, pAdj(ii));
    end
    
    %% --- Overall comparison plot (reshape to long; fix boxchart input) ---
    xCats = repelem(1:nSets, nTrials)';   %  [nSets*nTrials x 1]
    yVals = allErr(:);                    %  [nSets*nTrials x 1]
    
    figure; hold on;
    if exist('boxchart','file')
        boxchart(xCats, yVals, 'BoxFaceColor',[0.7 0.7 0.9]);
    else
        % Fallback for older MATLAB: use boxplot
        boxplot(yVals, xCats, 'Colors','k', 'Symbol','k.');
    end
    
    % jittered scatter for individual trials
    scatter(xCats + 0.02*randn(size(xCats)), yVals, 15, 'k', 'filled', 'MarkerFaceAlpha',0.4);
    
    set(gca, 'XTick', 1:nSets, 'XTickLabel', labels);
    ylabel('Mean Reconstruction Error (per trial)');
    title('Reconstruction Performance Across PC Sets');
    grid on;
    
    %% --- Per-subject paired plot (example: General_1 vs Specific_1) ---
    % If you have subject IDs per trial, average first; otherwise plots per-trial pairs.
    % OPTIONAL averaging by subject (uncomment if you have subjectID vector [nTrials x 1]):
    % subjIDs = ... % e.g., integers 1..nSubjects for each trial
    % A = splitapply(@mean, allErr(:,1), subjIDs); % General_1 per subject
    % B = splitapply(@mean, allErr(:,3), subjIDs); % Specific_1 per subject
    % lbl = 'Per-Subject';
    
    % Default: per-trial paired lines
    A = allErr(:,1);  % General_1
    B = allErr(:,3);  % Specific_1
    lbl = 'Per-Trial';
    
    figure; hold on;
    for i = 1:numel(A)
        plot([1 2], [A(i) B(i)], '-', 'Color',[0.6 0.6 0.6]);
    end
    meanA = mean(A); meanB = mean(B);
    plot([1 2], [meanA meanB], '-or', 'LineWidth',2, 'MarkerFaceColor','r');
    set(gca, 'XTick', [1 2], 'XTickLabel', {'General_1','Specific_1'});
    ylabel('Mean Reconstruction Error');
    title([lbl ' Reconstruction Error Comparison']);
    grid on;
    legend('Individual pairs','Mean difference','Location','best');
end


recon_err_sig()

% close all; clear all; clc;
% angle_src = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\angle\';
% time_src = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Participant_PC\PCA_times\';
% timetemp = dir(fullfile(time_src,'*.csv'));
% timefolder = {timetemp(~[timetemp.isdir]).name};
% maintemp = dir(fullfile(angle_src,'*'));
% mainfolder = setdiff({maintemp([maintemp.isdir]).name},{'.','..'});
% all_angles = [];
% for i = 1:numel(mainfolder)
%     timeTable = readtable(fullfile(time_src,timefolder{i}));
%     acttemp = dir(fullfile(angle_src,mainfolder{i},'*'));
%     actfolder = setdiff({acttemp([acttemp.isdir]).name},{'.','..'});
%     for l = 1:numel(actfolder)
%         subtemp = dir(fullfile(angle_src,mainfolder{i},actfolder{l},'*.csv'));
%         subfolder = {subtemp(~[subtemp.isdir]).name};
%         for j = 1:numel(subfolder)
%             file = fullfile(angle_src,mainfolder{i},actfolder{l},subfolder{j});
%             angle_arr = table2array(readtable(file));
%             activity = subfolder{j};
%             activity = activity(1:end-11);
%             grasp_time = timeTable.grasp_time(strcmp(timeTable.activities, activity));
%             try
%                 % all_angles = [all_angles;angle_arr((grasp_time-39):(grasp_time+10),:)];
%                 all_angles = [all_angles;angle_arr(grasp_time,:)];
%             catch
%                 continue
%             end
%         end
%     end
% end
% qf_bartlett(all_angles)