%% ================================================================
% Time-Resolved EEG Decoding (CSP + SVM)
% Description: Sliding-window decoding of motor imagery EEG
% ================================================================

% clear; clc; close all;

%% -------------------- Load or Simulate Data --------------------
% X: [n_channels x n_samples x n_trials]
% y: [n_trials x 1], e.g. 1 = left hand, 2 = right hand
% fs: sampling frequency (Hz)

% load('EEGdata.mat'); % <-- replace with your own file (must contain X, y, fs)

% X  = yourEEG;   % [channels x samples x trials]
X1 = ALLEEG(2).data;
X2 = ALLEEG(3).data;
X = cat(3,X1,X2);
class1 = 'G'; %Grasp
class2 = 'P'; %Pinch
% y  = yourLabels; % [trials x 1], integers (e.g., 1,2)
y = [repmat(class1,size(X1,3),1);repmat(class2,size(X2,3),1)];
fs = 1000;        % or whatever your sampling rate is
% save('EEGdata.mat', 'X', 'y', 'fs');

[nCh, nSamp, nTrials] = size(X);
fprintf('Loaded data: %d channels, %d samples, %d trials\n', nCh, nSamp, nTrials);

%% -------------------- Band-pass Filter --------------------
% filt_band = [8 30]; % mu + beta band
% [b, a] = butter(4, filt_band/(fs/2));
% 
% for i = 1:30
%     X1(:,:,i) = filtfilt(b, a, double(X1(:,:,i))')';
% end
% for i = 1:30
%     X2(:,:,i) = filtfilt(b, a, double(X2(:,:,i))')';
% end
% 
% X = cat(3,X1,X2);
% [nCh, nSamp, nTrials] = size(X);
% fprintf('Loaded data: %d channels, %d samples, %d trials\n', nCh, nSamp, nTrials);


%% -------------------- CSP Function --------------------
function [W] = computeCSP(X1,X2)
    class1 = X1;
    class2 = X2;

    cov1 = zeros(size(X1,1));
    cov2 = zeros(size(X2,1));

    for i = 1:size(class1,3)
        C = cov(class1(:,:,i)');
        cov1 = cov1 + C / trace(C);
    end
    for i = 1:size(class2,3)
        C = cov(class2(:,:,i)');
        cov2 = cov2 + C / trace(C);
    end

    cov1 = cov1 / size(class1,3);
    cov2 = cov2 / size(class2,3);
    composite = cov1 + cov2;

    [EVecs, EVals] = eig(cov1, composite);
    [~, ind] = sort(diag(EVals), 'descend');
    W = EVecs(:, ind);
end

%% -------------------- Sliding Window Parameters --------------------
win_len = 0.25;  % seconds
win_step = 0.05; % seconds
win_samp = round(win_len * fs);
step_samp = round(win_step * fs);

time_idx = 1:step_samp:(nSamp - win_samp);
acc_time = zeros(1, length(time_idx));

%% -------------------- Time-Resolved Decoding Loop --------------------
nCSP = 3; % number of CSP filters from each end
true_all = [];
pred_all = [];

for t = 1:length(time_idx)
    idx = time_idx(t):(time_idx(t)+win_samp-1);
    Xwin1 = X1(:, idx, :);
    Xwin2 = X2(:, idx, :);
    Xwin = X(:, idx, :);

    % Compute CSP filters for this time window
    W = computeCSP(Xwin1, Xwin2);
    Wf = [W(:,1:nCSP) W(:,end-nCSP+1:end)];

    % Extract CSP log-variance features
    feats = zeros(nTrials, 2*nCSP);
    for i = 1:nTrials
        Xf = Wf' * Xwin(:,:,i);
        var_feat = var(Xf,0,2);
        feats(i,:) = log(var_feat / sum(var_feat));
    end

    % Train/test split
    cv = cvpartition(y, 'HoldOut', 0.2);
    Xtrain = feats(training(cv), :);
    Ytrain = y(training(cv));
    Xtest  = feats(test(cv), :);
    Ytest  = y(test(cv));

    % Train SVM (linear)
    model = fitcsvm(Xtrain, Ytrain, ...
        'KernelFunction', 'linear', ...
        'Standardize', true, ...
        'BoxConstraint', 1);

    % Predict and compute accuracy
    Ypred = predict(model, Xtest);
    acc_time(t) = mean(Ypred == Ytest);

    pred_all = [pred_all; Ypred];
    true_all = [true_all; Ytest];
end

%% -------------------- Plot Confusion Matrix --------------------
figure;
confusionchart(true_all,pred_all,'Normalization','row-normalized');
titlestr = strcat("activity comparison");
title(titlestr)
fontsize(18,"points")

%% -------------------- Plot Decoding Accuracy --------------------
time_axis = time_idx / fs;
figure;
plot(time_axis, acc_time * 100, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Accuracy (%)');
title('Motor Classification Over Time (CSP + SVM)');
grid on;
ylim([40 100]);
yline(50, '--r', 'Chance Level');

%% -------------------- Summary --------------------
[~, bestIdx] = max(acc_time);
fprintf('Peak decoding accuracy: %.2f%% at %.2f s\n', acc_time(bestIdx)*100, time_axis(bestIdx));