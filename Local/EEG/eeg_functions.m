function time_resolved_decoding()
    win_len = 0.2;      % 200 ms window
    win_step = 0.05;    % 50 ms step
    fs = 250;           % Sampling rate
    win_samp = round(win_len * fs);
    step_samp = round(win_step * fs);
    
    nSamples = size(X,2);
    nTrials = size(X,3);
    
    times = 1:step_samp:(nSamples - win_samp);
    acc_time = zeros(1, length(times));
    
    for t = 1:length(times)
        idx = times(t):(times(t)+win_samp-1);
        Xwin = X(:, idx, :);
        
        % Extract CSP features
        W = csp(Xwin, y);
        Wf = [W(:,1:3) W(:,end-2:end)];
        
        feats = zeros(nTrials, 6);
        for i = 1:nTrials
            Xf = Wf' * Xwin(:,:,i);
            var_feat = var(Xf,0,2);
            feats(i,:) = log(var_feat / sum(var_feat));
        end
        
        % Train/test split
        cv = cvpartition(y, 'HoldOut', 0.2);
        model = fitcsvm(feats(training(cv),:), y(training(cv)), 'KernelFunction', 'linear');
        Ypred = predict(model, feats(test(cv),:));
        acc_time(t) = mean(Ypred == y(test(cv)));
    end
    
    % Plot decoding accuracy over time
    time_axis = times / fs;
    plot(time_axis, acc_time, 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('Accuracy');
    title('Time-resolved decoding accuracy');
    grid on;
end



function continuous_decoding()
    % X: EEG features [n_samples x n_channels]
    % y: continuous target [n_samples x 1]
    
    fs = 250;
    [b,a] = butter(4, [0.5 30]/(fs/2));
    Xf = filtfilt(b,a,X); % Filter EEG
    
    % Extract band power features in sliding windows
    win = 0.2 * fs; step = 0.05 * fs;
    for i = 1:step:(length(y)-win)
        seg = Xf(i:i+win-1,:);
        feats(i,:) = log(var(seg));   % log-bandpower
        targets(i) = mean(y(i:i+win-1));
    end
    
    % Fit Support Vector Regressor
    mdl = fitrsvm(feats, targets, 'KernelFunction', 'linear');
    
    % Predict over time
    yhat = predict(mdl, feats);
    plot(yhat); hold on; plot(targets); legend('Predicted','True');
    title('Continuous EEG decoding');
end


%% ================================================================
% Time-Resolved EEG Decoding (CSP + SVM)
% Author: ChatGPT (GPT-5)
% Description: Sliding-window decoding of motor imagery EEG
% ================================================================

clear; clc; close all;

%% -------------------- Load or Simulate Data --------------------
% X: [n_channels x n_samples x n_trials]
% y: [n_trials x 1], e.g. 1 = left hand, 2 = right hand
% fs: sampling frequency (Hz)

load('EEGdata.mat'); % <-- replace with your own file (must contain X, y, fs)

X  = yourEEG;   % [channels x samples x trials]
y  = yourLabels; % [trials x 1], integers (e.g., 1,2)
fs = 1000;        % or whatever your sampling rate is
save('EEGdata.mat', 'X', 'y', 'fs');

[nCh, nSamp, nTrials] = size(X);
fprintf('Loaded data: %d channels, %d samples, %d trials\n', nCh, nSamp, nTrials);

%% -------------------- Band-pass Filter --------------------
filt_band = [8 30]; % mu + beta band
[b, a] = butter(4, filt_band/(fs/2));

for i = 1:nTrials
    X(:,:,i) = filtfilt(b, a, double(X(:,:,i))')';
end

%% -------------------- CSP Function --------------------
function [W] = computeCSP(X, y)
    class1 = X(:,:,y==1);
    class2 = X(:,:,y==2);

    cov1 = zeros(size(X,1));
    cov2 = zeros(size(X,1));

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

for t = 1:length(time_idx)
    idx = time_idx(t):(time_idx(t)+win_samp-1);
    Xwin = X(:, idx, :);

    % Compute CSP filters for this time window
    W = computeCSP(Xwin, y);
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
end

%% -------------------- Plot Decoding Accuracy --------------------
time_axis = time_idx / fs;
figure;
plot(time_axis, acc_time * 100, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Accuracy (%)');
title('Time-Resolved Motor Imagery Decoding (CSP + SVM)');
grid on;
ylim([40 100]);
yline(50, '--r', 'Chance Level');

%% -------------------- Summary --------------------
[~, bestIdx] = max(acc_time);
fprintf('Peak decoding accuracy: %.2f%% at %.2f s\n', acc_time(bestIdx)*100, time_axis(bestIdx));














%% ================================
%% Version 2 with EEGlab process
% 1. Load and preprocess EEG data
% =================================

% Load EEGLAB
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

% Load your dataset (replace with your file)
EEG = pop_loadset('filename','mydata.set','filepath','C:\data\');

% Basic preprocessing
EEG = pop_eegfiltnew(EEG, 8, 30);        % bandpass filter 8–30 Hz (sensorimotor band)
EEG = pop_reref(EEG, []);                % common average reference
EEG = pop_runica(EEG, 'extended', 1);    % run ICA to clean artifacts
EEG = pop_subcomp(EEG, [1 3 7], 0);      % remove artifact components manually or automatically

% Extract data
data = EEG.data;       % [channels x samples]
fs   = EEG.srate;      % sampling rate

disp(['EEG data size: ', num2str(size(data))]);

%% ================================
% 2. Load or align continuous target
% =================================
% Target example: continuous variable (e.g. hand position or EMG)
load('target_signal.mat'); % should contain variable "target"

% Resample if needed to match EEG
if length(target) ~= size(data,2)
    target = resample(target, size(data,2), length(target));
end

% Normalize target
target = (target - mean(target)) / std(target);

%% ================================
% 3. Extract features over time
% =================================
% We'll use short-time log-bandpower features.

win  = 1;     % window length (seconds)
step = 0.1;   % step size (seconds)
nWin = floor((size(data,2)/fs - win)/step);
nChan = size(data,1);

features = zeros(nWin, nChan);
target_win = zeros(nWin,1);

for i = 1:nWin
    idx = round((i-1)*step*fs + 1) : round((i-1)*step*fs + win*fs);
    segment = data(:, idx);
    % compute bandpower for each channel
    bp = bandpower(segment', fs, [8 30]); % vector of bandpowers
    features(i,:) = log(bp + eps);
    target_win(i) = mean(target(idx));    % average target in same window
end

disp('Feature extraction complete.');

%% ================================
% 4. Split data into train/test
% =================================
nSamples = size(features,1);
idx = randperm(nSamples);
train_idx = idx(1:round(0.8*nSamples));
test_idx  = idx(round(0.8*nSamples)+1:end);

X_train = features(train_idx,:);
y_train = target_win(train_idx);
X_test  = features(test_idx,:);
y_test  = target_win(test_idx);

%% ================================
% 5. Train SVR decoder
% =================================
% Use Support Vector Regression to predict the continuous target
mdl = fitrsvm(X_train, y_train, ...
    'KernelFunction', 'rbf', ...
    'KernelScale', 'auto', ...
    'BoxConstraint', 1, ...
    'Epsilon', 0.1, ...
    'Standardize', true);

% Predict on test data
y_pred = predict(mdl, X_test);

% Evaluate
R = corr(y_pred, y_test);
RMSE = sqrt(mean((y_pred - y_test).^2));

fprintf('Decoding correlation: %.3f\n', R);
fprintf('RMSE: %.3f\n', RMSE);

%% ================================
% 6. Visualize results
% =================================
figure;
subplot(2,1,1);
plot(y_test, 'k'); hold on;
plot(y_pred, 'r');
xlabel('Time windows');
ylabel('Predicted value');
legend('True','Predicted');
title(sprintf('Continuous EEG Decoding (R = %.2f)', R));

subplot(2,1,2);
scatter(y_test, y_pred, 20, 'filled');
xlabel('True target'); ylabel('Predicted target');
title('Correlation plot');
grid on;











%% =================================
%% V3 with CSP elements
% 1. Load and preprocess EEG data
% =================================
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
EEG = pop_loadset('filename','mydata.set','filepath','C:\data\');

EEG = pop_eegfiltnew(EEG, 8, 30); % motor band
EEG = pop_reref(EEG, []);

data = EEG.data;      % [channels x samples]
fs   = EEG.srate;

load('target_signal.mat');  % target: continuous behavioral variable
if length(target) ~= size(data,2)
    target = resample(target, size(data,2), length(target));
end
target = (target - mean(target)) / std(target);


%% =================================
% 2. Create windowed epochs
% =================================
win  = 1; step = 0.1;
nWin = floor((size(data,2)/fs - win)/step);
nChan = size(data,1);

X = zeros(nChan, win*fs, nWin);
y = zeros(nWin, 1);

for i = 1:nWin
    idx = round((i-1)*step*fs + 1) : round((i-1)*step*fs + win*fs);
    X(:,:,i) = data(:, idx);
    y(i) = mean(target(idx));
end


%% =================================
% 3. Split into "low" vs "high" states for CSP
% =================================
th_low = prctile(y, 33);
th_high = prctile(y, 67);
idx_low = find(y <= th_low);
idx_high = find(y >= th_high);

X_low = X(:,:,idx_low);
X_high = X(:,:,idx_high);


%% =================================
% 4. Compute CSP filters
% =================================
C1 = zeros(nChan, nChan);
C2 = zeros(nChan, nChan);

for i = 1:size(X_low,3)
    covMat = cov(X_low(:,:,i)');
    C1 = C1 + covMat / trace(covMat);
end
C1 = C1 / size(X_low,3);

for i = 1:size(X_high,3)
    covMat = cov(X_high(:,:,i)');
    C2 = C2 + covMat / trace(covMat);
end
C2 = C2 / size(X_high,3);

% Generalized eigenvalue problem
[EVecs, EVals] = eig(C1, C1 + C2);
[~, ind] = sort(diag(EVals), 'descend');
W = EVecs(:, ind); % spatial filters

% Select top and bottom filters
nCSP = 3;
Wf = [W(:, 1:nCSP), W(:, end-nCSP+1:end)];


%% =================================
% 5. Apply CSP filters to all data and extract features
% =================================
nFeat = size(Wf,2);
features = zeros(nWin, nFeat);

for i = 1:nWin
    Xf = Wf' * X(:,:,i);
    features(i,:) = log(var(Xf, 0, 2)); % log-variance of CSP components
end


%% =================================
% 6. Continuous decoding with SVR
% =================================
nSamples = size(features,1);
idx = randperm(nSamples);
train_idx = idx(1:round(0.8*nSamples));
test_idx  = idx(round(0.8*nSamples)+1:end);

X_train = features(train_idx,:);
y_train = y(train_idx);
X_test  = features(test_idx,:);
y_test  = y(test_idx);

mdl = fitrsvm(X_train, y_train, ...
    'KernelFunction', 'rbf', ...
    'KernelScale', 'auto', ...
    'BoxConstraint', 1, ...
    'Epsilon', 0.1, ...
    'Standardize', true);

y_pred = predict(mdl, X_test);

R = corr(y_pred, y_test);
RMSE = sqrt(mean((y_pred - y_test).^2));

fprintf('CSP-based decoding: R = %.3f, RMSE = %.3f\n', R, RMSE);


%% =================================
% 7. Visualize results
% =================================
figure;
plot(y_test, 'k'); hold on;
plot(y_pred, 'r');
legend('True', 'Predicted');
xlabel('Time windows');
ylabel('Target');
title(sprintf('Continuous EEG Decoding using CSP + SVR (R = %.2f)', R));




%% Filterbank CSP variant
%% =================================
% 1. Load and preprocess EEG
% =================================
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
EEG = pop_loadset('filename','mydata.set','filepath','C:\data\');

data = EEG.data;    % [channels x samples]
fs   = EEG.srate;

load('target_signal.mat'); % target: continuous variable
if length(target) ~= size(data,2)
    target = resample(target, size(data,2), length(target));
end
target = (target - mean(target)) / std(target);

%% =================================
% 2. Define filter bank
% =================================
bands = [8 12; 12 16; 16 20; 20 30]; % μ + β bands
nBands = size(bands,1);

%% =================================
% 3. Window EEG and target
% =================================
win  = 1; step = 0.1;
nWin = floor((size(data,2)/fs - win)/step);
nChan = size(data,1);

X = zeros(nChan, win*fs, nWin);
y = zeros(nWin,1);

for i = 1:nWin
    idx = round((i-1)*step*fs + 1) : round((i-1)*step*fs + win*fs);
    X(:,:,i) = data(:, idx);
    y(i) = mean(target(idx));
end

%% =================================
% 4. Compute CSP features for each band
% =================================
nCSP = 3; % number of CSP filters per end
features = [];

for b = 1:nBands
    % Bandpass filter EEG
    [bFilt,aFilt] = butter(4, bands(b,:)/(fs/2));
    Xfilt = zeros(size(X));
    for iWin = 1:nWin
        Xfilt(:,:,iWin) = filtfilt(bFilt,aFilt, double(X(:,:,iWin))')';
    end
    
    % Discretize target for CSP: low vs high
    th_low = prctile(y,33); th_high = prctile(y,67);
    idx_low = find(y <= th_low);
    idx_high = find(y >= th_high);
    X_low  = Xfilt(:,:,idx_low);
    X_high = Xfilt(:,:,idx_high);
    
    % Compute CSP
    C1 = zeros(nChan); C2 = zeros(nChan);
    for i=1:size(X_low,3)
        C1 = C1 + cov(X_low(:,:,i)')/trace(cov(X_low(:,:,i)'));
    end
    C1 = C1 / size(X_low,3);
    for i=1:size(X_high,3)
        C2 = C2 + cov(X_high(:,:,i)')/trace(cov(X_high(:,:,i)'));
    end
    C2 = C2 / size(X_high,3);
    
    [EVecs, EVals] = eig(C1, C1+C2);
    [~, ind] = sort(diag(EVals),'descend');
    W = EVecs(:,ind);
    Wcsp = [W(:,1:nCSP), W(:,end-nCSP+1:end)];
    
    % Extract log-variance features for all windows
    feat_band = zeros(nWin, 2*nCSP);
    for iWin = 1:nWin
        Xcsp = Wcsp' * Xfilt(:,:,iWin);
        feat_band(iWin,:) = log(var(Xcsp,0,2));
    end
    
    % Concatenate features across bands
    features = [features feat_band];
end

disp(['Total feature size: ' num2str(size(features))]);

%% =================================
% 5. Train SVR for continuous decoding
% =================================
nSamples = size(features,1);
idx = randperm(nSamples);
train_idx = idx(1:round(0.8*nSamples));
test_idx  = idx(round(0.8*nSamples)+1:end);

X_train = features(train_idx,:);
y_train = y(train_idx);
X_test  = features(test_idx,:);
y_test  = y(test_idx);

mdl = fitrsvm(X_train, y_train, ...
    'KernelFunction','rbf', ...
    'KernelScale','auto', ...
    'BoxConstraint',1, ...
    'Epsilon',0.1, ...
    'Standardize',true);

y_pred = predict(mdl, X_test);

R = corr(y_pred, y_test);
RMSE = sqrt(mean((y_pred - y_test).^2));

fprintf('FBCSP decoding: R = %.3f, RMSE = %.3f\n', R, RMSE);

%% =================================
% 6. Visualize predictions
% =================================
figure;
plot(y_test,'k'); hold on; plot(y_pred,'r');
xlabel('Time windows'); ylabel('Target');
legend('True','Predicted');
title(sprintf('FBCSP + SVR Continuous EEG Decoding (R=%.2f)', R));