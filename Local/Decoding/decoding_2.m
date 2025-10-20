%% --- DATA ---
% src = 'C:\Users\czhe0008\Documents\EEG\3d\13_10_angle_trial\20251013T121715-122433';
% Y_raw = [];
% maintemp = dir(fullfile(src,'*.csv'));
% mainfolder = {maintemp(~[maintemp.isdir]).name};
% for mainfolder_i = 1:numel(mainfolder)
%     df = table2array(readtable(fullfile(src,mainfolder{mainfolder_i})));
%     Y_raw = cat(3,Y_raw,df);
% end
% Y = resample(Y_raw, 1000, 40);
% EEG_data = ALLEEG(6).data;

%% --- PCA ---
% PC_coeff = 'C:\Users\czhe0008\Documents\EEG\3d\13_10_PCA\PCA_coeffs\20251013T121715-122433.csv';
% PC_mean = 'C:\Users\czhe0008\Documents\EEG\3d\13_10_PCA\PCA_mean\20251013T121715-122433.csv';
% coeff = table2array(readtable(PC_coeff));
% mu = table2array(readtable(PC_mean));
% Y_centred = Y - mu;
% Y = Y_centred * coeff;

Projected = 'C:\Users\czhe0008\Documents\EEG\3d\13_10_PCA\Projection\20251013T134540-135116.csv';
df = table2array(readtable(Projected));
Y_raw = [];
for i = 1:(size(df,1)/120)
    Y_raw = cat(3,Y_raw,df((120*(i-1)+1):120*i,:));
end
Y = resample(Y_raw, 1000, 40);
EEG_data = ALLEEG(9).data;

%% --- PARAMETERS ---
Fs = 1000;                 % EEG sampling rate (Hz)
window_len = 0.5;          % window length (s)
step = 0.2;                % step size (s)
nTargets = size(Y,2);      % number of joint angles (targets)
nTrials = size(EEG_data,3);

% Define filterbank bands [low high] Hz
bands = [ 1 4;
          4 8;
          8 12;
          13 30;
          30 45 ];

%% --- DESIGN FILTERS ---
for b = 1:size(bands,1)
    fb(b).d = designfilt('bandpassiir', 'FilterOrder', 4, ...
        'HalfPowerFrequency1', bands(b,1), ...
        'HalfPowerFrequency2', bands(b,2), ...
        'SampleRate', Fs);
end

%% --- INITIALIZE FEATURE/ TARGET MATRICES ---
features_all = [];
targets_all  = [];

win_samp = round(window_len * Fs);
step_samp = round(step * Fs);

%% --- LOOP OVER TRIALS ---
for tr = 1:nTrials
    eeg_trial = EEG_data(:,:,tr);   % [channels × samples]
    y_trial   = Y(:,:,tr);          % [samples × nTargets]

    % Number of windows for this trial
    nWins = floor((size(eeg_trial,2) - win_samp) / step_samp);

    for i = 1:nWins
        idx = (1:win_samp) + (i-1)*step_samp;
        feats_per_band = [];

        % --- FILTERBANK FEATURE EXTRACTION ---
        for b = 1:length(fb)
            Xf = filtfilt(fb(b).d, eeg_trial')';   % filter each band
            Xwin = Xf(:, idx);

            % Simple log-variance features
            feat = log(var(Xwin,0,2));  % [channels × 1]
            feats_per_band = [feats_per_band; feat];
        end

        features_all = [features_all; feats_per_band'];
        targets_all  = [targets_all; mean(y_trial(idx,:),1)];
    end
end

%% --- TRAIN/TEST SPLIT (trial-wise or random) ---
nSamples = size(features_all,1);
nTrain = round(0.8 * nSamples);
Xtrain = features_all(1:nTrain,:);
Xtest  = features_all(nTrain+1:end,:);
Ytrain = targets_all(1:nTrain,:);
Ytest  = targets_all(nTrain+1:end,:);

%% --- TRAIN MULTI-OUTPUT SVR ---
nTargets = size(Ytrain,2);
mdl = cell(1,nTargets);
Ypred = zeros(size(Ytest));

for j = 1:nTargets
    mdl{j} = fitrsvm(Xtrain, Ytrain(:,j), ...
                     'KernelFunction','rbf', ...
                     'Standardize',true, ...
                     'KernelScale','auto');
    Ypred(:,j) = predict(mdl{j}, Xtest);
end

%% --- EVALUATION ---
R  = zeros(1,nTargets);
R2 = zeros(1,nTargets);

for j = 1:nTargets
    R(j)  = corr(Ypred(:,j), Ytest(:,j));
    R2(j) = R(j)^2;
    fprintf('Joint %d: R = %.3f, R² = %.3f\n', j, R(j), R2(j));
end

%% --- PLOTTING ---
t = (1:size(Ytest,1))*step;
figure;
% joints = [2,3,4,5,7,8];
% labels = ["TMCP flex","TIP flex","IMCP flex","IPIP flex","MMCP flex","MPIP flex"];
joints = 1:6;
labels = ["PC1","PC2","PC3","PC4","PC5","PC6"];
for j = 1:6
    subplot(6,1,j);
    plot(t, Ytest(:,joints(j)), 'k', 'LineWidth',1.2); hold on;
    plot(t, Ypred(:,joints(j)), 'r');
    ylabel(labels(j));
    legend('True','Pred');
end
xlabel('Time (s)');
sgtitle('Filterbank EEG Decoding of Side Pinch PCs');
