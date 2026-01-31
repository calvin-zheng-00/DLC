function res = decode_angles_trajectory_multioutput_ridge(EEG, Y, fsEEG, fsY, opts)
% Decode multiple joint-angle trajectories per trial using ONE multi-output ridge model.
%
% EEG: [nChan x nSamp x nTrials] @ fsEEG
% Y  : angles, variable format (see below), radians
%      Preferred: [nY x nAngles x nTrials]
%      Also accepted:
%        [nY x nTrials x nAngles]
%        [nTrials x nY x nAngles]
%        [nY x nTrials]  (single angle)
%
% Output:
%   res.Ypred: [nY x nAngles x nTrials]
%   res.Ytrue: [nY x nAngles x nTrials]
%   res.tY   : [nY x 1]
%   res.trialR2_perAngle: [nTrials x nAngles]
%   res.trialR2_mean: [nTrials x 1]

%% defaults
if nargin < 5, opts = struct(); end
if ~isfield(opts,'epochDur'),      opts.epochDur = 3.0; end
if ~isfield(opts,'winLenEEG'),     opts.winLenEEG = 0.25; end
if ~isfield(opts,'bands'),         opts.bands = [8 12; 13 20; 20 30; 30 40]; end
if ~isfield(opts,'lagsY'),         opts.lagsY = [-2 -1 0 1 2]; end
if ~isfield(opts,'lambda'),        opts.lambda = 300; end
if ~isfield(opts,'kfold'),         opts.kfold = 5; end
if ~isfield(opts,'standardizeX'),  opts.standardizeX = true; end
if ~isfield(opts,'exampleTrial'),  opts.exampleTrial = 1; end
if ~isfield(opts,'exampleAngles'), opts.exampleAngles = []; end % indices to plot; [] -> up to 3
if ~isfield(opts,'plotTrials'),        opts.plotTrials = []; end   % e.g., [1 3 7 12]
if ~isfield(opts,'maxPlotTrials'),     opts.maxPlotTrials = 10; end
if ~isfield(opts,'plotMode'),          opts.plotMode = "overlay"; end % "overlay" or "grid"

assert(ndims(EEG)==3, 'EEG must be [nChan x nSamp x nTrials].');
[~, nSamp, nTr] = size(EEG);

% --- Coerce Y into [nY x nAngles x nTrials] ---
[Y3, nY, nAngles, nTrY] = coerce_Y(Y, nTr);
assert(nTrY == nTr, 'Y trials must match EEG trials.');

tY = (0:nY-1)'/fsY;

% sanity checks
durEEG = nSamp/fsEEG;
durY = nY/fsY;
if abs(durEEG - opts.epochDur) > 0.1
    warning('EEG duration %.3fs differs from opts.epochDur %.3fs.', durEEG, opts.epochDur);
end
if abs(durY - opts.epochDur) > 0.1
    warning('Y duration %.3fs differs from opts.epochDur %.3fs.', durY, opts.epochDur);
end

%% Build per-trial design Xtrial
Fseq = eeg_features_sequence_aligned(EEG, fsEEG, tY, opts.winLenEEG, opts.bands); % [nY x nFeat x nTr]
Flag = add_lags_sequence(Fseq, opts.lagsY);                                       % [nY x nFeatLag x nTr]

[nY2, nFeatLag, ~] = size(Flag);
assert(nY2 == nY);
P = nY * nFeatLag;

Xtrial = zeros(nTr, P);
for tr = 1:nTr
    Xtrial(tr,:) = reshape(Flag(:,:,tr).', 1, []); % time-major flatten
end

%% Build multi-angle target: YtrialMat [nTr x (nY*nAngles)]
% Y3: [nY x nAngles x nTr]
YtrialMat = zeros(nTr, nY*nAngles);
for tr = 1:nTr
    Ymat = Y3(:,:,tr);              % [nY x nAngles]
    YtrialMat(tr,:) = reshape(Ymat, 1, []); % time-major: [1 x (nY*nAngles)]
end

%% Cross-validated multi-output ridge
cvp = cvpartition(nTr, 'KFold', opts.kfold);
YpredMat = nan(nTr, nY*nAngles);

for k = 1:opts.kfold
    trIdx = training(cvp,k);
    teIdx = test(cvp,k);

    Xtr = Xtrial(trIdx,:);
    Xte = Xtrial(teIdx,:);
    Ytr = YtrialMat(trIdx,:);

    % ---- Standardize X using training only ----
    if opts.standardizeX
        mu = mean(Xtr,1);
        sd = std(Xtr,0,1); sd(sd<1e-12)=1;
        XtrZ = (Xtr - mu)./sd;
        XteZ = (Xte - mu)./sd;
    else
        XtrZ = Xtr;
        XteZ = Xte;
    end

    % W = (XtrZ.'*XtrZ + opts.lambda*eye(size(XtrZ,2))) \ (XtrZ.'*Ytr);
    % b = mean(Ytr,1) - mean(XtrZ,1)*W;
    % YpredMat(teIdx,:) = XteZ*W + b;
    
    % ---- PCA on TRAIN only (no leakage) ----
    % Choose number of PCs:
    %   - opts.nPC can be a fixed number (e.g., 10)
    %   - or opts.varToKeep = 95 to keep 95% variance (recommended)
    [coeff, scoreTr, ~, ~, explained, muP] = pca(XtrZ, 'Centered', true);

    if isfield(opts,'varToKeep') && ~isempty(opts.varToKeep)
        cumExp = cumsum(explained);
        nPC = find(cumExp >= opts.varToKeep, 1, 'first');
    else
        nPC = opts.nPC;
    end

    % cap to what's possible given nTrain
    nPC = min(nPC, size(scoreTr,2));   % scoreTr has at most nTrain-1 cols

    XtrP = scoreTr(:,1:nPC);
    XteP = (XteZ - muP) * coeff(:,1:nPC);

    % ---- Multi-output ridge on reduced features ----
    W = (XtrP.'*XtrP + opts.lambda*eye(nPC)) \ (XtrP.'*Ytr);
    b = mean(Ytr,1) - mean(XtrP,1)*W;

    YpredMat(teIdx,:) = XteP*W + b;
end

%% Reshape predictions back to [nY x nAngles x nTr]
Ypred = zeros(nY, nAngles, nTr);
for tr = 1:nTr
    vec = YpredMat(tr,:);                 % [1 x (nY*nAngles)]
    Ypred(:,:,tr) = reshape(vec, [nY, nAngles]);
end

Ytrue = Y3;

%% Plots: mean predicted vs mean true for selected angles
if isempty(opts.exampleAngles)
    exampleAngles = 1:min(3, nAngles);
else
    exampleAngles = opts.exampleAngles(:)';
    exampleAngles = exampleAngles(exampleAngles>=1 & exampleAngles<=nAngles);
    if isempty(exampleAngles), exampleAngles = 1:min(3,nAngles); end
end

YtrueMean = mean(Ytrue, 3, 'omitnan'); % [nY x nAngles]
YpredMean = mean(Ypred, 3, 'omitnan'); % [nY x nAngles]

figure('Color','w');
for a = exampleAngles
    plot(tY, YtrueMean(:,a), 'LineWidth', 2); hold on;
    plot(tY, YpredMean(:,a), 'LineWidth', 2);
end
grid on;
xlabel('Time (s)'); ylabel('Angle (rad)');
title('Whole-trajectory decoder: mean predicted vs mean true (selected angles)');
leg = {};
for a = exampleAngles
    leg{end+1} = sprintf('True ang %d (mean)', a); %#ok<AGROW>
    leg{end+1} = sprintf('Pred ang %d (mean)', a); %#ok<AGROW>
end
legend(leg, 'Location','eastoutside');

% Plot multiple trials overlaid (true and predicted), for selected angles
if isempty(opts.plotTrials)
    % default: pick up to maxPlotTrials evenly spaced trials
    nToPlot = min(opts.maxPlotTrials, nTr);
    if nToPlot <= 1
        trialList = 1;
    else
        trialList = unique(round(linspace(1, nTr, nToPlot)));
    end
else
    trialList = unique(opts.plotTrials(:)');
    trialList = trialList(trialList>=1 & trialList<=nTr);
    if isempty(trialList)
        trialList = 1:min(opts.maxPlotTrials, nTr);
    end
end

figure('Color','w');
for a = exampleAngles
    subplot(numel(exampleAngles), 1, find(exampleAngles==a,1,'first'));
    hold on;

    % Overlay TRUE trials
    for tr = trialList
        plot(tY, Ytrue(:,a,tr), 'LineWidth', 1);
    end

    % Overlay PRED trials (dashed)
    for tr = trialList
        plot(tY, Ypred(:,a,tr), '--', 'LineWidth', 1);
    end

    grid on;
    xlabel('Time (s)');
    ylabel('Angle (rad)');
    title(sprintf('Angle %d: overlay %d trials (true=solid, pred=dashed)', a, numel(trialList)));
end


% % Example trial overlay
% ex = max(1, min(nTr, opts.exampleTrial));
% figure('Color','w');
% for a = exampleAngles
%     plot(tY, Ytrue(:,a,ex), 'LineWidth', 2); hold on;
%     plot(tY, Ypred(:,a,ex), 'LineWidth', 2);
% end
% grid on;
% xlabel('Time (s)'); ylabel('Angle (rad)');
% title(sprintf('Whole-trajectory decoder: trial %d predicted vs true (selected angles)', ex));
% legend(leg, 'Location','best');

%% Trial-wise R^2 per angle
trialR2_perAngle = nan(nTr, nAngles);
for tr = 1:nTr
    for a = 1:nAngles
        yt = Ytrue(:,a,tr);
        yp = Ypred(:,a,tr);
        ss_res = sum((yt-yp).^2);
        ss_tot = sum((yt-mean(yt)).^2) + 1e-12;
        trialR2_perAngle(tr,a) = 1 - ss_res/ss_tot;
    end
end
trialR2_mean = mean(trialR2_perAngle, 2, 'omitnan');

figure('Color','w');
histogram(trialR2_mean, 15); grid on;
xlabel('Trial-wise mean R^2 across angles'); ylabel('Count');
title('Whole-trajectory decoder: distribution of trial-wise mean R^2');

%% Outputs
res = struct();
res.tY = tY;
res.Ytrue = Ytrue;
res.Ypred = Ypred;
res.YtrueMean = YtrueMean;
res.YpredMean = YpredMean;
res.trialR2_perAngle = trialR2_perAngle;
res.trialR2_mean = trialR2_mean;
res.opts = opts;

end

%% ===================== helpers =====================

function [Y3, nY, nAngles, nTr] = coerce_Y(Y, nTrExpected)
% Convert various Y formats into [nY x nAngles x nTrials]
if ismatrix(Y)
    % [nY x nTrials] single angle
    if size(Y,2) == nTrExpected
        Y3 = reshape(Y, [size(Y,1), 1, size(Y,2)]);
    elseif size(Y,1) == nTrExpected
        Y = Y.'; % now [nY x nTrials]
        Y3 = reshape(Y, [size(Y,1), 1, size(Y,2)]);
    else
        error('2D Y must be [nY x nTrials] or [nTrials x nY].');
    end
elseif ndims(Y) == 3
    sz = size(Y);
    % try [nY x nAngles x nTrials]
    if sz(3) == nTrExpected
        Y3 = Y;
    % try [nY x nTrials x nAngles]
    elseif sz(2) == nTrExpected
        Y3 = permute(Y, [1 3 2]);
    % try [nTrials x nY x nAngles]
    elseif sz(1) == nTrExpected
        Y3 = permute(Y, [2 3 1]);
    else
        error('3D Y must contain nTrials matching EEG in one dimension.');
    end
else
    error('Y must be 2D or 3D.');
end

nY = size(Y3,1);
nAngles = size(Y3,2);
nTr = size(Y3,3);
end

function Fseq = eeg_features_sequence_aligned(EEG, fsEEG, tY, winLenSec, bands)
% Fseq: [nY x nFeat x nTr], nFeat = nChan*nBands
[nChan, nSamp, nTr] = size(EEG);
nY = numel(tY);
nBands = size(bands,1);
nFeat = nChan*nBands;

winSamp = max(16, round(winLenSec*fsEEG));
half = floor(winSamp/2);

Fseq = zeros(nY, nFeat, nTr);

for ti = 1:nY
    c = round(tY(ti)*fsEEG) + 1;
    c = min(max(c,1), nSamp);

    i0 = max(1, c - half);
    i1 = min(nSamp, i0 + winSamp - 1);
    i0 = max(1, i1 - winSamp + 1);
    idx = i0:i1;

    Xwin = EEG(:,idx,:);                  % [chan x win x tr]
    F = bandpower_features(Xwin, fsEEG, bands); % [tr x nFeat]
    Fseq(ti,:,:) = permute(F, [3 2 1]);   % [1 x nFeat x tr]
end
end

function Flag = add_lags_sequence(Fseq, lagsY)
% Fseq: [nY x nFeat x nTr]
% Flag: [nY x (nFeat*nLags) x nTr]
[nY, nFeat, nTr] = size(Fseq);
nL = numel(lagsY);
Flag = zeros(nY, nFeat*nL, nTr);

for li = 1:nL
    L = lagsY(li);
    src = (1:nY) - L;
    src(src < 1) = 1;
    src(src > nY) = nY;
    block = Fseq(src,:,:); % [nY x nFeat x nTr]
    c0 = (li-1)*nFeat + 1;
    c1 = li*nFeat;
    Flag(:,c0:c1,:) = block;
end
end

function F = bandpower_features(Xwin, fs, bands)
% Xwin: [chan x samp x tr] -> F: [tr x (chan*bands)] log bandpower
[nChan, ~, nTr] = size(Xwin);
nBands = size(bands,1);
F = zeros(nTr, nChan*nBands);

for tr = 1:nTr
    pos = 1;
    for ch = 1:nChan
        x = double(squeeze(Xwin(ch,:,tr)));
        [Pxx,Fhz] = pwelch(x, [], [], [], fs);
        for b = 1:nBands
            fidx = (Fhz >= bands(b,1) & Fhz <= bands(b,2));
            bp = mean(Pxx(fidx));
            F(tr,pos) = log(bp + 1e-12);
            pos = pos + 1;
        end
    end
end
end


opts = struct();
opts.epochDur = 3.0;
opts.winLenEEG = 0.25;
opts.bands = [8 12; 13 20; 20 30; 30 40];
opts.lagsY = [-4 -2 0 2 4];   % +/- 100 ms
opts.lambda = 500;            % usually larger for many outputs
opts.kfold = 5;
opts.varToKeep = 95;    % keep 95% of variance
opts.exampleAngles = 4;
opts.plotTrials    = 1:30; % which trials to overlay (optional)
% opts.maxPlotTrials = 10;          % used only if plotTrials is empty

EEG = ALLEEG(2).data;
% y: [nY x nJoints x nTrials] (radians)
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

src = 'C:\Users\czhe0008\Documents\EEG\3d\11_12_angles';
Y = [];

trialtemp = dir(fullfile(src,'*.csv'));
trialfolder = {trialtemp(~[trialtemp.isdir]).name};
for trialfolder_i = 1:numel(trialfolder)
    data = fullfile(src,trialfolder{trialfolder_i});
    T = table2array(readtable(data));
    temp = T((event_time(trialfolder_i)-39):(event_time(trialfolder_i)+80),:);
    Y = cat(3,Y,temp);
end

src_score = 'C:\Users\czhe0008\Documents\EEG\PCA\11_12\Projection\Projection.csv';
Y2 = [];
T = table2array(readtable(src_score));
for i = 1:30
    temp = T((((i-1)*120)+1):(i*120),1:8); % Adjust how many PCs to keep
    Y2 = cat(3,Y2,temp);
end

% res = decode_angles_trajectory_multioutput_ridge(EEG, Y, 1000, 40, opts);
% res2 = decode_angles_trajectory_multioutput_ridge(EEG, Y2, 1000, 40, opts);


score = res2.Ypred;
temp = res.Ypred;
baseline = mean(Y,3);
coeff = table2array(readtable('C:\Users\czhe0008\Documents\EEG\PCA\11_12\PCA_coeffs\coeffs.csv'));
mu = table2array(readtable('C:\Users\czhe0008\Documents\EEG\PCA\11_12\PCA_mean\mean.csv'));
Xstd = table2array(readtable('C:\Users\czhe0008\Documents\EEG\PCA\11_12\All\sigma.csv'));
Xmean = table2array(readtable('C:\Users\czhe0008\Documents\EEG\PCA\11_12\All\mu_global.csv'));
for i = 1:30
    Xz_recon = score(:,:,i) * coeff(2:end,1:8)' + mu;
    X_recon  = Xz_recon .* Xstd + Xmean;
    figure(2)
    hold on
    % plot(linspace(-1,2,120),X_recon(:,5),'Color','r')
    % plot(linspace(-1,2,120),temp(:,5,i),'Color','k')
    plot(linspace(-1,2,120),Y(:,5,i),'Color','b')
    plot(linspace(-1,2,120),baseline(:,5),'Color','g')
    % legend('Synergy Prediction','Joint Prediction','True Angle',Location='southwest')
    legend('True Angle','True Angle Mean',Location='southwest')
    title('Decoded Joint Angles Over Time')
    % title(strcat('Decoded Joint Angles Over Time Trial ',int2str(i)))
    xlabel('Time (s)')
    ylabel('Joint angle (Rad)')
    % for j = 1:15
    %     plot(linspace(-1,2,120),X_recon(:,j),'Color','r','HandleVisibility','off')
    %     plot(linspace(-1,2,120),temp(:,j,i),'Color','k','HandleVisibility','off')
    %     plot(linspace(-1,2,120),Y(:,j,i),'Color','b','HandleVisibility','off')
    % end
    % saveas(gcf,strcat('C:\Users\czhe0008\Documents\EEG\Figures\11_12\synergy\trial',int2str(i),'.png'))
    % close(2)
end