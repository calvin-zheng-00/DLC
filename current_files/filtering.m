close all; clear all; clc;

% src = 'C:\Users\USER\OneDrive\Documents\MATLAB\data\p26_3d';
% dest = 'C:\Users\USER\OneDrive\Documents\MATLAB\filtered\p26_filter';
src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\eeg_3d';
dest = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\eeg_filter';
threshold = 30;   % When to start flagging jumps
increment = 20;  % How fast to expand acceptable range

maintemp = dir(fullfile(src,'*'));
mainfolder = setdiff({maintemp([maintemp.isdir]).name},{'.','..'});



for mainfolder_i = 1:numel(mainfolder)
    % Extra loop incase of categorised 3d data
    extratemp = dir(fullfile(src,mainfolder{mainfolder_i},'*'));
    extrafolder = setdiff({extratemp([extratemp.isdir]).name},{'.','..'});
    for extrafolder_k = 1:numel(extrafolder)
        trialtemp = dir(fullfile(src,mainfolder{mainfolder_i},extrafolder{extrafolder_k},'*.csv'));
        trialfolder = {trialtemp(~[trialtemp.isdir]).name};
        for trialfolder_j = 1:numel(trialfolder)
            if strcmp(trialfolder{trialfolder_j}, "df_3d.csv")
                data = fullfile(src,mainfolder{mainfolder_i},extrafolder{extrafolder_k},trialfolder{trialfolder_j});
                df = readtable(data);
                df([1:5,end-4:end],:) = [];
                df(:,[1:3,79:end]) = [];
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
                    
                    % [z, p, k] = butter(2,fc/(fs/2));
                    % [sos,g] = zp2sos(z,p,k);
                    % filtered = filtfilt(sos,g,col);

                    [b, a] = butter(2,fc/(fs/2));
                    % filtered = filtfilt(b,a,col,"ctf");
                    filtered = filtfilt(b,a,col);

                    % b, a = signal.butter(2, 5, fs=40)
                    % if len(col) < 3 * max(len(a), len(b))
                    %     filtered = signal.filtfilt(b, a, col, padlen=0)
                    % else
                    %     filtered = signal.filtfilt(b, a, col)
                    colname = df.Properties.VariableNames{joint_idx};
                    if df2_check == 1
                        df2 = table(filtered);
                        df2.Properties.VariableNames = convertCharsToStrings(colname);
                        df2_check = 0;
                    else
                        df2.(colname) = filtered;
                    end
                end
                writetable(df2,fullfile(dest,strcat(extrafolder{extrafolder_k},"_filtered.csv")))
            end
        end
    end
end