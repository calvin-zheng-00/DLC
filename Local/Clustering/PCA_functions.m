function PCA_filt(src,dest,threshold,increment)
    % Filter input data to remove bad trials.
    % src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\3d';
    % dest = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\filter';
    % threshold = 30;   % When to start flagging jumps
    % increment = 20;  % How fast to expand acceptable range
    
    parttemp = dir(fullfile(src,'*'));
    partfolder = setdiff({parttemp([parttemp.isdir]).name},{'.','..'});

    for participant_i = 1:numel(partfolder)
        % Extra loop incase of categorised 3d data
        acttemp = dir(fullfile(src,partfolder{participant_i},'*'));
        actfolder = setdiff({acttemp([acttemp.isdir]).name},{'.','..'});
        for activity_k = 1:numel(actfolder)
            trialtemp = dir(fullfile(src,partfolder{participant_i},actfolder{activity_k},'*'));
            trialfolder = setdiff({trialtemp([trialtemp.isdir]).name},{'.','..'});
            for trial_j = 1:numel(trialfolder)
                if isfile(fullfile(src,partfolder{participant_i},actfolder{activity_k},trialfolder{trial_j},'df_3d.csv')) && isfile(fullfile(src,partfolder{participant_i},actfolder{activity_k},trialfolder{trial_j},'df_err.csv'))
                    % Getting 3d position values
                    data = fullfile(src,partfolder{participant_i},actfolder{activity_k},trialfolder{trial_j},'df_3d.csv');
                    df = readtable(data);
                    % Getting error data
                    % data_err = fullfile(src,partfolder{participant_i},actfolder{activity_k},trialfolder{trial_j},'df_err.csv');
                    % df_err = readtable(data_err);
                    % Removing the timestep ends since those sections have
                    % bad tracking.
                    df([1:5,end-4:end],:) = [];
                    % df_err([1:5,end-4:end],:) = [];
                    % Removing unused markers
                    df(:,[1:3,79:end]) = [];
                    % df_err = df_err(:,2:26);
                    df2_check = 1;
                    for joint_idx = 1:width(df)
                        col = df.(joint_idx);
                        % diff = threshold;
                        % check = 1;
                        col = fillmissing(col,'linear',1,'EndValues','nearest');
                        % % Checking for large jumps
                        % for i = 2:length(col)
                        %     if abs(col(i) - col(check)) > diff
                        %         diff = diff + increment;
                        %     else
                        %         grad = (col(i)-col(check))/(i-check);
                        %         for j = check+1:i
                        %             col(j) = col(j-1) + grad;
                        %         end
                        %         diff = threshold;
                        %         check = i;
                        %     end
                        % end
                        % Butterworth
                        fc = 2;
                        fs = 40;    
                        [b, a] = butter(2,fc/(fs/2));
                        filtered = filtfilt(b,a,col);

                        colname = df.Properties.VariableNames{joint_idx};
                        if df2_check == 1
                            df2 = table(filtered);
                            df2.Properties.VariableNames = convertCharsToStrings(colname);
                            df2_check = 0;
                        else
                            df2.(colname) = filtered;
                        end
                    end
                    if ~isfolder(fullfile(dest,partfolder{participant_i},actfolder{activity_k}))
                        mkdir(fullfile(dest,partfolder{participant_i},actfolder{activity_k}));
                    end
                    writetable(df2,fullfile(dest,partfolder{participant_i},actfolder{activity_k},strcat(actfolder{activity_k},trialfolder{trial_j},"_filtered.csv")))
                end
            end
        end
    end
end

function PCA_participant()%src,src_3d)
    % Finding PCs from individuals.
    % category = "All";
    src = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\angle';
    src_3d = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\filter';
    dest = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Participant_PC';
    
    allAngles = [];
    % grasp_angles = [];
    % p_sig_80 = [];
    % p_latent = [];
    % allExplained = [];
    % allCoeff = [];
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
                % thisTable = thisTable(:,1:21); % Removing arm
                % thisTable = thisTable(:,1:15); % Removing arm and keeping only flexion
                % thisTable = thisTable(:,[1:5,7,8,10,11,13,14]); % Removing arm and keeping only flexion and removing DIP
                % thisTable = thisTable(:,[1:5,7,8,10,11,13,14,16:21]); % Removing arm and removing DIP
        
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
        [coeff,score,latent,tsquared,explained,mu] = pca(zAngles, 'Rows', 'pairwise'); %pAngles
        % p_latent = [p_latent,sum(explained(find(latent>1)))];
        % p_sig_80 = [p_sig_80,find(cumsum(explained)>80,1)]; % Finds the amount of PCs needed to explain 80% of variance
        % curr_sig_80 = find(cumsum(explained)>80,1);
        % allExplained = [allExplained,explained];
        % allCoeff(:,:,i) = coeff;
        % save(strcat('Significance\',mainfolder{i},'_sig_80.mat'),'curr_sig_80');

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

function PCA_act()
    % Finding PCs from individuals.
    src = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\angle';
    src_3d = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\filter';
    dest = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\Activity_Hand_PC';
    
    % allAngles = [];
    % allExplained = [];
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
                thisTable = thisTable(:,1:21); % Removing arm
                % thisTable = thisTable(:,[1:5,7,8,10,11,13,14,16:21]); % Removing arm and removing DIP
        
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
                end_temp = min(reach+10,height(thisTable));                          % Try 10 frames after reach end rather than 40
                end_time = [end_time;reach_time(end)+min(10,height(thisTable)-reach)];
        
                %Normalizing angles
                categoryAngles = [categoryAngles;thisTable(max([1,reach-39]):end_temp,:)];
                % video_data = readtable(fullfile(src_3d,mainfolder{i},actfolder{l},posfolder{j}));
                % video_data = video_data(max([1,reach-39]):end_temp,:);
                % video_maker(video_data)
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

function video_maker(data)
    T = data;
    for i = 0:((width(T)/3)-1)
        T{:,i*3 + 1} = T{:,i*3 + 1} - data{:,67};
        T{:,i*3 + 2} = T{:,i*3 + 2} - data{:,68};
        T{:,i*3 + 3} = T{:,i*3 + 3} - data{:,69};
    end
    
    
    h = height(T); % Number of rows in table
    startFrame = 1;
    
    ThumbPlot = plot3([T.TT_z(startFrame),T.TIP_z(startFrame),T.TMCP_z(startFrame),T.TCMC_z(startFrame),T.RW_z(startFrame)], ...
        [T.TT_x(startFrame),T.TIP_x(startFrame),T.TMCP_x(startFrame),T.TCMC_x(startFrame),T.RW_x(startFrame)], ...
        [T.TT_y(startFrame),T.TIP_y(startFrame),T.TMCP_y(startFrame),T.TCMC_y(startFrame),T.RW_y(startFrame)], ...
        '-o', 'MarkerSize',3,'MarkerFaceColor','r');
    axis([-700 300 -600 200 -500 500])
    
    % Stitch together different view points
    
    hold on
    xlabel('z (mm)')
    ylabel('y (mm)')
    zlabel('x (mm)')
    title('Activity')
    
    % Setting camera angle (and reversing the y in relation to the x and z values)
    set(gca, 'XDir','reverse')
    set(gca, 'ZDir','reverse')
    set(gca, 'fontsize',14)
    set(gcf,'position',[200,200,600,570])
    % Default is back view
    %view([90 10]) %front view
    %view([150 30]) %left side view
    %view([0 90]) %top right side view
    view([50 30]) %e3v834c view
    
    % Coordinate system: Z is the vertical axis (up down), Y is the horizontal axis (left/right), and X is depth (forwards/backwards)
    % The Z and Y axis are flipped
    % The origin is between the chest and shoulder on the Y axis, just above
    % the elbow on the Z axis, and really far forwards in the X axis (X is more
    % negative forwards)
    
    myVideo = VideoWriter('PCA_test.avi');
    myVideo.FrameRate = 30;
    open(myVideo)
    
    IndexPlot = plot3([T.IT_z(startFrame),T.IDIP_z(startFrame),T.IPIP_z(startFrame),T.IMCP_z(startFrame),T.RW_z(startFrame)], ...
        [T.IT_x(startFrame),T.IDIP_x(startFrame),T.IPIP_x(startFrame),T.IMCP_x(startFrame),T.RW_x(startFrame)], ...
        [T.IT_y(startFrame),T.IDIP_y(startFrame),T.IPIP_y(startFrame),T.IMCP_y(startFrame),T.RW_y(startFrame)], ...
        '-o', 'MarkerSize',3,'MarkerFaceColor','m');
    MiddlePlot = plot3([T.MT_z(startFrame),T.MDIP_z(startFrame),T.MPIP_z(startFrame),T.MMCP_z(startFrame),T.RW_z(startFrame)], ...
        [T.MT_x(startFrame),T.MDIP_x(startFrame),T.MPIP_x(startFrame),T.MMCP_x(startFrame),T.RW_x(startFrame)], ...
        [T.MT_y(startFrame),T.MDIP_y(startFrame),T.MPIP_y(startFrame),T.MMCP_y(startFrame),T.RW_y(startFrame)], ...
        '-o', 'MarkerSize',3,'MarkerFaceColor','b');
    RingPlot = plot3([T.RT_z(startFrame),T.RDIP_z(startFrame),T.RPIP_z(startFrame),T.RMCP_z(startFrame),T.RW_z(startFrame)], ...
        [T.RT_x(startFrame),T.RDIP_x(startFrame),T.RPIP_x(startFrame),T.RMCP_x(startFrame),T.RW_x(startFrame)], ...
        [T.RT_y(startFrame),T.RDIP_y(startFrame),T.RPIP_y(startFrame),T.RMCP_y(startFrame),T.RW_y(startFrame)], ...
        '-o', 'MarkerSize',3,'MarkerFaceColor','c');
    LittlePlot = plot3([T.LT_z(startFrame),T.LDIP_z(startFrame),T.LPIP_z(startFrame),T.LMCP_z(startFrame),T.RW_z(startFrame)], ...
        [T.LT_x(startFrame),T.LDIP_x(startFrame),T.LPIP_x(startFrame),T.LMCP_x(startFrame),T.RW_x(startFrame)], ...
        [T.LT_y(startFrame),T.LDIP_y(startFrame),T.LPIP_y(startFrame),T.LMCP_y(startFrame),T.RW_y(startFrame)], ...
        '-o', 'MarkerSize',3,'MarkerFaceColor','g');
    RightArmPlot = plot3([T.C_z(startFrame),T.RS_z(startFrame),T.RE_z(startFrame),T.RW_z(startFrame)], ...
         [T.C_x(startFrame),T.RS_x(startFrame),T.RE_x(startFrame),T.RW_x(startFrame)], ...
         [T.C_y(startFrame),T.RS_y(startFrame),T.RE_y(startFrame),T.RW_y(startFrame)], ...
         '-o', 'MarkerSize',3,'MarkerFaceColor',	"#D95319");
    
    % saveas(gcf,'reach_length.png')
    
    frame = getframe(gcf);
    writeVideo(myVideo, frame);
    
    % Animating the plot
    for i = (startFrame+1):(h-3)
        set(ThumbPlot, 'ZData', [T.TT_y(i),T.TIP_y(i),T.TMCP_y(i),T.TCMC_y(i),T.RW_y(i)], ...
            'XData', [T.TT_z(i),T.TIP_z(i),T.TMCP_z(i),T.TCMC_z(i),T.RW_z(i)], ...
            'YData', [T.TT_x(i),T.TIP_x(i),T.TMCP_x(i),T.TCMC_x(i),T.RW_x(i)]);
        set(IndexPlot,'ZData', [T.IT_y(i),T.IDIP_y(i),T.IPIP_y(i),T.IMCP_y(i),T.RW_y(i)], ...
            'XData', [T.IT_z(i),T.IDIP_z(i),T.IPIP_z(i),T.IMCP_z(i),T.RW_z(i)], ...
            'YData', [T.IT_x(i),T.IDIP_x(i),T.IPIP_x(i),T.IMCP_x(i),T.RW_x(i)]);
        set(MiddlePlot, 'ZData', [T.MT_y(i),T.MDIP_y(i),T.MPIP_y(i),T.MMCP_y(i),T.RW_y(i)], ...
            'XData', [T.MT_z(i),T.MDIP_z(i),T.MPIP_z(i),T.MMCP_z(i),T.RW_z(i)], ...
            'YData', [T.MT_x(i),T.MDIP_x(i),T.MPIP_x(i),T.MMCP_x(i),T.RW_x(i)]);
        set(RingPlot, 'ZData', [T.RT_y(i),T.RDIP_y(i),T.RPIP_y(i),T.RMCP_y(i),T.RW_y(i)], ...
            'XData', [T.RT_z(i),T.RDIP_z(i),T.RPIP_z(i),T.RMCP_z(i),T.RW_z(i)], ...
            'YData', [T.RT_x(i),T.RDIP_x(i),T.RPIP_x(i),T.RMCP_x(i),T.RW_x(i)]);
        set(LittlePlot, 'ZData', [T.LT_y(i),T.LDIP_y(i),T.LPIP_y(i),T.LMCP_y(i),T.RW_y(i)], ...
            'XData', [T.LT_z(i),T.LDIP_z(i),T.LPIP_z(i),T.LMCP_z(i),T.RW_z(i)], ...
            'YData', [T.LT_x(i),T.LDIP_x(i),T.LPIP_x(i),T.LMCP_x(i),T.RW_x(i)]);
        set(RightArmPlot, 'ZData', [T.C_y(i),T.RS_y(i),T.RE_y(i),T.RW_y(i)], ...
            'XData', [T.C_z(i),T.RS_z(i),T.RE_z(i),T.RW_z(i)], ...
            'YData', [T.C_x(i),T.RS_x(i),T.RE_x(i),T.RW_x(i)]);
    
        frame = getframe(gcf);
        writeVideo(myVideo, frame);
    end
    fprintf('done')
    close(myVideo)
end

function allExplained = PCA_test()%src,src_3d)
    % Finding PCs from individuals.
    category = "All";
    src = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\angle';
    src_3d = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\filter';
    
    allAngles = [];
    grasp_angles = [];
    p_sig_80 = [];
    p_latent = [];
    allExplained = [];
    allCoeff = [];
    maintemp = dir(fullfile(src,'*'));
    mainfolder = setdiff({maintemp([maintemp.isdir]).name},{'.','..'});
    % Looping through participants
    for i = 1:numel(mainfolder)
        
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
            categoryAngles = [];
            
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
                end_temp = min(reach+10,height(thisTable));                                 % Try 10 frames after reach end rather than 40
                end_time = [end_time;reach_time(end)+min(40,height(thisTable)-reach)];
        
                categoryAngles = [categoryAngles;thisTable(max([1,reach-39]):end_temp,:)];
                % video_data = readtable(fullfile(src_3d,mainfolder{i},actfolder{l},posfolder{j}));
                % video_data = video_data(max([1,reach-39]):end_temp,:);
                % video_maker(video_data)
            end
            % do PCA here
            
            [coeff,score,latent,tsquared,explained,mu] = pca(table2array(categoryAngles), 'Rows', 'pairwise');
            allExplained = [allExplained,explained];
        end
        % allAngles = [allAngles;categoryAngles];
        % 
        % % Keeping tracking of the start and end times of each activity after
        % % concatenation
        % table_times = table(activities,start_time,reach_time,end_time,grasp_time);
        % writetable(table_times, strcat('PCA_times\',mainfolder{i},'.csv'))
        % 
        % pAngles = table2array(categoryAngles);
        % 
        % %% PCA uncomment after variance
        % [coeff,score,latent,tsquared,explained,mu] = pca(pAngles, 'Rows', 'pairwise'); %pAngles
        % p_latent = [p_latent,sum(explained(find(latent>1)))];
        % p_sig_80 = [p_sig_80,find(cumsum(explained)>80,1)]; % Finds the amount of PCs needed to explain 80% of variance
        % curr_sig_80 = find(cumsum(explained)>80,1);
        % allExplained = [allExplained,explained];
        % allCoeff(:,:,i) = coeff;
        % save(strcat('Significance\',mainfolder{i},'_sig_80.mat'),'curr_sig_80');
        % 
        % % The score matrix contains the projected data
        % table_projection = array2table(score);
        % writetable(table_projection, strcat('Projection\',mainfolder{i},'.csv'))
        % % The mu matrix contains the estimated mean of each variable, needed
        % % for reprojection
        % table_mean = array2table(mu);
        % writetable(table_mean, strcat('PCA_mean\',mainfolder{i},'.csv'))
        % % Saving the raw data
        % writetable(categoryAngles, strcat('Unprojected\',mainfolder{i},'.csv'))
        % 
        % % The PCA coefficients
        % table_coeff = array2table(coeff,'VariableNames',string(1:length(coeff)));
        % 
        % figure
        % bar(explained, 'w');
        % ylabel('Percentage of variance explained')
        % xlabel('Principal Component')
        % titlestr = strcat("PCA of grasp type: ", mainfolder{i});
        % title(titlestr, 'Interpreter', 'none')
        % set(gca,'fontsize',14, 'TickDir', 'out')
        % ylim([0 80]);
        % box off
        % writetable(table_coeff, strcat('PCA_coeffs\',mainfolder{i},'.csv'))
        % saveas(gcf,strcat('PCA_images\',mainfolder{i},'_PCA.png'))
    
    end
end

function move_files()
    % Filter input data to find bad sessions
    src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\3d';
    src2 = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\angles_participant';
    dest = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\angle_data';
    
    parttemp = dir(fullfile(src,'*'));
    partfolder = setdiff({parttemp([parttemp.isdir]).name},{'.','..'});
    
    total_count = [];
    for participant_i = 1:numel(partfolder)
        fprintf("participant %d\n",participant_i)
        % Extra loop incase of categorised 3d data
        par_count = [];
        acttemp = dir(fullfile(src,partfolder{participant_i},'*'));
        actfolder = setdiff({acttemp([acttemp.isdir]).name},{'.','..'});
        for activity_k = 1:numel(actfolder)
            filtertemp = dir(fullfile(src2,partfolder{participant_i},actfolder{activity_k},'*.csv'));
            filterfolder = {filtertemp(~[filtertemp.isdir]).name};
            trialtemp = dir(fullfile(src,partfolder{participant_i},actfolder{activity_k},'*'));
            trialfolder = setdiff({trialtemp([trialtemp.isdir]).name},{'.','..'});
            trial_count = 0;
            adjust = 0;
            for trial_j = 1:numel(trialfolder)
                if isfile(fullfile(src,partfolder{participant_i},actfolder{activity_k},trialfolder{trial_j},'df_3d.csv')) && isfile(fullfile(src,partfolder{participant_i},actfolder{activity_k},trialfolder{trial_j},'df_err.csv'))
                    data_err = fullfile(src,partfolder{participant_i},actfolder{activity_k},trialfolder{trial_j},'df_err.csv');
                    df_err = readtable(data_err);
                    df_err([1:5,end-4:end],:) = [];
                    df_err = df_err(:,2:26);
                    df_err = table2array(df_err);
                    df_err_comp = df_err;
                    df_err_comp_size = sum(sum(~isnan(df_err)));
                    df_err(df_err>30) = NaN;
                    df_err_size = sum(sum(~isnan(df_err)));
                    if df_err_size > df_err_comp_size*0.9
                        trial_count = trial_count + 1;
                        if ~isfolder(fullfile(dest,partfolder{participant_i},actfolder{activity_k}))
                            mkdir(fullfile(dest,partfolder{participant_i},actfolder{activity_k}));
                        end
                        copyfile(fullfile(src2,partfolder{participant_i},actfolder{activity_k},filterfolder{trial_j-adjust}), fullfile(dest,partfolder{participant_i},actfolder{activity_k},filterfolder{trial_j-adjust}))
                    elseif ismember(partfolder{participant_i},['p14','p15','p16','p17'])
                        if ~isfolder(fullfile(dest,partfolder{participant_i},actfolder{activity_k}))
                            mkdir(fullfile(dest,partfolder{participant_i},actfolder{activity_k}));
                        end
                        copyfile(fullfile(src2,partfolder{participant_i},actfolder{activity_k},filterfolder{trial_j-adjust}), fullfile(dest,partfolder{participant_i},actfolder{activity_k},filterfolder{trial_j-adjust}))
                    end
                else
                    adjust = adjust + 1;
                end
            end
            par_count = [par_count;trial_count];
        end
        total_count = [total_count,par_count];
    end
end

function [rating,part_rating,task_rating] = PC_rating(explained,PC_num,limit)
    rating = zeros(size(explained,2),1);
    part_rating = zeros(18,1);
    part_tracker = 1;
    task_rating = zeros(80,1);
    task_tracker = 1;
    for i = 1:size(explained,2)
        if sum(explained(1:PC_num,i)) > limit
            rating(i) = 1;
            part_rating(part_tracker) = part_rating(part_tracker) + 1;
            task_rating(task_tracker) = task_rating(task_tracker) + 1;
        end
        part_tracker = part_tracker + 1;
        if part_tracker > 18
            part_tracker = 1;
        end
        task_tracker = task_tracker + 1;
        if task_tracker > 80
            task_tracker = 1;
        end
    end
end

close all; clear all; clc;
PCA_act();
% src = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\old\';
% dest = 'C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\filter_new\';
% PCA_filt(src,dest,30,20);
% [rating2,part_rating2,task_rating2] = PC_rating(ex,3,80);
% src = 'C:\Users\czhe0008\Documents\MATLAB\Openpose\per_par\PCA_coeffs';
% allPCs = [];
% parttemp = dir(fullfile(src,'*.csv'));
% partfolder = {parttemp(~[parttemp.isdir]).name};
% for participant_i = 1:numel(partfolder)
%     data = fullfile(src,partfolder{participant_i});
%     df = readtable(data);
%     df = table2array(df);
%     df = df(2:end,1:4);
%     allPCs = [allPCs,df];
% end
% PCA_cluster(allPCs');

% PCA_filt('C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\3d','C:\Users\czhe0008\Documents\DLCprojects\MATLAB\Clustering\filter3',30,20)