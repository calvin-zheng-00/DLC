close all; clear all; clc;

colnames = categorical(["TCMC","Thumb knuckle","Thumb 1","Thumb tip","Index knuckle","Index 1","Index 2","Index tip","Middle knuckle","Middle 1","Middle 2","Middle tip",...
    "Ring knuckle","Ring 1","Ring 2","Ring tip","Little knuckle","Little 1","Little 2","Little tip","Chest","Shoulder","Elbow","Wrist"]);
%colnames = categorical(["Thumb tip","Index tip","Middle tip","Ring tip","Little tip"]);
colnames = reordercats(colnames,string(colnames));
joints = ["Thumb carpormetacarpal flexion","Thumb metacarpophalangeal flexion","Thumb interphalangeal flexion","Index metacarpophalangeal flexion",...
    "Index proximal interphalangeal flexion","Index distal interphalangeal flexion","Middle metacarpophalangeal flexion","Middle proximal interphalangeal flexion",...
    "Middle distal interphalangeal flexion","Ring metacarpophalangeal flexion","Ring proximal interphalangeal flexion","Ring distal interphalangeal flexion",...
    "Little metacarpophalangeal flexion","Little proximal interphalangeal flexion","Little distal interphalangeal flexion","Elbow flexion",...
    "Index abduction","Middle abduction","Ring abduction","Little abduction","Thumb abduction","Thumb rotation", "Wrist flexion","Wrist Abduction",...
    "Wrist rotation","Shoulder flexion","Shoulder abduction","Shoulder rotation"];
grasp_cat_num = [1,3,4,5,8,9,11,12,13,14,15,18,19];
grasp_cat = dictionary(grasp_cat_num,["power_palm_5_tabd","power_pad_3_tabd","power_pad_4_tabd","power_pad_5_tabd",...
    "precision_pad_2_tabd","precision_pad_3_tabd","precision_pad_5_tabd","precision_side_3_tabd","power_palm_noindex_tadd","power_palm_5_tadd",...
    "int_side_2_tadd","finger_press","palm_press"]);
angcol = categorical(joints);
angcol = reordercats(angcol,string(angcol));

reach_speed = 10;
magic = 6;
post_grasp_time = 40;
focus_joint = 26; %26 for wrist, 23 for chest
focus_joint_3 = 76; % 76 for wrist, 67 for chest
contact_point = 76;
contact_angle = 16; % Elbow 16, MDIP 9, IDIP 6
graph_height = 2; % 200 for variance, 2 for SEM

grasptypes = readtable('C:\Users\czhe0008\Documents\MATLAB\Openpose\grasp_type.csv');
grasps = dictionary(grasptypes.Var1,grasptypes.Var2);

src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\filter_participant';
angles_src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\angles_participant';
err_src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\3d';

maintemp = dir(fullfile(src,'*'));
mainfolder = setdiff({maintemp([maintemp.isdir]).name},{'.','..'});
acttemp = dir(fullfile(src,'p23','*'));
actfolder = setdiff({acttemp([acttemp.isdir]).name},{'.','..'});
grasp_var = [];
reach_var = [];
lift_var = [];
interaction_var = [];
elbow = [];
grasp_order = [];
% multi_joint = [];

% err_table = array2table(zeros(80,14));
% err_table.Properties.VariableNames = {'p10','p11','p12','p13','p14',...
%             'p15','p16','p17','p18','p19','p20','p21','p22','p23'};
% All joints, one person, one activity. Variance between trials.

%% Looping through activities
for l = 1:numel(actfolder)
    % key_var = [];
    % key_var_angle = [];
    grasp_order = cat(2,grasp_order,grasps(cellstr(actfolder{l})));
    % Looping through participants
    grasp_temp = [];
    reach_temp = [];
    lift_temp = [];
    interaction_temp = [];
    for k = 1:numel(mainfolder)
        % if ismember(k,[4,5,6,7,8,9,10,11,12])
        %     continue
        % end
        reach_start = [];
        lift_start = [];
        grasp_start = [];
        post_grasp = [];
        files = dir(fullfile(src,mainfolder{k},actfolder{l},'*.csv'));
        anglefiles = dir(fullfile(angles_src,mainfolder{k},actfolder{l},'*.csv'));
        if numel(files) < 1
            grasp_start = cat(1,grasp_start,NaN(1,28));
            lift_start = cat(1,lift_start,NaN(1,28));
            reach_start = cat(3,reach_start,NaN(39,28));
            post_grasp = cat(3,post_grasp,NaN(post_grasp_time,28));
        end
        for j = 1:numel(files)
            baseFileName = files(j).name;
            %% Reading files
            T = fullfile(src,mainfolder{k},actfolder{l},baseFileName);
            T = readtable(T);
            %Looking through angle files
            angleFileName = anglefiles(j).name;
            A = fullfile(angles_src,mainfolder{k},actfolder{l},angleFileName);
            A = readtable(A);
            A(end-9:end,:) = [] ;
            A(1:10,:) = [] ;
            angle_arr = table2array(A);
            angle_arr = angle_arr.*(180/pi);
    
            % Clipping off first and last few frames.
            T(end-9:end,:) = [] ;
            T(1:10,:) = [] ;
            arr = table2array(T);
            % Finding first reduction in velocity after the first rise.
            % This should ideally be right after the reach phase.
            vel = arr(1:end-1,:)-arr(2:end,:);
            vel = sqrt((vel(:,1:3:end)).^2 + (vel(:,2:3:end)).^2 + (vel(:,3:3:end)).^2);
            [~, vel_rise] = max(vel(1:end-1, focus_joint) > reach_speed & vel(2:end, focus_joint) > reach_speed);
            [~, vel_fall] = max(vel(vel_rise:end-1, focus_joint) < reach_speed & vel(vel_rise+1:end, focus_joint) < reach_speed);
            reach = vel_fall + vel_rise + magic;
            if reach < 40
                grasp_start = cat(1,grasp_start,NaN(1,28));
                lift_start = cat(1,lift_start,NaN(1,28));
                reach_start = cat(3,reach_start,NaN(39,28));
                post_grasp = cat(3,post_grasp,NaN(post_grasp_time,28));
                continue
            end
            %% Error checking
            err_folder = baseFileName(1:22);
            err = fullfile(err_src,mainfolder{k},actfolder{l},err_folder,'df_err.csv');
            err = readtable(err);
            err(end-9:end,:) = [] ;
            err(1:10,:) = [] ;
            err(:,[1,22,27:46]) = [];
            cutoff = (sum(sum(table2array(err>40)))+sum(sum(isnan(table2array(err)))))/(width(err)*height(err));
            if cutoff > 0.15
                grasp_start = cat(1,grasp_start,NaN(1,28));
                lift_start = cat(1,lift_start,NaN(1,28));
                reach_start = cat(3,reach_start,NaN(39,28));
                post_grasp = cat(3,post_grasp,NaN(post_grasp_time,28));
                continue
            end
            % over_40 = sum(sum(table2array(err(:,{'C','RE','RS','RW'})>40)));
            % over_30 = sum(sum(table2array(err(:,{'TCMC','TMCP','TIP','TT','IMCP','IPIP','IDIP','IT',...
            %     'MMCP','MPIP','MDIP','MT','RMCP','RPIP','RDIP','RT','LMCP','LPIP','LDIP','LT'})>30)));
            % cutoff = (over_40+over_30+sum(sum(isnan(table2array(err)))))/(width(err)*height(err));
            % if cutoff < 0.15
            %     err_table{l,k} = err_table{l,k} + 1;
            % else
            %     continue
            % end
            
            % key_var = cat(2, key_var, vel_fall);
    
            %% Variance of a multiple joints from start to 0.5 seconds after grasp
            % focus_angle = angle_arr;
            % if isempty(key_var_angle)
            %     key_var_angle = focus_angle(1:reach+20,:);
            % else
            %     if height(key_var_angle) <= height(focus_angle(1:reach+20,:))
            %         key_var_angle = cat(3,key_var_angle,focus_angle((reach+21-height(key_var_angle)):reach+20,:));
            %     else
            %         key_var_angle = cat(3,key_var_angle(height(key_var_angle)-reach-19:end,:,:),focus_angle(1:reach+20,:));
            %     end
            % end

            %% Alternate method of getting a slice of time
            grasp_start = cat(1,grasp_start,angle_arr(reach,:));
            lift_start = cat(1,lift_start,angle_arr(reach - vel_fall,:));
            reach_start = cat(3,reach_start,angle_arr(reach-39:reach-1,:));
            post_grasp = cat(3,post_grasp,angle_arr(reach+1:reach+post_grasp_time,:));
            %reach_start = cat(1,reach_start,angle_arr(reach+4,:)); %vel_rise
        end
        %% Find variance for each joint per participant (per activity)
        grasp_temp = cat(1,grasp_temp,var(grasp_start,0,1,"omitnan"));
        lift_temp = cat(1,lift_temp,var(lift_start,0,1,"omitnan"));
        reach_temp = cat(3,reach_temp,var(reach_start,0,3,"omitnan"));
        interaction_temp = cat(3,interaction_temp,var(post_grasp,0,3,"omitnan"));
        elbow_temp = cat(1,squeeze(reach_start(:,16,:)),grasp_start(:,16)',squeeze(post_grasp(:,16,:)));
        elbow = cat(2,elbow,elbow_temp);
    end
    %% Per activity
    grasp_var = cat(3,grasp_var,grasp_temp);
    lift_var = cat(3,lift_var,lift_temp);
    reach_var = cat(4,reach_var,reach_temp);
    interaction_var = cat(4,interaction_var,interaction_temp);
    %% Find variance of individual joints in the activity.
    % if isempty(key_var_angle)
    %     continue
    % end
    % variance = var(key_var_angle,0,3,"omitnan");
    % if isempty(multi_joint)
    %     multi_joint = variance;
    % else
    %     if height(multi_joint) <= height(variance)
    %         multi_joint = cat(3,multi_joint,variance(end-height(multi_joint)+1:end,:));
    %     else
    %         multi_joint = cat(3,multi_joint(end-height(variance)+1:end,:,:),variance);
    %     end
    % end

    %% Find variance between joints in the activity alternative method.
    % grasp_var = cat(1,grasp_var,var(grasp_start,0,1,"omitnan"));
    % reach_var = cat(1,reach_var,var(reach_start,0,1,"omitnan"));
end
% output = mean(multi_joint,2,"omitnan");
% figure;
% hold on;
% plot((1:height(output))/40,multi_joint, 'k')
% plot((1:height(output))/40,output, 'r')
% xline((length(output)-20)/40, '-', {'Grasp'})
% xline((length(output)-36)/40, '-', {'Reach Start'})
% xlabel('time (s)')
% titlestr = strcat("Average variance of each participant for each activity with ranksum p value: ", num2str(p,'%.E'));
% title(titlestr)
% ylabel('variance')
% %ylim([0 0.5])
% box off
% set(gcf,'position',[200,200,800,570])
%saveas(gcf,strcat('C:\Users\czhe0008\Documents\DLCprojects\openpose\figures\matlab_gen\variance_by_cat\variance_',joints(i),'.png'))

%writetable(err_table,'C:\Users\czhe0008\Documents\DLCprojects\openpose\figures\matlab_gen\errors.csv')

% writematrix(reach_var,"variance\reach_perpar_peract.txt")
% writematrix(interaction_var,"variance\interaction_perpar_peract.txt")
% writematrix(grasp_var,"variance\grasp_perpar_peract.txt")
save("variance\grasp_order.mat","grasp_order")
save("variance\reach_whole.mat","reach_var")
save("variance\interaction_whole.mat","interaction_var")
save("variance\grasp_whole.mat","grasp_var")
save("variance\lift_whole.mat","lift_var")
save("variance\elbow_whole.mat","elbow")

%% Ranksum

%% Printing variance overtime
% total_rank = [];
% for i = 1:28
%     [p, h] = ranksum(reach_var(:,i),grasp_var(:,i));
%     m = mean(reach_var(:,i))-mean(grasp_var(:,i));
%     total_rank = cat(1,total_rank, [p, m]);
    % figure;
    % hold on;
    % average = mean(multi_joint(:,i,:),3,"omitnan");
    % plot((1:height(average))/40,squeeze(multi_joint(:,i,:)), 'k')
    % plot((1:height(average))/40,average, 'r')
    % xline((length(average)-20)/40, '-', {'Grasp'})
    % %xline((length(average)-40)/40, '-', {'Reach Start'})
    % xlabel('time (s)')
    % set(0, 'DefaultTextInterpreter', 'none')
    % titlestr = strcat("Average ",joints(i)," variance across all activities with ranksum p value: ", num2str(p,'%.E'));
    % title(titlestr)
    % ylabel('variance')
    % %ylim([0 0.5])
    % box off
    % set(gcf,'position',[200,200,800,570])
    % saveas(gcf,strcat('C:\Users\czhe0008\Documents\DLCprojects\openpose\figures\matlab_gen\variance_by_joint_40\variance_',joints(i),'.png'))
% end
% for i = 1:length(grasp_cat_num)
%     temp = find(grasp_order==grasp_cat_num(i));
%     figure;
%     hold on;
% 
%     %% sort by joint type
%     for j = 1:28
%         average = mean(multi_joint(:,j,temp),3,"omitnan");
%         plot((1:height(average))/40,squeeze(multi_joint(:,j,temp)), 'k')
%         plot((1:height(average))/40,average, 'r')
%         xline((length(average)-20)/40, '-', {'Grasp'})
%         xline((length(average)-40)/40, '-', {'Reach Start'})
%         xlabel('time (s)')
%         set(0, 'DefaultTextInterpreter', 'none')
%         titlestr = strcat("Average ",joints(j)," variance across activity: ", grasp_cat(grasp_cat_num(i)));
%         title(titlestr)
%         ylabel('variance')
%         ylim([0 0.5])
%         box off
%         set(gcf,'position',[200,200,800,570])
%         saveas(gcf,strcat('C:\Users\czhe0008\Documents\DLCprojects\openpose\figures\matlab_gen\variance_by_cat_and_joint\variance_',joints(i),'.png'))
%         %flag = 1;
%     end
% end

%% Ranksum old
% results = [];
% i = 1:28
% for i = 1:28 %1:length(grasp_cat_num)
% 
%     average = mean(multi_joint(:,i,:),3,"omitnan");
%     time1 = squeeze(multi_joint((length(average)-36),i,:)); % -40
%     time2 = squeeze(multi_joint((length(average)-24),i,:)); % -20
%     [p, h] = ranksum(time1,time2);
%     results = cat(2,results,[p; h]);

    % temp = find(grasp_order,grasp_cat_num(i));
    % % figure;
    % % hold on;
    % average = mean(multi_joint(:,:,temp),3,"omitnan");
    % %average2 = mean(average,2,"omitnan");
    % act1 = squeeze(average((length(average)-40),:));
    % act2 = squeeze(average((length(average)-20),:));
    % [p, h] = ranksum(act1,act2);
    % results = cat(2,results,[p; h]);


    % figure;
    % hold on;
    % plot((1:height(average))/40,average, 'k')
    % xline((length(average)-24)/40, '-', {'Grasp'})
    % xline((length(average)-40)/40, '-', {'Reach Start'})
    % xlabel('time (s)')
    % set(0, 'DefaultTextInterpreter', 'none')
    % titlestr = strcat("Average ",joints(i)," variance with p value: ", string(p));
    % title(titlestr)
    % ylabel('variance')
    % %ylim([0 0.45])
    % box off
    % set(gcf,'position',[200,200,800,570])
    % saveas(gcf,strcat('C:\Users\czhe0008\Documents\DLCprojects\openpose\figures\matlab_gen\ranksum_test\ranksum_',joints(i),'.png'))
    %flag = 1;
% end