close all; clear all; clc;
colnames = categorical(["Thumb cmc flexion","Thumb mcp flexion","Thumb ip flexion","Index mcp flexion",...
    "Index pip flexion","Index dip flexion","Middle mcp flexion","Middle pip flexion",...
    "Middle dip flexion","Ring mcp flexion","Ring pip flexion","Ring dip flexion",...
    "Little mcp flexion","Little pip flexion","Little dip flexion",...
    "Index abduction","Middle abduction","Ring abduction","Little abduction","Thumb abduction","Thumb rotation",...
    "Wrist flexion","Wrist Abduction","Wrist rotation",...
    "Elbow flexion","Shoulder flexion","Shoulder abduction","Shoulder rotation"]);
colnames = reordercats(colnames,string(colnames));
grasp_cat_num = [1,3,4,5,8,9,11,12,13,14,15,18,19];
grasp_cat = dictionary(grasp_cat_num,["power_palm_5_tabd","power_pad_3_tabd","power_pad_4_tabd","power_pad_5_tabd",...
    "precision_pad_2_tabd","precision_pad_3_tabd","precision_pad_5_tabd","precision_side_3_tabd","power_palm_noindex_tadd","power_palm_5_tadd",...
    "int_side_2_tadd","finger_press","palm_press"]);
% grasp = importdata("variance\grasp_perpar_peract.txt");
% reach = importdata("variance\reach_perpar_peract.txt");
grasp = load("variance\grasp_whole.mat");
grasp = grasp.grasp_var;
% lift = load("variance\lift_perpar_peract.mat");
lift = load("variance\lift_whole.mat");
lift = lift.lift_var;
% reach = load("variance\reach_perpar_peract.mat");
reach = load("variance\reach_whole.mat");
reach = reach.reach_var;
% interaction = load("variance\interaction_perpar_peract.mat");
interaction = load("variance\interaction_whole.mat");
interaction = interaction.interaction_var;
elbow = load("variance\elbow_whole.mat");
elbow = elbow.elbow;
grasp = sqrt(grasp);
lift = sqrt(lift);
reach = sqrt(reach);
interaction = sqrt(interaction);

order = load("variance\grasp_order.mat");
order = order.grasp_order;

%Reordering to move elbow column
grasp = [grasp(:,1:15,:),grasp(:,17:25,:),grasp(:,16,:),grasp(:,26:28,:)]; % participants, joints, activities
lift = [lift(:,1:15,:),lift(:,17:25,:),lift(:,16,:),lift(:,26:28,:)];
reach = [reach(:,1:15,:,:),reach(:,17:25,:,:),reach(:,16,:,:),reach(:,26:28,:,:)]; % Time points, joints, participants, activities
interaction = [interaction(:,1:15,:,:),interaction(:,17:25,:,:),interaction(:,16,:,:),interaction(:,26:28,:,:)];
% grasp = [grasp(:,1:15),grasp(:,17:25),grasp(:,16),grasp(:,26:28)];
% reach = [reach(:,1:15),reach(:,17:25),reach(:,16),reach(:,26:28)];

%% Activities (numbering goes like: 1,10,11,2,3,...)
% for i = 1:size(reach,4)
%     grasp_mean = mean(grasp(:,:,i),1,"omitnan");
%     reach_mean = mean(reach(:,:,:,i),3,"omitnan");
%     int_mean = mean(interaction(:,:,:,i),3,"omitnan");
%     all = [reach_mean;grasp_mean;int_mean];
% 
%     figure;
%     hold on;
%     average = mean(all,2,"omitnan");
%     plot((1:height(average))/40,all,'k')
%     plot((1:height(average))/40,average,'r')
%     xline((length(average)-40)/40, '-', {'Grasp'})
%     xlabel('time (s)')
%     titlestr = strcat("Average standard deviation across activity: ", string(i));
%     title(titlestr)
%     ylabel('standard deviation')
%     %ylim([0 60])
%     box off
%     set(gcf,'position',[200,200,800,570])
%     % saveas(gcf,strcat('C:\Users\czhe0008\Documents\DLCprojects\openpose\figures\matlab_gen\best_activity\variance_',string(i),'.png'))
% end

%% STD difference boxplot for each joint
% grasp_mean = squeeze(mean(grasp,1,"omitnan"))';
% lift_mean = squeeze(mean(lift,1,"omitnan"))';
% std_diff = lift_mean-grasp_mean;
% total_rank = [];
% for i = 1:28
%     [p, h] = ranksum(lift_mean(:,i),grasp_mean(:,i));
%     total_rank = cat(1,total_rank, p);
% end
% figure
% hold on
% temp = total_rank<0.05;
% sig = 1:28;
% sig = sig(temp);
% boxchart(std_diff)
% plot(sig, ones(length(sig),1)*(max(max(std_diff)+3)), "*k")
% yline(0)
% xline(3.5,"--","thumb flex",'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
% xline(6.5,"--","index flex",'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
% xline(9.5,"--","middle flex",'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
% xline(12.5,"--","ring flex",'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
% xline(15.5,"--","little flex",'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
% xline(21.5,"--","finger abduction",'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
% xline(24.5,"--","wrist",'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
% xline(28.5,"--","arm",'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
% xticklabels(colnames)
% legend('difference', 'significant','Location','northeastoutside')
% xlabel('Joint')
% ylabel('standard deviation (degrees)')
% set(gcf, 'Position', get(0, 'Screensize'));

%% Normalized difference boxplot for each joint
% figure;
% hold on;
% grasp_mean = squeeze(mean(grasp,1,"omitnan"))';
% lift_mean = squeeze(mean(lift,1,"omitnan"))';
% grasp_norm = normalize(grasp_mean);
% lift_norm = normalize(lift_mean);
% total_rank = [];
% for i = 1:28
%     [p, h] = ranksum(lift_mean(:,i),grasp_mean(:,i));
%     total_rank = cat(1,total_rank, p);
% end
% temp = total_rank<0.05;
% sig = 1:28;
% sig = sig(temp);
% boxchart(grasp_norm)
% boxchart(lift_norm)
% plot(sig, ones(length(sig),1)*(max(max([grasp_norm,lift_norm])+1)), "*k")
% yline(0)
% xline(3.5,"--","thumb flex",'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
% xline(6.5,"--","index flex",'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
% xline(9.5,"--","middle flex",'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
% xline(12.5,"--","ring flex",'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
% xline(15.5,"--","little flex",'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
% xline(21.5,"--","finger abduction",'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
% xline(24.5,"--","wrist",'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
% xline(28.5,"--","arm",'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
% xticklabels(colnames)
% legend('grasp', 'reach','Location','northeastoutside')
% xlabel('Joint')
% ylabel('standard deviation (degrees)')
% set(gcf, 'Position', get(0, 'Screensize'));

%% Normalized 2.0
% figure;
% hold on;
% grasp_mean = squeeze(mean(grasp,1,"omitnan"))';
% lift_mean = squeeze(mean(lift,1,"omitnan"))';
% grasp_mean_2 = squeeze(mean(grasp_mean,1,"omitnan"));
% grasp_std = squeeze(std(grasp_mean,1,"omitnan"));
% % grasp_norm = normalize(grasp_mean);
% lift_norm = (lift_mean-grasp_mean_2)./grasp_std;
% total_rank = [];
% for i = 1:28
%     [p, h] = ranksum(lift_mean(:,i),grasp_mean(:,i));
%     total_rank = cat(1,total_rank, p);
% end
% temp = total_rank<0.05;
% sig = 1:28;
% sig = sig(temp);
% boxchart(lift_norm)
% plot(sig, ones(length(sig),1)*(max(max(lift_norm)+1)), "*k")
% yline(0)
% xline(3.5,"--","thumb flex",'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
% xline(6.5,"--","index flex",'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
% xline(9.5,"--","middle flex",'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
% xline(12.5,"--","ring flex",'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
% xline(15.5,"--","little flex",'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
% xline(21.5,"--","finger abduction",'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
% xline(24.5,"--","wrist",'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
% xline(28.5,"--","arm",'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
% xticklabels(colnames)
% legend('reach','significant','Location','northeastoutside')
% xlabel('Joint')
% ylabel('standard deviation (degrees)')
% set(gcf, 'Position', get(0, 'Screensize'));

%% Comparison of general variance at different reach times
grasp_mean = squeeze(mean(mean(grasp,1,"omitnan"),2,"omitnan"))';
reach_mean = squeeze(mean(mean(reach,3,"omitnan"),2,"omitnan"));
int_mean = squeeze(mean(mean(interaction,3,"omitnan"),2,"omitnan"));
whole = [reach_mean;grasp_mean;int_mean];
% std_diff = reach_mean-grasp_mean;
% total_rank = [];
% for i = 1:size(reach_mean,1)
%     [p, h] = ranksum(squeeze(reach_mean(i,:)),grasp_mean);
%     total_rank = cat(1,total_rank, p);
% end
% figure
% hold on
% temp = total_rank<0.05;
% sig = 1:size(reach_mean,1);
% sig = sig(temp);
% plot((1:size(std_diff,1))/40,std_diff,"k")
% plot((1:size(std_diff,1))/40,mean(std_diff,2),"r")
% plot(sig/40, ones(length(sig),1)*(max(max(std_diff)+3)), "*b")
% legend('std difference', 'significant')
% xlabel('Time (s)')
% ylabel('standard deviation difference (degrees)')
% set(gcf, 'Position', get(0, 'Screensize'));


% [h,p] = kstest(grasp_mean);

test2 = mean(elbow,2,"omitnan");
test2 = (test2 - min(test2))/2;


for i = grasp_cat_num
    filter = order == i;
    temp = whole .* filter;
    temp( :, all(~temp,1) ) = [];
    figure
    hold on
    plot((1:size(temp,1))/40-1,test2,"b")
    plot((1:size(temp,1))/40-1,temp,"k")
    set(0, 'DefaultTextInterpreter', 'none')
    titlestr = strcat("Average variance of the hand for activity group: ", grasp_cat(i));
    title(titlestr)
    xline(0,"k","grasp")
    legend('elbow','std')
    xlabel('Time (s)')
    ylabel('standard deviation difference (degrees)')
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf,strcat('C:\Users\czhe0008\Documents\DLCprojects\openpose\figures\matlab_gen\boxchart\raw_std_overtime\raw_std_',grasp_cat(i),'.png'))
end

%% STD of all joints
% grasp_mean = squeeze(mean(grasp,1,"omitnan"))';
% lift_mean = squeeze(mean(lift,1,"omitnan"))';
% std_diff = squeeze(mean(lift_mean,2,"omitnan"))-squeeze(mean(grasp_mean,2,"omitnan"));
% figure
% hold on
% boxchart(std_diff)
% xlabel("Average hand standard deviation difference")
% ylabel('standard deviation (degrees)')
% set(gca, 'xticklabel', {[]});
% set(gcf, 'Position', get(0, 'Screensize'));

%% STD of both reach and grasp
% total_rank = [];
% for i = 1:28
%     [p, h] = ranksum(reach(:,i),grasp(:,i));
%     m = mean(reach(:,i))-mean(grasp(:,i));
%     total_rank = cat(1,total_rank, [p, m]);
% end
% figure
% hold on
% temp = total_rank(:,1)<0.05;
% sig = 1:28;
% sig = sig(temp);
% boxchart(grasp)
% boxchart(reach)
% plot(sig, ones(length(sig),1)*55, "*k")
% xline(3.5,"--","thumb flex",'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
% xline(6.5,"--","index flex",'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
% xline(9.5,"--","middle flex",'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
% xline(12.5,"--","ring flex",'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
% xline(15.5,"--","little flex",'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
% xline(21.5,"--","finger abduction",'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
% xline(24.5,"--","wrist",'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
% xline(28.5,"--","arm",'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
% xticklabels(colnames)
% legend('grasp', 'reach', 'significant','Location','northeastoutside')
% xlabel('Joint')
% ylabel('standard deviation (degrees)')
% set(gcf, 'Position', get(0, 'Screensize'));
