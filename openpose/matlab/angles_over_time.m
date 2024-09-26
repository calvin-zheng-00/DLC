close all; clear all; clc;

% Reading in data
File = readtable('C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\data\high_thresh\filtered\p4_2_high_angles\20240424T112749-112821_angles.csv');
joints = ["Thumb carpormetacarpal flexion","Thumb metacarpophalangeal flexion","Thumb interphalangeal flexion","Index metacarpophalangeal flexion",...
    "Index proximal interphalangeal flexion","Index distal interphalangeal flexion","Middle metacarpophalangeal flexion","Middle proximal interphalangeal flexion",...
    "Middle distal interphalangeal flexion","Ring metacarpophalangeal flexion","Ring proximal interphalangeal flexion","Ring distal interphalangeal flexion",...
    "Little metacarpophalangeal flexion","Little proximal interphalangeal flexion","Little distal interphalangeal flexion","Elbow flexion",...
    "Index abduction","Middle abduction","Ring abduction","Little abduction","Thumb abduction","Thumb rotation", "Wrist flexion","Wrist Abduction",...
    "Wrist rotation","Shoulder flexion","Shoulder abduction","Shoulder rotation"];
titles = ["Thumb","Index","Middle","Ring","Little","Elbow","Finger abduction","Thumb base","Wrist","Shoulder"];
thumbtitles = ["Thumb carpormetacarpal flexion","Thumb metacarpophalangeal flexion","Thumb interphalangeal flexion","Thumb metacarpophalangeal abduction","Thumb carpormetacarpal rotation",];
%colours = ['b','g','m','r','c','k','y'];
colours = ["#0012A7","#00A5F3","#00F3FF","#A4FCFF","#0B5ACC"];
counter = 1;
freq = 40;
phase = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
phase_count = 1;
num_trial = 1;

angles = table2array(File);
% angles = angles(1:240,:);
% 1-240, 280-520, 600-800, 880-1080, 1160-end
%trials = [1,240;305,545;600,800;880,1080;1160,1395];

%% Figuring out phases (og)
% Using E_flex to determine which phase we are in
% < 1.5 = resting, > 2 = grasping
% boundaries = [1.5,2];
% flag = 1;
% skip = false;
% for i = 1:height(angles)
%     if skip
%         skip = false;
%         continue;
%     end
%     if flag == 1
%         if File.RE_flex(i) > boundaries(1)
%             phase(phase_count) = i;
%             phase_count = phase_count + 1;
%             flag = 2;
%             skip = true;
%         end
%     elseif flag == 2
%         if File.RE_flex(i) > boundaries(2)
%             phase(phase_count) = i;
%             phase_count = phase_count + 1;
%             flag = 3;
%             skip = true;
%         end
%     elseif flag == 3
%         if File.RE_flex(i) < boundaries(2)
%             phase(phase_count) = i;
%             phase_count = phase_count + 1;
%             flag = 4;
%             skip = true;
%         end
%     elseif flag == 4
%         if File.RE_flex(i) < boundaries(1)
%             phase(phase_count) = i;
%             phase_count = phase_count + 1;
%             flag = 1;
%             skip = true;
%         end
%     end
% end

%% Test phases
%trials = [phase(1),phase(4);phase(6)-10,phase(8);phase(9),phase(12);phase(13),phase(16);phase(17),phase(20)] + [-40,40];
% phase(1) = phase(1) - 16;
% phase(2) = phase(2) - 20;
% phase(3) = phase(3) - 10;
% phase(4) = phase(4) - 10;
% phase = phase/freq;
% phase = [0.850000000000000	1.22500000000000	4.05000000000000	4.35000000000000	5.77500000000000	9.07500000000000	12.2000000000000	12.5750000000000	16.3000000000000	16.6500000000000	19.2000000000000	19.4750000000000	23.4000000000000	23.9250000000000	26.3250000000000	26.6750000000000	30.1500000000000	31.1250000000000	31.7500000000000	33.4750000000000	34.9000000000000];
% trials = [10	224;313	543; 612	819;896	1107;1166	1379];

%% Separating trials
T = readtable('C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\data\high_thresh\filtered\p4_2_high_filt\20240424T112749-112821_filtered.csv');

h = height(T); % Number of rows in table
start_buffer = 20;

% Calculating different phases
count = 0;
stage = 0;
start_trial = [];
end_trial = [];
for i=1:h
    if stage == 0
        if (T.RW_x(i) > (T.RW_x(1)+start_buffer)) && (T.RW_y(i) < (T.RW_y(1)-start_buffer)) && (T.RW_z(i) < (T.RW_z(1)-start_buffer))
            count = count + 1;
            if count == 40
                count = 0;
                stage = 1;
                start_trial = [start_trial, i-40];
            end
        else
            count = 0;
        end
    else
        if (T.RW_x(i) < (T.RW_x(1)+start_buffer)) && (T.RW_y(i) > (T.RW_y(1)-start_buffer)) && (T.RW_z(i) > (T.RW_z(1)-start_buffer))
            count = count + 1;
            if count == 40
                count = 0;
                stage = 0;
                end_trial = [end_trial, i-40];
            end
        else
            count = 0;
        end
    end
end

valid = 0;
if length(start_trial) == length(end_trial)
    if length(start_trial) == 5
        valid = 1;
    elseif length(start_trial) == 10
        valid = 2;
    end
end

%% Angles
angles = (angles*180)/pi;
%trials = [1,length(angles);1,1];
trials = [start_trial;end_trial]';
num_trial = length(start_trial);
phase = [0,0,0,0];
%figure size = [100,300,1600,400]

% Thumb flex
for i = 1:3
    figure
    hold on
    for j = 1:num_trial
        plot((1:(trials(j,2)-trials(j,1)+1))/freq,angles(trials(j,1):trials(j,2),i), 'Color', colours(j), 'LineWidth', 2)
    end
    titlestr = strcat(thumbtitles(i)," angles over time");
    title(titlestr)
    legend(["trial 1","trial 2","trial 3","trial 4","trial 5"],"AutoUpdate","off")
    ylabel('angle (degs)')
    xlabel('time (s)')
    xlim([0 (trials(1,2)-trials(1,1)+1)/freq])
    ylim([0 180])
    set(gcf,'position',[100,300,400,300])
    xline(phase(1:4),'-')
end
%Expected angles
% y1 = ones(1,(trials(1,2)-trials(1,1)+1))*0.9;
% x1 = (1:(trials(1,2)-trials(1,1)+1))/freq;
% plot(x1, y1, 'k', 'LineWidth', 5)


%     ,{'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase'})

% Thumb abd and rot
for i = 21:22
    figure
    hold on
    for j = 1:num_trial
        plot((1:(trials(j,2)-trials(j,1)+1))/freq,angles(trials(j,1):trials(j,2),i), 'Color', colours(j), 'LineWidth', 2)
    end
    titlestr = strcat(thumbtitles(i-17)," angles over time");
    title(titlestr)
    legend(["trial 1","trial 2","trial 3","trial 4","trial 5"],"AutoUpdate","off")
    ylabel('angle (degs)')
    xlabel('time (s)')
    xlim([0 (trials(1,2)-trials(1,1)+1)/freq])
    ylim([0 180])
    set(gcf,'position',[100,300,400,300])
    xline(phase(1:4),'-')
end

% All thumb
figure
hold on
for j = 1:num_trial
    for i = [1,2,3,21,22]
        if i > 4
            colour = colours(i-17);
        else
            colour = colours(i);
        end
        plot((1:(trials(j,2)-trials(j,1)+1))/freq,angles(trials(j,1):trials(j,2),i), 'Color', colour, 'LineWidth', 2)
    end
end
%Expected angles
% y1 = ones(1,(trials(1,2)-trials(1,1)+1))*0.9;
% x1 = (1:(trials(1,2)-trials(1,1)+1))/freq;
% plot(x1, y1, 'k', 'LineWidth', 5)

titlestr = strcat(titles(1)," angles over time");
title(titlestr)
legend(joints([1,2,3,21,22]),"AutoUpdate","off")
ylabel('angle (degs)')
xlabel('time (s)')
xlim([0 (trials(1,2)-trials(1,1)+1)/freq])
ylim([0 180])
set(gcf,'position',[100,300,400,300])
xline(phase(1:4),'-')

%Index
figure
hold on
for j = 1:num_trial
    for i = 4:6
        plot((1:(trials(j,2)-trials(j,1)+1))/freq,angles(trials(j,1):trials(j,2),i), 'Color', colours(i-3), 'LineWidth', 2)
    end
end
%Expected angles
% y1 = ones(1,(trials(1,2)-trials(1,1)+1))*0.9;
% x1 = (1:(trials(1,2)-trials(1,1)+1))/freq;
% plot(x1, y1, 'k', 'LineWidth', 5)

titlestr = strcat(titles(2)," angles over time");
title(titlestr)
legend(joints(4:6),"AutoUpdate","off")
ylabel('angle (degs)')
xlabel('time (s)')
xlim([0 (trials(1,2)-trials(1,1)+1)/freq])
ylim([0 180])
set(gcf,'position',[100,300,400,300])
xline(phase(1:4),'-')

% Middle
figure
hold on
for j = 1:num_trial
    for i = 7:9
        plot((1:(trials(j,2)-trials(j,1)+1))/freq,angles(trials(j,1):trials(j,2),i), 'Color', colours(i-6), 'LineWidth', 2)
    end
end
%Expected angles
% y1 = ones(1,(trials(1,2)-trials(1,1)+1))*0.9;
% x1 = (1:(trials(1,2)-trials(1,1)+1))/freq;
% plot(x1, y1, 'k', 'LineWidth', 5)

titlestr = strcat(titles(3)," angles over time");
title(titlestr)
legend(joints(7:9),"AutoUpdate","off")
ylabel('angle (degs)')
xlabel('time (s)')
xlim([0 (trials(1,2)-trials(1,1)+1)/freq])
ylim([0 180])
set(gcf,'position',[100,300,400,300])
xline(phase(1:4),'-')

% Ring
figure
hold on
for j = 1:num_trial
    for i = 10:12
        plot((1:(trials(j,2)-trials(j,1)+1))/freq,angles(trials(j,1):trials(j,2),i), 'Color', colours(i-9), 'LineWidth', 2)
    end
end
titlestr = strcat(titles(4)," angles over time");
title(titlestr)
legend(joints(10:12),"AutoUpdate","off")
ylabel('angle (degs)')
xlabel('time (s)')
xlim([0 (trials(1,2)-trials(1,1)+1)/freq])
ylim([0 180])
set(gcf,'position',[100,300,400,300])
xline(phase(1:4),'-')

% Little
figure
hold on
for j = 1:num_trial
    for i = 13:15
        plot((1:(trials(j,2)-trials(j,1)+1))/freq,angles(trials(j,1):trials(j,2),i), 'Color', colours(i-12), 'LineWidth', 2)
    end
end
titlestr = strcat(titles(5)," angles over time");
title(titlestr)
legend(joints(13:15),"AutoUpdate","off")
ylabel('angle (degs)')
xlabel('time (s)')
xlim([0 (trials(1,2)-trials(1,1)+1)/freq])
ylim([0 180])
set(gcf,'position',[100,300,400,300])
xline(phase(1:4),'-')

% Elbow
figure
hold on
for j = 1:num_trial
    for i = 16:16
        plot((1:(trials(j,2)-trials(j,1)+1))/freq,angles(trials(j,1):trials(j,2),i), 'Color', colours(i-15), 'LineWidth', 2)
    end
end
titlestr = strcat(titles(6)," angles over time");
title(titlestr)
%legend(joints(16),"AutoUpdate","off")
ylabel('angle (degs)')
xlabel('time (s)')
xlim([0 (trials(1,2)-trials(1,1)+1)/freq])
ylim([0 180])
set(gcf,'position',[100,300,400,300])
set(gca,'fontsize',14)
xline(phase(1:4),'-')

% Finger abduction
figure
hold on
for j = 1:num_trial
    for i = 17:20
        plot((1:(trials(j,2)-trials(j,1)+1))/freq,angles(trials(j,1):trials(j,2),i), 'Color', colours(i-16), 'LineWidth', 2)
    end
end
titlestr = strcat(titles(7)," angles over time");
title(titlestr)
legend(joints(17:20),"AutoUpdate","off")
ylabel('angle (degs)')
xlabel('time (s)')
xlim([0 (trials(1,2)-trials(1,1)+1)/freq])
ylim([0 180])
set(gcf,'position',[100,300,400,300])
xline(phase(1:4),'-')

% Wrist
figure
hold on
for j = 1:num_trial
    for i = 23:25
        plot((1:(trials(j,2)-trials(j,1)+1))/freq,angles(trials(j,1):trials(j,2),i), 'Color', colours(i-22), 'LineWidth', 2)
    end
end
titlestr = strcat(titles(9)," angles over time");
title(titlestr)
legend(joints(23:25),"AutoUpdate","off")
ylabel('angle (degs)')
xlabel('time (s)')
xlim([0 (trials(1,2)-trials(1,1)+1)/freq])
ylim([0 180])
set(gcf,'position',[100,300,400,300])
xline(phase(1:4),'-')

% Shoulder
figure
hold on
for j = 1:num_trial
    for i = 26:28
        plot((1:(trials(j,2)-trials(j,1)+1))/freq,angles(trials(j,1):trials(j,2),i), 'Color', colours(i-25), 'LineWidth', 2)
    end
end
titlestr = strcat(titles(10)," angles over time");
title(titlestr)
legend(joints(26:28),"AutoUpdate","off")
ylabel('angle (degs)')
xlabel('time (s)')
xlim([0 (trials(1,2)-trials(1,1)+1)/freq])
ylim([0 180])
set(gcf,'position',[100,300,400,300])
xline(phase(1:4),'-')

% for i = 1:width(angles)
%     figure
%     hold on
%     plot((1:height(angles))/freq,angles(:,i)/pi)
%     titlestr = strcat(joints(counter)," angle over time");
%     title(titlestr)
%     ylabel('angle (rads)')
%     xlabel('time (s)')
%     ylim([0 1])
%     set(gcf,'position',[100,300,1600,400])
%     counter = counter + 1;
% end