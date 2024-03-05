close all; clear all; clc;

% Reading in data
File = readtable('activity4_angles.csv');
joints = ["Thumb carpormetacarpal flexion","Thumb metacarpophalangeal flexion","Thumb interphalangeal flexion","Index metacarpophalangeal flexion",...
    "Index proximal interphalangeal flexion","Index distal interphalangeal flexion","Middle metacarpophalangeal flexion","Middle proximal interphalangeal flexion",...
    "Middle distal interphalangeal flexion","Ring metacarpophalangeal flexion","Ring proximal interphalangeal flexion","Ring distal interphalangeal flexion",...
    "Little metacarpophalangeal flexion","Little proximal interphalangeal flexion","Little distal interphalangeal flexion","Elbow flexion",...
    "Index and middle abduction","Middle and ring abduction","Ring and little abduction","Thumb and index abduction","Wrist flexion","Wrist Abduction",...
    "Wrist rotation","Shoulder flexion","Shoulder rotation"];
titles = ["Thumb","Index","Middle","Ring","Little","Elbow","Inter-digit abduction","Wrist","Shoulder"];
colours = ['b','g','m','r','c','k','y'];
counter = 1;
freq = 40;
phase = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
phase_count = 1;

angles = table2array(File);
% angles = angles(1:240,:);
% 1-240, 280-520, 600-800, 880-1080, 1160-end
trials = [1,240;280,520;600,800;880,1080;1160,1395];

% Using E_flex to determine which phase we are in
% < 1.5 = resting, > 2 = grasping
boundaries = [1.5,2];
flag = 1;
skip = false;
for i = 1:height(angles)
    if skip
        skip = false;
        continue;
    end
    if flag == 1
        if File.RE_flex(i) > boundaries(1)
            phase(phase_count) = i;
            phase_count = phase_count + 1;
            flag = 2;
            skip = true;
        end
    elseif flag == 2
        if File.RE_flex(i) > boundaries(2)
            phase(phase_count) = i;
            phase_count = phase_count + 1;
            flag = 3;
            skip = true;
        end
    elseif flag == 3
        if File.RE_flex(i) < boundaries(2)
            phase(phase_count) = i;
            phase_count = phase_count + 1;
            flag = 4;
            skip = true;
        end
    elseif flag == 4
        if File.RE_flex(i) < boundaries(1)
            phase(phase_count) = i;
            phase_count = phase_count + 1;
            flag = 1;
            skip = true;
        end
    end
end
phase = phase/freq;
angles = (angles*180)/pi;
%figure size = [100,300,1600,400]

% Thumb
for i = 1:3
    figure
    hold on
    for j = 1:5
        plot((1:(trials(j,2)-trials(j,1)+1))/freq,angles(trials(j,1):trials(j,2),i), colours(i))
    end
    titlestr = strcat(titles(1)," angles over time");
    title(titlestr)
    legend(joints(i),"AutoUpdate","off")
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

%Index
figure
hold on
for j = 1:5
    for i = 4:6
        plot((1:(trials(j,2)-trials(j,1)+1))/freq,angles(trials(j,1):trials(j,2),i), colours(i-3))
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
%     ,{'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase'})

figure
hold on
for j = 1:5
    for i = 7:9
        plot((1:(trials(j,2)-trials(j,1)+1))/freq,angles(trials(j,1):trials(j,2),i), colours(i-6))
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
%     ,{'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase'})

figure
hold on
for j = 1:5
    for i = 10:12
        plot((1:(trials(j,2)-trials(j,1)+1))/freq,angles(trials(j,1):trials(j,2),i), colours(i-9))
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
%     ,{'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase'})

figure
hold on
for j = 1:5
    for i = 13:15
        plot((1:(trials(j,2)-trials(j,1)+1))/freq,angles(trials(j,1):trials(j,2),i), colours(i-12))
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
%     ,{'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase'})

figure
hold on
for j = 1:5
    for i = 16:16
        plot((1:(trials(j,2)-trials(j,1)+1))/freq,angles(trials(j,1):trials(j,2),i), colours(i-15))
    end
end
titlestr = strcat(titles(6)," angles over time");
title(titlestr)
legend(joints(16),"AutoUpdate","off")
ylabel('angle (degs)')
xlabel('time (s)')
xlim([0 (trials(1,2)-trials(1,1)+1)/freq])
ylim([0 180])
set(gcf,'position',[100,300,400,300])
xline(phase(1:4),'-')
%     ,{'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase'})

figure
hold on
for j = 1:5
    for i = 17:20
        plot((1:(trials(j,2)-trials(j,1)+1))/freq,angles(trials(j,1):trials(j,2),i), colours(i-16))
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
%     ,{'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase'})

figure
hold on
for j = 1:5
    for i = 21:23
        plot((1:(trials(j,2)-trials(j,1)+1))/freq,angles(trials(j,1):trials(j,2),i), colours(i-20))
    end
end
titlestr = strcat(titles(8)," angles over time");
title(titlestr)
legend(joints(21:23),"AutoUpdate","off")
ylabel('angle (degs)')
xlabel('time (s)')
xlim([0 (trials(1,2)-trials(1,1)+1)/freq])
ylim([0 180])
set(gcf,'position',[100,300,400,300])
xline(phase(1:4),'-')
%     ,{'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase'})

figure
hold on
for j = 1:5
    for i = 23:24
        plot((1:(trials(j,2)-trials(j,1)+1))/freq,angles(trials(j,1):trials(j,2),i), colours(i-22))
    end
end
titlestr = strcat(titles(9)," angles over time");
title(titlestr)
legend(joints(23:24),"AutoUpdate","off")
ylabel('angle (degs)')
xlabel('time (s)')
xlim([0 (trials(1,2)-trials(1,1)+1)/freq])
ylim([0 180])
set(gcf,'position',[100,300,400,300])
xline(phase(1:4),'-')
%     ,{'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase',...
%     'reaching phase','grasping phase','releasing phase','resting phase'})

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