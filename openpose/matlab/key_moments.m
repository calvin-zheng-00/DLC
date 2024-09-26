close all; clear all; clc;

colnames = categorical(["TCMC","TMCP","TIP","TT","IMCP","IPIP","IDIP","IT","MMCP","MPIP","MDIP","MT",...
    "RMCP","RPIP","RDIP","RT","LMCP","LPIP","LDIP","LT","C","RS","RE","RW"]);
colnames = reordercats(colnames,string(colnames));

reach_speed = 3;

src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\3d\p11_3d_sorted';
participant_src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\3d';

participanttemp = dir(fullfile(participant_src,'*'));
participantfolder = setdiff({participanttemp([participanttemp.isdir]).name},{'.','..'});

maintemp = dir(fullfile(src,'*'));
mainfolder = setdiff({maintemp([maintemp.isdir]).name},{'.','..'});
% Looping through activities
for k = 1:numel(mainfolder)
    alignment = [];
    for p = 1:numel(participantfolder)
        subtemp = dir(fullfile(participant_src,participantfolder{p},mainfolder{k},'*'));
        subfolder = setdiff({subtemp([subtemp.isdir]).name},{'.','..'});
        % subtemp = dir(fullfile(src,mainfolder{k},'*'));
        % subfolder = setdiff({subtemp([subtemp.isdir]).name},{'.','..'});
        % starttimes = [];
        % endtimes = [];
        % tablelist = [];
        for j = 1:numel(subfolder)
            File = fullfile(participant_src,participantfolder{p},mainfolder{k},subfolder{j},'df_3d.csv');
            File = readtable(File);
            T = File;
            for i = 0:((width(T)/3)-1)
                T{:,i*3 + 1} = T{:,i*3 + 1} - File{:,76};
                T{:,i*3 + 2} = T{:,i*3 + 2} - File{:,77};
                T{:,i*3 + 3} = T{:,i*3 + 3} - File{:,78};
            end
            T(end-9:end,:) = [] ;
            T(1:10,:) = [] ;
            arr = table2array(T);
            vel = arr(1:end-1,:)-arr(2:end,:);
            [~, vel_rise] = max(vel(1:end-1, 78) > reach_speed & vel(2:end, 78) > reach_speed);
            [~, vel_fall] = max(vel(vel_rise:end-1, 78) < reach_speed & vel(vel_rise+1:end, 78) < reach_speed);
            vel_fall = vel_fall + vel_rise;
            alignment = [alignment;arr(vel_fall,:)];
        end
    end
    %Removing unused variables from alignment
    alignment(:,[1:3,64:66,79:end]) = [];

    average = mean(alignment,1,"omitnan");
    variance = var(alignment,0,1,"omitnan");
    standard_dev = sqrt(variance);
    figure;
    hold on;
    scatter3(average(3:3:end),average(1:3:end),average(2:3:end))
    xlabel('x (mm)')
    ylabel('y (mm)')
    zlabel('z (mm)')
    plot3([average(3:3:end);average(3:3:end)], [average(1:3:end);average(1:3:end)], [average(2:3:end);average(2:3:end)]+[-standard_dev(2:3:end);standard_dev(2:3:end)], '-r')
    plot3([average(3:3:end);average(3:3:end)], [average(1:3:end);average(1:3:end)]+[-standard_dev(1:3:end);standard_dev(1:3:end)], [average(2:3:end);average(2:3:end)], '-r')
    plot3([average(3:3:end);average(3:3:end)]+[-standard_dev(3:3:end);standard_dev(3:3:end)], [average(1:3:end);average(1:3:end)], [average(2:3:end);average(2:3:end)], '-r')
    set(gca, 'XDir','reverse')
    set(gca, 'ZDir','reverse')
    set(gca, 'fontsize',14)
    set(gcf,'position',[200,200,600,570])
    view([50 30])
    figure;
    bar(colnames,variance(3:3:end))
    xlabel('joint')
    ylabel('variance (mm^2)')
    title('x variance during grab')
    figure;
    bar(colnames,variance(1:3:end))
    xlabel('joint')
    ylabel('variance (mm^2)')
    title('y variance during grab')
    figure;
    bar(colnames,variance(2:3:end))
    xlabel('joint')
    ylabel('variance (mm^2)')
    title('z variance during grab')
    flag = 1;
end


