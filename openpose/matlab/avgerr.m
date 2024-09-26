close all; clear all; clc;
main_dir = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\p10_12\p10\3d';
headers = {'TCMC','TMCP','TIP','TT','IMCP','IPIP','IDIP','IT','MMCP','MPIP','MDIP','MT',...
    'RMCP','RPIP','RDIP','RT','LMCP','LPIP','LDIP','LT','C','RS','RE','RW'};
colnames = categorical(["TCMC","TMCP","TIP","TT","IMCP","IPIP","IDIP","IT","MMCP","MPIP","MDIP","MT",...
    "RMCP","RPIP","RDIP","RT","LMCP","LPIP","LDIP","LT","C","RS","RE","RW"]);
colnames = reordercats(colnames,string(colnames));
high_err_thresh = 25; % In mm

%% Individual cameras

% tot_err_0e = [];
% tot_err_1f = [];
% tot_err_37 = [];
% tot_err_3c = [];
% tot_err_4c = [];
% tot_err_80 = [];
% tot_err_89 = [];
% 
% tot_drop_0e = [];
% tot_drop_1f = [];
% tot_drop_37 = [];
% tot_drop_3c = [];
% tot_drop_4c = [];
% tot_drop_80 = [];
% tot_drop_89 = [];
% cd(main_dir)
% trials = dir(main_dir);
% for k = 3:length(trials)
%     currT = trials(k).name;
%     files = dir(currT);
%     cd(currT)
%     err_0e = readtable('e3v830e.csv');
%     err_0e = err_0e(:,headers);
%     drop_0e = sum(isnan(table2array(err_0e)));
%     high_err = sum(table2array(err_0e) > high_err_thresh);
%     drop_0e = (drop_0e+high_err)/height(table2array(err_0e));
%     err_0e = mean(table2array(err_0e), "omitnan");
%     
%     err_1f = readtable('e3v831f.csv');
%     err_1f = err_1f(:,headers);
%     drop_1f = sum(isnan(table2array(err_1f)));
%     high_err = sum(table2array(err_1f) > high_err_thresh);
%     drop_1f = (drop_1f+high_err)/height(table2array(err_1f));
%     err_1f = mean(table2array(err_1f), "omitnan");
%     
%     err_37 = readtable('e3v8337.csv');
%     err_37 = err_37(:,headers);
%     drop_37 = sum(isnan(table2array(err_37)));
%     high_err = sum(table2array(err_37) > high_err_thresh);
%     drop_37 = (drop_37+high_err)/height(table2array(err_37));
%     err_37 = mean(table2array(err_37), "omitnan");
%     
%     err_3c = readtable('e3v833c.csv');
%     err_3c = err_3c(:,headers);
%     drop_3c = sum(isnan(table2array(err_3c)));
%     high_err = sum(table2array(err_3c) > high_err_thresh);
%     drop_3c = (drop_3c+high_err)/height(table2array(err_3c));
%     err_3c = mean(table2array(err_3c), "omitnan");
%     
%     err_4c = readtable('e3v834c.csv');
%     err_4c = err_4c(:,headers);
%     drop_4c = sum(isnan(table2array(err_4c)));
%     high_err = sum(table2array(err_4c) > high_err_thresh);
%     drop_4c = (drop_4c+high_err)/height(table2array(err_4c));
%     err_4c = mean(table2array(err_4c), "omitnan");
%     
%     err_80 = readtable('e3v8380.csv');
%     err_80 = err_80(:,headers);
%     drop_80 = sum(isnan(table2array(err_80)));
%     high_err = sum(table2array(err_80) > high_err_thresh);
%     drop_80 = (drop_80+high_err)/height(table2array(err_80));
%     err_80 = mean(table2array(err_80), "omitnan");
%     
%     err_89 = readtable('e3v8389.csv');
%     err_89 = err_89(:,headers);
%     drop_89 = sum(isnan(table2array(err_89)));
%     high_err = sum(table2array(err_89) > high_err_thresh);
%     drop_89 = (drop_89+high_err)/height(table2array(err_89));
%     err_89 = mean(table2array(err_89), "omitnan");
%     
%     if size(err_0e, 2) == 24
%         tot_err_0e = [tot_err_0e;err_0e];
%         tot_err_1f = [tot_err_1f;err_1f];
%         tot_err_37 = [tot_err_37;err_37];
%         tot_err_3c = [tot_err_3c;err_3c];
%         tot_err_4c = [tot_err_4c;err_4c];
%         tot_err_80 = [tot_err_80;err_80];
%         tot_err_89 = [tot_err_89;err_89];
%     end
% 
%     if size(drop_0e, 2) == 24
%         tot_drop_0e = [tot_drop_0e;drop_0e];
%         tot_drop_1f = [tot_drop_1f;drop_1f];
%         tot_drop_37 = [tot_drop_37;drop_37];
%         tot_drop_3c = [tot_drop_3c;drop_3c];
%         tot_drop_4c = [tot_drop_4c;drop_4c];
%         tot_drop_80 = [tot_drop_80;drop_80];
%         tot_drop_89 = [tot_drop_89;drop_89];
%     end
% 
% %     for j = 5:11
% %         error_name = files(j).name;
% %         error = readtable(error_name);
% %         error = error(:,headers);
% %         drops = sum(isnan(table2array(error)));
% %         drops = drops/height(table2array(error));
% %         error = mean(table2array(error), "omitnan");
% %     end
%     cd('..')
% end
% %pixels = readtable('p2_avgerr.csv');
% 
% total_avg_err = [mean(tot_err_0e);mean(tot_err_1f);mean(tot_err_37);mean(tot_err_3c);mean(tot_err_4c);mean(tot_err_80);mean(tot_err_89)];
% total_avg_drop = [mean(tot_drop_0e);mean(tot_drop_1f);mean(tot_drop_37);mean(tot_drop_3c);mean(tot_drop_4c);mean(tot_drop_80);mean(tot_drop_89)];
% 
% total_err = [tot_err_0e;tot_err_1f;tot_err_37;tot_err_3c;tot_err_4c;tot_err_80;tot_err_89];
% total_drop = [tot_drop_0e;tot_drop_1f;tot_drop_37;tot_drop_3c;tot_drop_4c;tot_drop_80;tot_drop_89];
% 
% figure
% plot(colnames, total_avg_err)
% title("Average Reprojection Error")
% ylabel('Reprojection Error (pixels)')
% xlabel('Joints')
% legend(["e3v830e","e3v831f","e3v8337","e3v833c","e3v834c","e3v8380","e3v8389"],"AutoUpdate","off")
% set(gca,'fontsize',14)
% 
% figure
% plot(colnames, total_avg_drop)
% title("Average Dropout")
% ylabel('Dropout %')
% xlabel('Joints')
% legend(["e3v830e","e3v831f","e3v8337","e3v833c","e3v834c","e3v8380","e3v8389"],"AutoUpdate","off")
% set(gca,'fontsize',14)

%% All cameras
% Quality rules:
% Thumb tracking drop out < 0.25 is green
% Thumb tracking drop out < 0.5 is yellow
% No more than 3 finger joints drop out > 0.05 is green
% No more than 3 finger joints drop out > 0.1 is yellow

cd(main_dir)
trials = dir(main_dir);

output = [];
for k = 3:length(trials)
    currT = trials(k).name;
    files = dir(currT);
    cd(currT)

    %% Separating trials
    
    T = readtable('df_3d.csv');
    h = height(T); % Number of rows in table
    start_buffer = 30;
    
    % Calculating different phases
    count = 0;
    stage = 0;
    start_trial = [];
    end_trial = [];
    for i=1:h
        if stage == 0
            if (T.RW_x(i) > (mean(T.RW_x(2:6))+start_buffer)) && (T.RW_y(i) < (mean(T.RW_y(2:6))-start_buffer)) && (T.RW_z(i) < (mean(T.RW_z(2:6))-2*start_buffer))
                count = count + 1;
                if count == 20
                    count = 0;
                    stage = 1;
                    start_trial = [start_trial, i-30];
                end
            else
                count = 0;
            end
        else
            if (T.RW_x(i) < (mean(T.RW_x(2:6))+start_buffer)) && (T.RW_y(i) > (mean(T.RW_y(2:6))-start_buffer)) && (T.RW_z(i) > (mean(T.RW_z(2:6))-2*start_buffer))
                count = count + 1;
                if count == 20
                    count = 0;
                    stage = 0;
                    end_trial = [end_trial, i-30];
                end
            else
                count = 0;
            end
        end
    end
    
    if length(start_trial) < length(end_trial)
        end_trial(end)=[];
    end
    if length(start_trial) > length(end_trial)
        start_trial(end)=[];
    end
    valid = 0;
    if length(start_trial) == 5
        valid = 1;
    elseif length(start_trial) == 10
        valid = 2;
    end

    %% Error
    err = readtable('df_err.csv');
    err = err(:,headers);
    err = table2array(err);

    if valid == 1
        loop = 1:5; 
    elseif valid == 2
        loop = 1:2:10;
    elseif valid == 0
        loop = 1:length(start_trial);
    end

    temp_output = [];
    for j = loop
        err2 = err(start_trial(j):end_trial(j),:);
        drop = sum(isnan(err2));
        high_err = sum(err2 > high_err_thresh);
        drop = (drop+high_err)/height(err2);
        if size(drop, 2) == 24
            count = 0;
            if sum(drop(1:4) < 0.25) == 4
            elseif sum(drop(1:4) < 0.5) == 4
                count = count + 1;
            else
                count = 10;
            end
            if sum(drop(5:20) < 0.05) >= 13
            elseif sum(drop(5:20) < 0.1) >= 13
                count = count + 1;
            else
                count = 10;
            end
    
            if count == 0
                temp_output = [temp_output,"3"];
            elseif count == 10
                temp_output = [temp_output,"1"];
            else
                temp_output = [temp_output,"2"];
            end
    
        end
    end

    total_drop = sum(isnan(err));
    high_err = sum(err > high_err_thresh);
    total_drop = (total_drop+high_err)/height(err);
    total_out = "0";
    if size(total_drop, 2) == 24
        % figure
        % plot(colnames, drop)
        % title("Average Dropout")
        % ylabel('Dropout %')
        % xlabel('Joints')
        % set(gca,'fontsize',14)
        count = 0;
        if sum(total_drop(1:4) < 0.25) == 4
        elseif sum(total_drop(1:4) < 0.5) == 4
            count = count + 1;
        else
            count = 10;
        end
        if sum(total_drop(5:20) < 0.05) >= 13
        elseif sum(total_drop(5:20) < 0.1) >= 13
            count = count + 1;
        else
            count = 10;
        end

        if count == 0
            total_out = "3";
        elseif count == 10
            total_out = "1";
        else
            total_out = "2";
        end

    end

    while length(temp_output) < 5
    temp_output = [temp_output, "0"];
    end
    if length(temp_output) == 5
        output = [output;[currT, [temp_output, total_out]]];
    else
        %output = [output;[currT, "0","0","0","0","0"]];
        output = [output;[currT, [temp_output(1:5), "0"]]];
    end
    cd('..')
end
% green = sum(output(:,2) == "Green");
% yellow = sum(output(:,2) == "Yellow");
% red = sum(output(:,2) == "Red");

temp_out = str2double(output(:,2:6));
summation = [sum(temp_out(:) == 0),sum(temp_out(:) == 1),sum(temp_out(:) == 2),sum(temp_out(:) == 3)];
output = array2table(output, 'VariableNames',{'Activity time', 'Trial 1', 'Trial 2', 'Trial 3', 'Trial 4', 'Trial 5', 'Total'});

%% Previous work

% figure
% bar(colnames, err_0e)
% title("Average Reprojection Error")
% ylabel('Reprojection Error (pixels)')
% xlabel('Joints')
% set(gca,'fontsize',14)
% figure
% bar(colnames, table2array(pixels))
% title("Average Reprojection Error")
% ylabel('Reprojection Error (pixels)')
% xlabel('Joints')
% set(gca,'fontsize',14)
% 
% figure
% plot([1,2,3,3,4,5,6,7,7,8,9],[0,0,0,1,1,1,1,1,0,0,0])
% ylim([0 1.5])
% title("Boxcar filter (window = 5)")
% ylabel('Amplitude')
% xlabel('Frames')
% set(gca,'fontsize',14)
% 
% exp_len =  [45,40,31,28,95,42,24,22,90,44,28,24,88,39,26,24,87,32,22,22,180,290,245];
% %meas_len = [41,31,31,28,98,41,26,23,94,44,30,24,90,38,26,23,85,29,19,18,154,254,251];
% meas_len = [45.1416601819741,38.4582932126603,31.8674249985028,25.8229433417104,96.3862612676963,34.3750626749038,23.3801493065490,...
%     20.9001908983792,94.1444604435723,36.8327002977649,27.7555270789662,23.3106208678708,88.9111249709211,33.6769328324028,...
%     24.1048859838439,20.9439158123359,83.8070199216481,26.4544451069235,16.8428194467693,15.7192462593969,179.673531357003,...
%     256.881737691240,248.140364924508];
% 
% both_lens = [exp_len;meas_len];
% lengths = categorical(["W-TCMC","TCMC-TMCP","TMCP-TIP","TIP-TT","W-IMCP","IMCP-IPIP","IPIP-IDIP","IDIP-IT",...
%     "W-MMCP","MMCP-MPIP","MPIP-MDIP","MDIP-MT","W-RMCP","RMCP-RPIP","RPIP-RDIP","RDIP-RT",...
%     "W-LMCP","LMCP-LPIP","LPIP-LDIP","LDIP-LT","C-S","S-E","E-W"]);
% lengths = reordercats(lengths,string(lengths));
% figure
% hold on
% bar(lengths,both_lens')
% title("Limb lengths")
% ylabel('Length (mm)')
% xlabel('Limbs')
% legend(["Expected","Measured"],"AutoUpdate","off")
% set(gca,'fontsize',14)
% 
% len_err = (abs(exp_len-meas_len)./exp_len)*100;
% figure
% hold on
% bar(lengths,len_err)
% title("Limb length error")
% ylabel('% error')
% xlabel('Limbs')
% set(gca,'fontsize',14)
% 
% figure
% PCS = [11,8,6,9,10,10,12,10,10,9,10];
% cats = categorical(["All","Bathing","Cleaning","Cooking","Digital","Dining","Dressing","General","Mobility","Office","Physical"]);
% bar(cats,PCS)
% title("PCs accounting for over 80% variance in each category")
% ylabel('Number of PCs')
% xlabel('Categories')
% set(gca,'fontsize',12)