close all; clear all; clc;

%src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\data\high_thresh\angles';
src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\categorized_angles';
src_3d = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\categorized';
pc_src = 'C:\Users\czhe0008\Documents\MATLAB\Openpose\PCA_coeffs\';
mean_src = 'C:\Users\czhe0008\Documents\MATLAB\Openpose\PCA_mean\';
pro_src = 'C:\Users\czhe0008\Documents\MATLAB\Openpose\Projection\';
time_src = 'C:\Users\czhe0008\Documents\MATLAB\Openpose\PCA_times\';
unpro_src = 'C:\Users\czhe0008\Documents\MATLAB\Openpose\Unprojected\';
sig_src = 'C:\Users\czhe0008\Documents\MATLAB\Openpose\Significance\';
joints = ["Thumb carpormetacarpal flexion","Thumb metacarpophalangeal flexion","Thumb interphalangeal flexion","Index metacarpophalangeal flexion",...
    "Index proximal interphalangeal flexion","Index distal interphalangeal flexion","Middle metacarpophalangeal flexion","Middle proximal interphalangeal flexion",...
    "Middle distal interphalangeal flexion","Ring metacarpophalangeal flexion","Ring proximal interphalangeal flexion","Ring distal interphalangeal flexion",...
    "Little metacarpophalangeal flexion","Little proximal interphalangeal flexion","Little distal interphalangeal flexion","Elbow flexion",...
    "Index abduction","Middle abduction","Ring abduction","Little abduction","Thumb abduction","Thumb rotation", "Wrist flexion","Wrist Abduction",...
    "Wrist rotation","Shoulder abduction","Shoulder flexion","Shoulder rotation"];
% Fixed joints by swapping shoulder abduction with flexion
colnames = categorical(joints);
colnames = reordercats(colnames,string(colnames));

allpcarr = load("All\pconcat_coeff.mat");
allpcarr = allpcarr.coeff;
all_pro = load("All\pconcat_project.mat");
all_pro = all_pro.score;
all_mean = load("All\pconcat_mean.mat");
all_mean = all_mean.mu;
pctemp = dir(fullfile(pc_src,'*'));
pcfolder = setdiff({pctemp(~[pctemp.isdir]).name},{'.','..'});
meantemp = dir(fullfile(mean_src,'*'));
meanfolder = setdiff({meantemp(~[meantemp.isdir]).name},{'.','..'});
protemp = dir(fullfile(pro_src,'*'));
profolder = setdiff({protemp(~[protemp.isdir]).name},{'.','..'});
timetemp = dir(fullfile(time_src,'*'));
timefolder = setdiff({timetemp(~[timetemp.isdir]).name},{'.','..'});
unprotemp = dir(fullfile(unpro_src,'*'));
unprofolder = setdiff({unprotemp(~[unprotemp.isdir]).name},{'.','..'});
sigtemp = dir(fullfile(sig_src,'*'));
sigfolder = setdiff({sigtemp(~[sigtemp.isdir]).name},{'.','..'});
sig_all = matfile('All\pAll_sig_80.mat');
sig_all = sig_all.significant_80;
maintemp = dir(fullfile(src_3d,'*'));
mainfolder = setdiff({maintemp([maintemp.isdir]).name},{'.','..'});

total_mean_corr = [];
grasp_times = [];
final_start = 0;

% Looping through categories
for i = 1:numel(profolder)
    % Loading files
    category = profolder{i};
    file = fullfile(pro_src,profolder{i});
    fprintf(1, 'Now reading %s\n', profolder{i});
    proTable = readtable(file);
    proarr = table2array(proTable);             % Projected joint angles
    file = fullfile(pc_src,pcfolder{i});
    pcTable = readtable(file);
    pcarr = table2array(pcTable(2:end,:));      % PCA coefficients
    file = fullfile(mean_src,meanfolder{i});
    meanTable = readtable(file);
    meanarr = table2array(meanTable);           % The mean of each variable, needed for reprojection
    file = fullfile(time_src,timefolder{i});
    timeTable = readtable(file);                % Start and end time of each activity
    file = fullfile(unpro_src,unprofolder{i});
    unproTable = readtable(file);
    unproarr = table2array(unproTable);         % Unprojected joint angles
    file = fullfile(sig_src,sigfolder{i});
    sig = matfile(file);
    sig = sig.curr_sig_80;                      % Significant PCs for each participant

    subtemp = dir(fullfile(src_3d,mainfolder{i},'*.csv')); % Getting original 3d joint locations
    subfolder = {subtemp(~[subtemp.isdir]).name};

    %% Finding grasp time for min and max of each PC
    grasp_times = cat(1,grasp_times,timeTable.start_time + timeTable.grasp_time - 1 + final_start);
    final_start = timeTable.end_time(end);
    
    %combined = [pcarr(:,1:sig), allpcarr(:,1:sig_all)];
    %[U,S,V] = svd(combined);
    %temp = proarr * U;

    %% PC grasp category coefficients
    % for j = 1:10 %sig
    %     figure
    %     hold on
    %     bar(colnames,pcarr(:,j))
    %     titlestr = strcat("Coefficients of PC ", int2str(j));
    %     title(titlestr)
    %     xlabel('Joint')
    %     ylabel('Weight (-1 to 1)')
    %     box off
    %     %saveas(gcf,strcat('PCA_images\PC_min_max_test\coefficients_of_PC',int2str(j),'.png'))
    %     saveas(gcf,strcat('ABEC_images\',category,'_PC',int2str(j),'.png'))
    % end

    %% PC signal reconstruction
    % This code visualises how the combinations of PCs reconstruct the original signal

    % unproarr(any(isnan(unproarr), 2), :) = [];
    % temp_corr = [];
    % for j = 8:8         % Looping through joints
    %     for k = 1:28    % Looping through each PC
    %         reproarr = proarr(:,1:k)*pcarr(:,1:k)' + repmat(meanarr,height(proarr),1);
    %         reproarr(any(isnan(reproarr), 2), :) = [];
    %         % figure 
    %         % hold on
    %         % plot((timeTable.start_time(1):timeTable.end_time(1))./40, reproarr(timeTable.start_time(1):timeTable.end_time(1),j).*(180/pi));
    %         % plot((timeTable.start_time(1):timeTable.end_time(1))./40, unproarr(timeTable.start_time(1):timeTable.end_time(1),j).*(180/pi));
    %         % legend('reproject','original')
    %         % xlabel('time (seconds)')
    %         % ylabel('angle (degrees)')
    %         % %ylim([40 180]);
    %         % title(strcat(joints(j)," PC", int2str(k)))
    %         % saveas(gcf,strcat('PCA_images\middle\middle',int2str(k),'.png'))
    %         correlations = corrcoef(unproarr(timeTable.start_time(1):timeTable.end_time(1),j),reproarr(timeTable.start_time(1):timeTable.end_time(1),j));
    %         temp_corr = [temp_corr, correlations(1,2)];
    %     end
    % end
    % figure
    % bar(1:28,temp_corr)
    % title("Correlation between original and reprojected data")
    % ylabel('Correlation')
    % saveas(gcf,strcat('PCA_images\middle\middle_corr.png'))

    %% Reprojection
    %total_corr = [];
    %mean_corr = [];
    unproTable=rmmissing(unproTable);
    unproarr(any(isnan(unproarr), 2), :) = [];
    num_pcs = width(unproarr);
    %unproarr = normalize(unproarr,1); % Seems to be messing things up.

    % Looping through PCs to reproject
    for j = 0:(num_pcs-1)
        reproarr = proarr(:,1:(num_pcs-j))*pcarr(:,1:(num_pcs-j))' + repmat(meanarr,height(proarr),1);
        reproarr(any(isnan(reproarr), 2), :) = [];
        %correlations = corrcoef(unproarr,reproarr);
        %total_corr = [total_corr, correlations(1,2)];


        if j == sig %sig
            trial = 1;
            %% Getting limb lengths
            T = fullfile(src_3d,mainfolder{i},subfolder{trial});
            T = readtable(T);
            original = table2array(T);
            original_table = T(:,[4:63,67:78]);
            original = original(:,[4:63,67:78]);            
            grasp_time = timeTable.grasp_time(trial);
            lengths = pdist([original(grasp_time,1),original(grasp_time,2),original(grasp_time,3);original(grasp_time,70),original(grasp_time,71),original(grasp_time,72)]);
            % Finding the vectors to describe [TCMC, IMCP, MMCP, RMCP, LMCP]
            knuckles = [original(grasp_time,1),original(grasp_time,2),original(grasp_time,3)]-[original(grasp_time,70),original(grasp_time,71),original(grasp_time,72)];
            % Realigning original data
            alignment = original(grasp_time,:) - repmat(original(grasp_time,70:72),1,24);
            for k = 0:22
                if (mod(k+1,4) == 0) && (k < 19)
                    limb = [original(grasp_time,(k+1)*3+1),original(grasp_time,(k+1)*3+2),original(grasp_time,(k+1)*3+3);original(grasp_time,70),original(grasp_time,71),original(grasp_time,72)];
                    knuckles = cat(1,knuckles,[original(grasp_time,(k+1)*3+1),original(grasp_time,(k+1)*3+2),original(grasp_time,(k+1)*3+3)]-[original(grasp_time,70),original(grasp_time,71),original(grasp_time,72)]);
                elseif (mod(k+1,4) == 0) && (k == 19)
                        continue
                else
                    limb = [original(grasp_time,k*3+1),original(grasp_time,k*3+2),original(grasp_time,k*3+3);original(grasp_time,k*3+4),original(grasp_time,k*3+5),original(grasp_time,k*3+6)];
                end
                lengths = cat(1,lengths,pdist(limb));
            end

            %% Adjusting angles
            palm_plane = cross(knuckles(2,:),knuckles(5,:));
            z_axis = [0 0 1];
            % Determine the angle between the vector and the z-axis
            theta = acos(dot(palm_plane, z_axis)/(norm(palm_plane)*norm(z_axis)));
            % Determine the axis of rotation
            axis = cross(palm_plane, z_axis)/norm(cross(palm_plane, z_axis));
            % Construct the rotation matrix using Rodrigues' formula
            K = [0 -axis(3) axis(2); axis(3) 0 -axis(1); -axis(2) axis(1) 0];
            original_R = eye(3) + sin(theta)*K + (1-cos(theta))*K*K;
            % We are rotating the palm_plane by angle theta around vector K
            palm_plane = (-original_R*palm_plane')'.*[1,-1,-1];
            palm_plane = palm_plane./norm(palm_plane);
            knuckles = transpose(-original_R*knuckles');
            knuckles = knuckles.*[1,-1,1];

            K = [0 -palm_plane(3) palm_plane(2); palm_plane(3) 0 -palm_plane(1); -palm_plane(2) palm_plane(1) 0];
            R = eye(3) + sin(unproTable.IMCP_abd(grasp_time)-pi)*K + (1-cos(unproTable.IMCP_abd(grasp_time)-pi))*K*K;
            IPIP = R*knuckles(2,:)'.*(lengths(6)/lengths(5)) + knuckles(2,:)';
            R = eye(3) + sin(unproTable.MMCP_abd(grasp_time)-pi)*K + (1-cos(unproTable.MMCP_abd(grasp_time)-pi))*K*K;
            MPIP = R*knuckles(3,:)'.*(lengths(10)/lengths(9)) + knuckles(3,:)';
            R = eye(3) + sin(unproTable.RMCP_abd(grasp_time)-pi)*K + (1-cos(unproTable.RMCP_abd(grasp_time)-pi))*K*K;
            RPIP = R*knuckles(4,:)'.*(lengths(14)/lengths(13)) + knuckles(4,:)';
            R = eye(3) + sin(unproTable.LMCP_abd(grasp_time)-pi)*K + (1-cos(unproTable.LMCP_abd(grasp_time)-pi))*K*K;
            LPIP = R*knuckles(5,:)'.*(lengths(18)/lengths(17)) + knuckles(5,:)';

            new_plane = cross(knuckles(2,:),palm_plane);
            new_plane = new_plane/norm(new_plane);
            K = [0 -new_plane(3) new_plane(2); new_plane(3) 0 -new_plane(1); -new_plane(2) new_plane(1) 0];
            R = eye(3) + sin(unproTable.IMCP_flex(grasp_time)-pi)*K + (1-cos(unproTable.IMCP_flex(grasp_time)-pi))*K*K;
            IPIP = R*(IPIP-knuckles(2,:)') + knuckles(2,:)';
            R = eye(3) + sin(unproTable.IPIP_flex(grasp_time)-pi)*K + (1-cos(unproTable.IPIP_flex(grasp_time)-pi))*K*K;
            IDIP = R*(IPIP-knuckles(2,:)').*(lengths(7)/lengths(6)) + IPIP;
            R = eye(3) + sin(unproTable.IDIP_flex(grasp_time)-pi)*K + (1-cos(unproTable.IDIP_flex(grasp_time)-pi))*K*K;
            IT = R*(IDIP-IPIP).*(lengths(8)/lengths(7)) + IDIP;

            new_plane = cross(knuckles(3,:),palm_plane);
            new_plane = new_plane/norm(new_plane);
            K = [0 -new_plane(3) new_plane(2); new_plane(3) 0 -new_plane(1); -new_plane(2) new_plane(1) 0];
            R = eye(3) + sin(-(unproTable.MMCP_flex(grasp_time)-pi))*K + (1-cos(-(unproTable.MMCP_flex(grasp_time)-pi)))*K*K;
            MPIP = R*(MPIP-knuckles(3,:)') + knuckles(3,:)';
            R = eye(3) + sin(unproTable.MPIP_flex(grasp_time)-pi)*K + (1-cos(unproTable.MPIP_flex(grasp_time)-pi))*K*K;
            MDIP = R*(MPIP-knuckles(3,:)').*(lengths(11)/lengths(10)) + MPIP;
            R = eye(3) + sin(unproTable.MDIP_flex(grasp_time)-pi)*K + (1-cos(unproTable.MDIP_flex(grasp_time)-pi))*K*K;
            MT = R*(MDIP-MPIP).*(lengths(12)/lengths(11)) + MDIP;

            new_plane = cross(knuckles(4,:),palm_plane);
            new_plane = new_plane/norm(new_plane);
            K = [0 -new_plane(3) new_plane(2); new_plane(3) 0 -new_plane(1); -new_plane(2) new_plane(1) 0];
            R = eye(3) + sin(-(unproTable.RMCP_flex(grasp_time)-pi))*K + (1-cos(-(unproTable.RMCP_flex(grasp_time)-pi)))*K*K;
            RPIP = R*(RPIP-knuckles(4,:)') + knuckles(4,:)';
            R = eye(3) + sin(unproTable.RPIP_flex(grasp_time)-pi)*K + (1-cos(unproTable.RPIP_flex(grasp_time)-pi))*K*K;
            RDIP = R*(RPIP-knuckles(4,:)').*(lengths(15)/lengths(14)) + RPIP;
            R = eye(3) + sin(unproTable.RDIP_flex(grasp_time)-pi)*K + (1-cos(unproTable.RDIP_flex(grasp_time)-pi))*K*K;
            RT = R*(RDIP-RPIP).*(lengths(16)/lengths(15)) + RDIP;

            new_plane = cross(knuckles(5,:),palm_plane);
            new_plane = new_plane/norm(new_plane);
            K = [0 -new_plane(3) new_plane(2); new_plane(3) 0 -new_plane(1); -new_plane(2) new_plane(1) 0];
            R = eye(3) + sin(-(unproTable.LMCP_flex(grasp_time)-pi))*K + (1-cos(-(unproTable.LMCP_flex(grasp_time)-pi)))*K*K;
            LPIP = R*(LPIP-knuckles(5,:)') + knuckles(5,:)';
            R = eye(3) + sin(unproTable.LPIP_flex(grasp_time)-pi)*K + (1-cos(unproTable.LPIP_flex(grasp_time)-pi))*K*K;
            LDIP = R*(LPIP-knuckles(5,:)').*(lengths(19)/lengths(18)) + LPIP;
            R = eye(3) + sin(-(unproTable.LDIP_flex(grasp_time)-pi))*K + (1-cos(-(unproTable.LDIP_flex(grasp_time)-pi)))*K*K;
            LT = R*(LDIP-LPIP).*(lengths(20)/lengths(19)) + LDIP;

            thumb_rot_axis = knuckles(2,:) - knuckles(1,:);
            thumb_plane = cross(knuckles(1,:),knuckles(2,:));
            thumb_plane = (thumb_plane./norm(thumb_plane));
            K = [0 -thumb_plane(3) thumb_plane(2); thumb_plane(3) 0 -thumb_plane(1); -thumb_plane(2) thumb_plane(1) 0];
            R = eye(3) + sin(-unproTable.TMCP_abd(grasp_time))*K + (1-cos(-unproTable.TMCP_abd(grasp_time)))*K*K;
            TMCP = (R*thumb_rot_axis')'.*lengths(2)/norm(thumb_rot_axis);
            vmag = dot(TMCP,thumb_rot_axis)/(norm(thumb_rot_axis)^2);
            v_inline = vmag*thumb_rot_axis;
            v_perp = TMCP - v_inline;
            thumb_rot_axis_unit = thumb_rot_axis/norm(thumb_rot_axis);
            K = [0 -thumb_rot_axis_unit(3) thumb_rot_axis_unit(2); thumb_rot_axis_unit(3) 0 -thumb_rot_axis_unit(1); -thumb_rot_axis_unit(2) thumb_rot_axis_unit(1) 0];
            R = eye(3) + sin(pi/2-unproTable.TCMC_rot(grasp_time))*K + (1-cos(pi/2-unproTable.TCMC_rot(grasp_time)))*K*K;
            TMCP = R*v_perp' + v_inline' + knuckles(1,:)';
            thumb_plane2 = R*thumb_plane';

            K = [0 -thumb_plane2(3) thumb_plane2(2); thumb_plane2(3) 0 -thumb_plane2(1); -thumb_plane2(2) thumb_plane2(1) 0];
            R = eye(3) + sin(pi-unproTable.TMCP_flex(grasp_time))*K + (1-cos(pi-unproTable.TMCP_flex(grasp_time)))*K*K;
            TIP = R*(TMCP-knuckles(1,:)').*(lengths(3)/lengths(2)) + TMCP;
            R = eye(3) + sin(pi-unproTable.TIP_flex(grasp_time))*K + (1-cos(pi-unproTable.TIP_flex(grasp_time)))*K*K;
            TT = R*(TIP-TMCP).*(lengths(4)/lengths(3)) + TIP;

            % Elbow starts in the direction of MMCP to W
            E = -knuckles(3,:).*lengths(23)/lengths(9);
            vmag = dot(E,palm_plane)/(norm(palm_plane)^2);
            v_inline = vmag*palm_plane;
            v_perp = E - v_inline;
            wrist_plane = cross(palm_plane,v_perp);
            wrist_plane = wrist_plane/norm(wrist_plane);
            K = [0 -wrist_plane(3) wrist_plane(2); wrist_plane(3) 0 -wrist_plane(1); -wrist_plane(2) wrist_plane(1) 0];
            R = eye(3) + sin(unproTable.W_flex(grasp_time)-pi)*K + (1-cos(unproTable.W_flex(grasp_time)-pi))*K*K;
            E = R*E';
            K = [0 -palm_plane(3) palm_plane(2); palm_plane(3) 0 -palm_plane(1); -palm_plane(2) palm_plane(1) 0];
            R = eye(3) + sin(unproTable.W_abd(grasp_time))*K + (1-cos(unproTable.W_abd(grasp_time)))*K*K;
            E = transpose(R*E);

            % Edit the below palm back to knuckles(5,:) once I change the
            % python code
            palm_plane2 = cross(knuckles(2,:),knuckles(4,:));
            palm_plane2 = -palm_plane2/norm(palm_plane2);
            
            mid_norm = knuckles(3,:)/norm(knuckles(3,:));
            K = [0 -mid_norm(3) mid_norm(2); mid_norm(3) 0 -mid_norm(1); -mid_norm(2) mid_norm(1) 0];
            R = eye(3) + sin(-unproTable.W_rot(grasp_time))*K + (1-cos(-unproTable.W_rot(grasp_time)))*K*K;
            arm_plane = R*palm_plane2';
            K = [0 -arm_plane(3) arm_plane(2); arm_plane(3) 0 -arm_plane(1); -arm_plane(2) arm_plane(1) 0];
            R = eye(3) + sin(pi-unproTable.RE_flex(grasp_time))*K + (1-cos(pi-unproTable.RE_flex(grasp_time)))*K*K;
            S = transpose(R*E'.*lengths(22)/lengths(23)) + E;
            
            %% Drawing reprojection
            figure
            hold on;
            % New hand
            % palmPlot = plot3([0,palm_plane(1)*30], ...
            %     [0,palm_plane(2)*30], ...
            %     [0,palm_plane(3)*30], ...
            %     '-o', 'MarkerSize',3,'MarkerFaceColor',	'k', 'Color','k');
            ThumbPlot = plot3([0,knuckles(1,1),TMCP(1),TIP(1),TT(1)], ...
                [0,knuckles(1,2),TMCP(2),TIP(2),TT(2)], ...
                [0,knuckles(1,3),TMCP(3),TIP(3),TT(3)], ...
                '-o', 'MarkerSize',3,'MarkerFaceColor',	'r', 'Color','r');
            IndexPlot = plot3([0,knuckles(2,1),IPIP(1),IDIP(1),IT(1)], ...
                [0,knuckles(2,2),IPIP(2),IDIP(2),IT(2)], ...
                [0,knuckles(2,3),IPIP(3),IDIP(3),IT(3)], ...
                '-o', 'MarkerSize',3,'MarkerFaceColor',	'm', 'Color','m');
            MiddlePlot = plot3([0,knuckles(3,1),MPIP(1),MDIP(1),MT(1)], ...
                [0,knuckles(3,2),MPIP(2),MDIP(2),MT(2)], ...
                [0,knuckles(3,3),MPIP(3),MDIP(3),MT(3)], ...
                '-o', 'MarkerSize',3,'MarkerFaceColor',	'b', 'Color','b');
            RingPlot = plot3([0,knuckles(4,1),RPIP(1),RDIP(1),RT(1)], ...
                [0,knuckles(4,2),RPIP(2),RDIP(2),RT(2)], ...
                [0,knuckles(4,3),RPIP(3),RDIP(3),RT(3)], ...
                '-o', 'MarkerSize',3,'MarkerFaceColor',	'c', 'Color','c');
            LittlePlot = plot3([0,knuckles(5,1),LPIP(1),LDIP(1),LT(1)], ...
                [0,knuckles(5,2),LPIP(2),LDIP(2),LT(2)], ...
                [0,knuckles(5,3),LPIP(3),LDIP(3),LT(3)], ...
                '-o', 'MarkerSize',3,'MarkerFaceColor',	'g', 'Color','g');
            RightArmPlot = plot3([0,E(1),S(1)], ...
                [0,E(2),S(2)], ...
                [0,E(3),S(3)], ...
                '-o', 'MarkerSize',3,'MarkerFaceColor',"#D95319", 'Color',"#D95319");
            % xlim([-50 50])
            % ylim([0 100])
            % zlim([-50 50])
            xlabel('x (mm)')
            ylabel('y (mm)')
            zlabel('z (mm)')
            flag = 1;


            %% Drawing original
            alignment_transformed = [];
            x_axis = [1 0 0];
            K = [0 -x_axis(3) x_axis(2); x_axis(3) 0 -x_axis(1); -x_axis(2) x_axis(1) 0];
            R1 = eye(3) + sin(-pi/2)*K + (1-cos(-pi/2))*K*K;
            y_axis = [0 1 0];
            K = [0 -y_axis(3) y_axis(2); y_axis(3) 0 -y_axis(1); -y_axis(2) y_axis(1) 0];
            R2 = eye(3) + sin(-pi/2)*K + (1-cos(-pi/2))*K*K;
            for l = 1:24
                temp = transpose(-original_R*alignment(:,(l*3-2):l*3)').*[1,-1,1];
                temp = transpose(R1*temp');
                temp = transpose(R2*temp');
                alignment_transformed = cat(2,alignment_transformed,temp);
            end
            % figure
            % hold on
            ThumbPlot = plot3([alignment_transformed(12),alignment_transformed(9),alignment_transformed(6),alignment_transformed(3),alignment_transformed(72)], ...
                [alignment_transformed(10),alignment_transformed(7),alignment_transformed(4),alignment_transformed(1),alignment_transformed(70)], ...
                [alignment_transformed(11),alignment_transformed(8),alignment_transformed(5),alignment_transformed(2),alignment_transformed(71)], ...
                '-o', 'MarkerSize',3,'MarkerFaceColor','r', 'Color','k');
            IndexPlot = plot3([alignment_transformed(24),alignment_transformed(21),alignment_transformed(18),alignment_transformed(15),alignment_transformed(72)], ...
                [alignment_transformed(22),alignment_transformed(19),alignment_transformed(16),alignment_transformed(13),alignment_transformed(70)], ...
                [alignment_transformed(23),alignment_transformed(20),alignment_transformed(17),alignment_transformed(14),alignment_transformed(71)], ...
                '-o', 'MarkerSize',3,'MarkerFaceColor','m', 'Color','k');
            MiddlePlot = plot3([alignment_transformed(36),alignment_transformed(33),alignment_transformed(30),alignment_transformed(27),alignment_transformed(72)], ...
                [alignment_transformed(34),alignment_transformed(31),alignment_transformed(28),alignment_transformed(25),alignment_transformed(70)], ...
                [alignment_transformed(35),alignment_transformed(32),alignment_transformed(29),alignment_transformed(26),alignment_transformed(71)], ...
                '-o', 'MarkerSize',3,'MarkerFaceColor','b', 'Color','k');
            RingPlot = plot3([alignment_transformed(48),alignment_transformed(45),alignment_transformed(42),alignment_transformed(39),alignment_transformed(72)], ...
                [alignment_transformed(46),alignment_transformed(43),alignment_transformed(40),alignment_transformed(37),alignment_transformed(70)], ...
                [alignment_transformed(47),alignment_transformed(44),alignment_transformed(41),alignment_transformed(38),alignment_transformed(71)], ...
                '-o', 'MarkerSize',3,'MarkerFaceColor','c', 'Color','k');
            LittlePlot = plot3([alignment_transformed(60),alignment_transformed(57),alignment_transformed(54),alignment_transformed(51),alignment_transformed(72)], ...
                [alignment_transformed(58),alignment_transformed(55),alignment_transformed(52),alignment_transformed(49),alignment_transformed(70)], ...
                [alignment_transformed(59),alignment_transformed(56),alignment_transformed(53),alignment_transformed(50),alignment_transformed(71)], ...
                '-o', 'MarkerSize',3,'MarkerFaceColor','g', 'Color','k');
            RightArmPlot = plot3([alignment_transformed(63),alignment_transformed(66),alignment_transformed(69),alignment_transformed(72)], ...
                 [alignment_transformed(61),alignment_transformed(64),alignment_transformed(67),alignment_transformed(70)], ...
                 [alignment_transformed(62),alignment_transformed(65),alignment_transformed(68),alignment_transformed(71)], ...
                 '-o', 'MarkerSize',3,'MarkerFaceColor',"#D95319", 'Color','k');
            % scatter3(mean(alignment_transformed(:,3:3:end),1,"omitnan"),mean(alignment_transformed(:,1:3:end),1,"omitnan"),mean(alignment_transformed(:,2:3:end),1,"omitnan"))
            xlabel('x (mm)')
            ylabel('y (mm)')
            zlabel('z (mm)')
            % titlestr = strcat("Average hand position of grasp type: ", mainfolder{k});
            % title(titlestr, 'Interpreter', 'none')
            % set(gca, 'XDir','reverse')
            % set(gca, 'ZDir','reverse')
            % set(gca, 'fontsize',14)
            % set(gcf,'position',[200,200,700,600])
            % view([210 40])
            flag = 1;
        end

        % meanreproarr = mean_pro(:,1:(num_pcs-j))*allpcarr(:,1:(num_pcs-j))' + repmat(meanarr,height(mean_pro),1);
        % meanreproarr(any(isnan(meanreproarr), 2), :) = [];
        % correlations = corrcoef(unproarr,meanreproarr);
        % mean_corr = [mean_corr, correlations(1,2)];

    end

    % total_mean_corr = [total_mean_corr;mean_corr];
    % 
    % figure
    % bar(flip(total_corr))
    % titlestr = strcat("Correlation between original and reprojected data for ", category);
    % title(titlestr)
    % xlabel('Number of principal components used in reprojection')
    % ylabel('Correlation')
    % 
    % figure
    % bar(flip(mean_corr))
    % titlestr = strcat("Correlation between original and generalised PC reprojected data for ", category);
    % title(titlestr)
    % xlabel('Number of principal components used in reprojection')
    % ylabel('Correlation')
    
end

%% Also need to find the average joint angle during the moment of grasp for all activities

%% Finding grasp time for min and max of each PC
% Min max are mirrored, so don't need both
for j = 1:10 %sig_all
    % allreproarr = all_pro(:,j)*allpcarr(:,j)';% + repmat(all_mean,height(all_pro),1);
    % antiallreproarr = all_pro(:,j)*allpcarr(:,j)'*-1;% + repmat(all_mean,height(all_pro),1);
    % allrepromean = mean(allreproarr(grasp_times,:),1,"omitnan");
    % antiallrepromean = mean(antiallreproarr(grasp_times,:),1,"omitnan");
    % figure
    % hold on
    % bar(colnames,allrepromean)
    % %plot(colnames,allrepromean)
    % %plot(colnames,antiallrepromean)
    % titlestr = strcat("The effect of PC ", int2str(j));
    % title(titlestr)
    % xlabel('Joint')
    % ylabel('Angle (radians)')
    % %ylim([0 3.15]);
    % box off
    % %saveas(gcf,strcat('PCA_images\PC_min_max_test\reproject_diff_PC',int2str(j),'.png'))
    % saveas(gcf,strcat('grasp_time\images\reproject_diff_PC',int2str(j),'.png'))
    

    % Plotting coefficients
    figure
    hold on
    bar(colnames,allpcarr(:,j))
    titlestr = strcat("Coefficients of PC ", int2str(j));
    title(titlestr)
    xlabel('Joint')
    ylabel('Weight (-1 to 1)')
    box off
    %saveas(gcf,strcat('PCA_images\PC_min_max_test\coefficients_of_PC',int2str(j),'.png'))
    saveas(gcf,strcat('ABEC_images_all\coefficients_of_PC',int2str(j),'.png'))
end


%% Looking at mean projection
% for j = 0:27
%     meanreproarr = mean_pro(:,1:(num_pcs-j))*allpcarr(:,1:(num_pcs-j))' + repmat(meanarr,height(mean_pro),1);
%     meanreproarr(any(isnan(meanreproarr), 2), :) = [];
%     correlations = corrcoef(unproarr,meanreproarr);
%     mean_corr = [mean_corr, correlations(1,2)];
% end

% figure
% hold on
% mean_var = flip(mean(total_mean_corr.^2, 1));
% bar(mean_var, 'w')
% % SEM = std(total_mean_corr.^2,0,1)/sqrt(size(total_mean_corr.^2,1));
% % er = errorbar(mean_var, flip(SEM));
% % er.Color = [0 0 0];
% % er.LineStyle = 'none';
% %title("Average correlation between original and generalised PC reprojected data")
% %title("Reprojection accuracy of generalised principal components")
% xlabel('Number of principal components used in reprojection')
% ylabel('Percentage of variance explained in reprojection')
% set(gca,'fontsize',14, 'TickDir', 'out')
% box off
% ylim([0 1])

% pc_path = strcat(pc_src,coeff);
% 
% pc = readtable(pc_path);
% pc = table2array(pc);
% pc = pc(2:end,:);

% maintemp = dir(fullfile(src,'*'));
% mainfolder = setdiff({maintemp([maintemp.isdir]).name},{'.','..'});
% % for i = 1:numel(mainfolder)
% for i = 5:5
%     if contains(bad_p,mainfolder{i})
%         continue
%     end
%     participantAngles = [];
%     subtemp = dir(fullfile(src,mainfolder{i},'*.csv'));
%     subfolder = {subtemp(~[subtemp.isdir]).name};
%     %for j = 1:numel(subfolder)
%     for j = 33:33
%         file = fullfile(src,mainfolder{i},subfolder{j});
%         fprintf(1, 'Now reading %s\n', subfolder{j});
%         thisTable = readtable(file);
%         participantAngles = [participantAngles;thisTable];
%     end
%     anglesarr = table2array(participantAngles);
%     repo_angles = anglesarr * pc;
% 
%     figure
%     hold on
%     plot(1:height(anglesarr), anglesarr(:,16));
%     legend('flex','abd','rot')
%     ylabel('Angle (radians)')
%     xlabel('timestep')
%     title(mainfolder{i})
%     set(gca,'fontsize',14, 'TickDir', 'out')
%     box off
% 
%     figure
%     hold on
%     plot(1:height(anglesarr), repo_angles(:,5));
%     ylabel('PC activation level')
%     xlabel('timestep')
%     title(mainfolder{i})
%     set(gca,'fontsize',14, 'TickDir', 'out')
%     box off
% end