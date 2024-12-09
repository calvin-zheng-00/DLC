close all; clear all; clc;

%src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\data\high_thresh\angles';
src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\categorized_angles_2';
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
    % unproTable=rmmissing(unproTable);
    % unproarr(any(isnan(unproarr), 2), :) = [];
    num_pcs = width(unproarr);
    %unproarr = normalize(unproarr,1); % Seems to be messing things up.
    if i == 1
        figure
        myVideo = VideoWriter('PCA.avi');
        myVideo.FrameRate = 1;
        % 1-240, 280-520, 600-800, 880-1080, 1160-end
        open(myVideo)
    end

    % Looping through PCs to reproject
    for j = 0:(num_pcs-1)
        reproarr = proarr(:,1:(num_pcs-j))*pcarr(:,1:(num_pcs-j))' + repmat(meanarr,height(proarr),1);
        reproarr(any(isnan(reproarr), 2), :) = [];
        %correlations = corrcoef(unproarr,reproarr);
        %total_corr = [total_corr, correlations(1,2)];


        %if j == num_pcs-1 %sig
        if i == 1% && ismember(j,0:num_pcs-1)
            % trial = 1;
            % unproTable = fullfile(src,mainfolder{i},subfolder{trial});
            % unproTable = readtable(unproTable);
            %% Getting limb lengths
            T = fullfile(src_3d,mainfolder{i},subfolder{1});
            T = readtable(T);
            grasp_time_3d = timeTable.grasp_time(1); %164, 120
            grasp_time = timeTable.end_time(1);
            original = table2array(T);
            % original_table = T(:,[4:63,67:78]);
            original = original(:,[4:63,67:78]);            
            lengths = pdist([original(grasp_time_3d,1),original(grasp_time_3d,2),original(grasp_time_3d,3);original(grasp_time_3d,70),original(grasp_time_3d,71),original(grasp_time_3d,72)]);
            % Finding the vectors to describe [TCMC, IMCP, MMCP, RMCP, LMCP]
            knuckles = [original(grasp_time_3d,1),original(grasp_time_3d,2),original(grasp_time_3d,3)]-[original(grasp_time_3d,70),original(grasp_time_3d,71),original(grasp_time_3d,72)];
            % Realigning original data
            alignment = original(grasp_time_3d,:) - repmat(original(grasp_time_3d,70:72),1,24);
            for k = 0:22
                if (mod(k+1,4) == 0) && (k < 19)
                    limb = [original(grasp_time_3d,(k+1)*3+1),original(grasp_time_3d,(k+1)*3+2),original(grasp_time_3d,(k+1)*3+3);original(grasp_time_3d,70),original(grasp_time_3d,71),original(grasp_time_3d,72)];
                    knuckles = cat(1,knuckles,[original(grasp_time_3d,(k+1)*3+1),original(grasp_time_3d,(k+1)*3+2),original(grasp_time_3d,(k+1)*3+3)]-[original(grasp_time_3d,70),original(grasp_time_3d,71),original(grasp_time_3d,72)]);
                elseif (mod(k+1,4) == 0) && (k == 19)
                        continue
                else
                    limb = [original(grasp_time_3d,k*3+1),original(grasp_time_3d,k*3+2),original(grasp_time_3d,k*3+3);original(grasp_time_3d,k*3+4),original(grasp_time_3d,k*3+5),original(grasp_time_3d,k*3+6)];
                end
                lengths = cat(1,lengths,pdist(limb));
            end
            
            %% Adjusting angles
            palm_plane = -cross(knuckles(2,:),knuckles(5,:));
            z_axis = [0 0 1];
            % Determine the angle between the vector and the z-axis
            theta = acos(dot(palm_plane, z_axis)/(norm(palm_plane)*norm(z_axis)));
            % Determine the axis of rotation
            axis = cross(palm_plane, z_axis)/norm(cross(palm_plane, z_axis));
            % Construct the rotation matrix using Rodrigues' formula
            K = [0 -axis(3) axis(2); axis(3) 0 -axis(1); -axis(2) axis(1) 0];
            R1 = eye(3) + sin(theta)*K + (1-cos(theta))*K*K;
            % We are rotating the palm_plane by angle theta around vector K
            knuckles_t = transpose(R1*knuckles');
            
            y_axis = [0 1 0];
            % Determine the angle between the vector and the z-axis
            theta = acos(dot(knuckles_t(2,:), y_axis)/(norm(knuckles_t(2,:))*norm(y_axis)));
            % Determine the axis of rotation
            axis = cross(knuckles_t(2,:), y_axis)/norm(cross(knuckles_t(2,:), y_axis));
            % Construct the rotation matrix using Rodrigues' formula
            K = [0 -axis(3) axis(2); axis(3) 0 -axis(1); -axis(2) axis(1) 0];
            R2 = eye(3) + sin(theta)*K + (1-cos(theta))*K*K;
            
            original_R = R2*R1;
            % We are rotating the palm_plane by angle theta around vector K
            palm_plane = (original_R*palm_plane')';
            palm_plane = palm_plane./norm(palm_plane);
            knuckles = transpose(original_R*knuckles');
            
            PIP = [];
            DIP = [];
            Tip = [];
            
            for k = 2:5
                m = (k-1)*4+2;
                v2_np = palm_plane.*((knuckles(k,1).*palm_plane(1) + knuckles(k,2).*palm_plane(2) + knuckles(k,3).*palm_plane(3)));   % Vector component normal to the palm plane
                v2_p = [knuckles(k,1)-v2_np(1),knuckles(k,2)-v2_np(2),knuckles(k,3)-v2_np(3)];                                       % Vector component on the palm plane
                v2mag = sqrt(v2_p(1).^2 + v2_p(2).^2 + v2_p(3).^2);
                v2norm = v2_p./v2mag;
                
                x_axis = [1 0 0];
                z_axis = [0 0 1];
                % Determine the angle between the vector and the x-axis
                theta = acos(dot(x_axis, v2norm)/(norm(x_axis)*norm(v2norm)));
                % Determine the axis of rotation
                axis = cross(x_axis, v2norm)/norm(cross(x_axis, v2norm));
                % Construct the rotation matrix using Rodrigues' formula
                K = [0 -axis(3) axis(2); axis(3) 0 -axis(1); -axis(2) axis(1) 0];
                R1 = eye(3) + sin(theta)*K + (1-cos(theta))*K*K;
                new_z_axis = (R1*z_axis')';
                
                % Determine the angle between the vector and the y-axis
                theta = acos(dot(new_z_axis, palm_plane)/(norm(new_z_axis)*norm(palm_plane)));
                % Determine the axis of rotation
                axis = cross(new_z_axis, palm_plane)/norm(cross(new_z_axis, palm_plane));
                % Construct the rotation matrix using Rodrigues' formula
                K = [0 -axis(3) axis(2); axis(3) 0 -axis(1); -axis(2) axis(1) 0];
                R2 = eye(3) + sin(theta)*K + (1-cos(theta))*K*K;
                R_knuckle = R2*R1;
                [x,y,z] = sph2cart(reproarr(grasp_time,k+16),reproarr(grasp_time,(k-1)*3+1),lengths(m));
                PIP_temp = [x,y,z];
                PIP_temp = (R_knuckle*PIP_temp')';
            
                finger_axis_np = palm_plane.*(PIP_temp(1).*palm_plane(1) + PIP_temp(2).*palm_plane(2) + PIP_temp(3).*palm_plane(3));  % Vector component normal to the thumb plane
                finger_axis_p = [PIP_temp(1)-finger_axis_np(1),PIP_temp(2)-finger_axis_np(2),PIP_temp(3)-finger_axis_np(3)];                                    % Vector component on the thumb plane
                finger_axismag = sqrt(finger_axis_p(1).^2 + finger_axis_p(2).^2 + finger_axis_p(3).^2);
                finger_axis = finger_axis_p./finger_axismag;
                K = [0 -palm_plane(3) palm_plane(2); palm_plane(3) 0 -palm_plane(1); -palm_plane(2) palm_plane(1) 0];
                R = eye(3) + sin(-pi/2)*K + (1-cos(-pi/2))*K*K;
                finger_axis = (R*finger_axis')';
                PIP = [PIP;PIP_temp + knuckles(k,:)];
                K = [0 -finger_axis(3) finger_axis(2); finger_axis(3) 0 -finger_axis(1); -finger_axis(2) finger_axis(1) 0];
                R = eye(3) + sin(reproarr(grasp_time,(k-1)*3+2))*K + (1-cos(reproarr(grasp_time,(k-1)*3+2)))*K*K;
                DIP = [DIP;(R*(PIP(k-1,:)-knuckles(k,:))')'.*(lengths(m+1)/lengths(m)) + PIP(k-1,:)];
                R = eye(3) + sin(reproarr(grasp_time,(k-1)*3+3))*K + (1-cos(reproarr(grasp_time,(k-1)*3+3)))*K*K;
                Tip = [Tip;(R*(DIP(k-1,:)-PIP(k-1,:))')'.*(lengths(m+2)/lengths(m+1)) + DIP(k-1,:)];
            end
            
            %% Fix definition
            %TCMC_flex is the angle between TCMC-IMCP and TCMC-TMCP along the thumb plane (RW-TCMC-IMCP)
            %TMCP_abd is the angle between TCMC-IMCP and TCMC-TMCP out of the thumb plane (along the tabd plane described by the cross product between the thumb pane normal and TCMC-IMCP)
            %TCMC_rot is the angle between the component of the thumb rotation plane (TCMC-TMCP-TIP) that's orthogonal to TCMC-IMCP, and the tabd plane
            thumb_rot_axis = knuckles(2,:) - knuckles(1,:);
            thumb_rot_axis = thumb_rot_axis./norm(thumb_rot_axis);
            thumb_plane = cross(knuckles(1,:),knuckles(2,:));
            thumb_plane = thumb_plane./norm(thumb_plane);
            thumb_abd_axis = cross(thumb_plane,thumb_rot_axis);
            thumb_abd_axis = (thumb_abd_axis./norm(thumb_abd_axis));
            
            x_axis = [1 0 0];
            y_axis = [0 1 0];
            % Determine the angle between the vector and the z-axis
            theta = acos(dot(x_axis, thumb_rot_axis)/(norm(x_axis)*norm(thumb_rot_axis)));
            % Determine the axis of rotation
            axis = cross(x_axis, thumb_rot_axis)/norm(cross(x_axis, thumb_rot_axis));
            % Construct the rotation matrix using Rodrigues' formula
            K = [0 -axis(3) axis(2); axis(3) 0 -axis(1); -axis(2) axis(1) 0];
            R1 = eye(3) + sin(theta)*K + (1-cos(theta))*K*K;
            % We are rotating the palm_plane by angle theta around vector K
            new_y_axis = (R1*y_axis')';
            
            % Determine the angle between the vector and the z-axis
            theta = acos(dot(new_y_axis, thumb_abd_axis)/(norm(new_y_axis)*norm(thumb_abd_axis)));
            % Determine the axis of rotation
            axis = cross(new_y_axis, thumb_abd_axis)/norm(cross(new_y_axis, thumb_abd_axis));
            % Construct the rotation matrix using Rodrigues' formula
            K = [0 -axis(3) axis(2); axis(3) 0 -axis(1); -axis(2) axis(1) 0];
            R2 = eye(3) + sin(theta)*K + (1-cos(theta))*K*K;
            % We are rotating the palm_plane by angle theta around vector K
            R_thumb = R2*R1;
            
            [x,y,z] = sph2cart(reproarr(grasp_time,1),reproarr(grasp_time,16),lengths(2));
            TMCP = [x,y,z];
            TMCP = (R_thumb*TMCP')';
            TMCP_mag = sqrt(TMCP(:,1).^2 + TMCP(:,2).^2 + TMCP(:,3).^2);
            
            thumb_plane_np = TMCP.*((thumb_plane(:,1).*TMCP(:,1) + thumb_plane(:,2).*TMCP(:,2) + thumb_plane(:,3).*TMCP(:,3))./TMCP_mag.^2);  % thumb plane component along v1
            thumb_plane_p = thumb_plane - thumb_plane_np;                               % thumb plane component orthogonal to v1
            thumb_plane2 = thumb_plane_p./norm(thumb_plane_p);
            TMCP_norm = TMCP./TMCP_mag;
            K = [0 -TMCP_norm(3) TMCP_norm(2); TMCP_norm(3) 0 -TMCP_norm(1); -TMCP_norm(2) TMCP_norm(1) 0];
            R = eye(3) + sin(reproarr(grasp_time,17))*K + (1-cos(reproarr(grasp_time,17)))*K*K;
            thumb_plane2 = (R*thumb_plane2')';
            thumb_plane2 = -thumb_plane2./norm(thumb_plane2);
            TMCP = TMCP + knuckles(1,:);
            
            K = [0 -thumb_plane2(3) thumb_plane2(2); thumb_plane2(3) 0 -thumb_plane2(1); -thumb_plane2(2) thumb_plane2(1) 0];
            R = eye(3) + sin(reproarr(grasp_time,2))*K + (1-cos(reproarr(grasp_time,2)))*K*K;
            TIP = (R*(TMCP-knuckles(1,:))')'.*(lengths(3)/lengths(2)) + TMCP;
            R = eye(3) + sin(reproarr(grasp_time,3))*K + (1-cos(reproarr(grasp_time,3)))*K*K;
            TT = (R*(TIP-TMCP)')'.*(lengths(4)/lengths(3)) + TIP;
            
            % Elbow starts in the direction of MMCP to W
            x_axis = [1 0 0];
            z_axis = [0 0 1];
            wrist_forward_np = palm_plane.*((knuckles(3,1).*palm_plane(1) + knuckles(3,2).*palm_plane(2) + knuckles(3,3).*palm_plane(3)));   % Vector component normal to the palm plane
            wrist_forward_p = [knuckles(3,1)-wrist_forward_np(1),knuckles(3,2)-wrist_forward_np(2),knuckles(3,3)-wrist_forward_np(3)];                                       % Vector component on the palm plane
            wrist_forwardmag = sqrt(wrist_forward_p(1).^2 + wrist_forward_p(2).^2 + wrist_forward_p(3).^2);
            wrist_forwardnorm = wrist_forward_p./wrist_forwardmag;
            wrist_fe_axis = cross(palm_plane, wrist_forwardnorm,2);
            % Determine the angle between the vector and the x-axis
            theta = acos(dot(x_axis, wrist_fe_axis)/(norm(x_axis)*norm(wrist_fe_axis)));
            % Determine the axis of rotation
            axis = cross(x_axis, wrist_fe_axis)/norm(cross(x_axis, wrist_fe_axis));
            % Construct the rotation matrix using Rodrigues' formula
            K = [0 -axis(3) axis(2); axis(3) 0 -axis(1); -axis(2) axis(1) 0];
            R1 = eye(3) + sin(theta)*K + (1-cos(theta))*K*K;
            new_z_axis = (R1*z_axis')';
            
            % Determine the angle between the vector and the y-axis
            theta = acos(dot(new_z_axis, palm_plane)/(norm(new_z_axis)*norm(palm_plane)));
            % Determine the axis of rotation
            axis = cross(new_z_axis, palm_plane)/norm(cross(new_z_axis, palm_plane));
            % Construct the rotation matrix using Rodrigues' formula
            K = [0 -axis(3) axis(2); axis(3) 0 -axis(1); -axis(2) axis(1) 0];
            R2 = eye(3) + sin(theta)*K + (1-cos(theta))*K*K;
            R_Elbow = R2*R1;
            [x,y,z] = sph2cart(reproarr(grasp_time,22),reproarr(grasp_time,23),lengths(23));
            E = [x,y,z];
            E = (R_Elbow*E')';
            
            v1_np = palm_plane.*(knuckles(3,1).*palm_plane(:,1) + knuckles(3,2).*palm_plane(:,2) + knuckles(3,3).*palm_plane(:,3));   % Vector component normal to the palm plane
            v1_p = [knuckles(3,1)-v1_np(:,1),knuckles(3,2)-v1_np(:,2),knuckles(3,3)-v1_np(:,3)];                                       % Vector component on the palm plane
            v1mag = sqrt(v1_p(:,1).^2 + v1_p(:,2).^2 + v1_p(:,3).^2);
            wrist_forward = v1_p./v1mag;
            wrist_fe_axis = cross(palm_plane,wrist_forward,2);  % Wrist flexion axis
            wrist_fe_axis_mag = sqrt(wrist_fe_axis(:,1).^2 + wrist_fe_axis(:,2).^2 + wrist_fe_axis(:,3).^2);
            wrist_fe_axis = wrist_fe_axis./wrist_fe_axis_mag;
            forearm_np = E.*((wrist_fe_axis(:,1).*E(:,1) + wrist_fe_axis(:,2).*E(:,2) + wrist_fe_axis(:,3).*E(:,3))./(lengths(23).^2));   % Vector component normal to the palm plane
            forearm_p = [wrist_fe_axis(:,1)-forearm_np(:,1),wrist_fe_axis(:,2)-forearm_np(:,2),wrist_fe_axis(:,3)-forearm_np(:,3)];        % Vector component on the palm plane
            forearmmag = sqrt(forearm_p(:,1).^2 + forearm_p(:,2).^2 + forearm_p(:,3).^2);
            forearm_ref = forearm_p./forearmmag;
            
            E_rot = E/lengths(23);
            K = [0 -E_rot(3) E_rot(2); E_rot(3) 0 -E_rot(1); -E_rot(2) E_rot(1) 0];
            R = eye(3) + sin(-reproarr(grasp_time,24))*K + (1-cos(-reproarr(grasp_time,24)))*K*K;
            E_axis = (R*forearm_ref')';
            K = [0 -E_axis(3) E_axis(2); E_axis(3) 0 -E_axis(1); -E_axis(2) E_axis(1) 0];
            R = eye(3) + sin(reproarr(grasp_time,25))*K + (1-cos(reproarr(grasp_time,25)))*K*K;
            S = (R*E')'.*lengths(22)./lengths(23) + E;
            
            s_third_axis = cross(E_axis,(S-E),2);
            x_axis = [1 0 0];
            z_axis = [0 0 1];
            % Determine the angle between the vector and the z-axis
            theta = acos(dot(x_axis, S-E)/(norm(x_axis)*norm(S-E)));
            % Determine the axis of rotation
            axis = cross(x_axis, S-E)/norm(cross(x_axis, S-E));
            % Construct the rotation matrix using Rodrigues' formula
            K = [0 -axis(3) axis(2); axis(3) 0 -axis(1); -axis(2) axis(1) 0];
            R1 = eye(3) + sin(theta)*K + (1-cos(theta))*K*K;
            % We are rotating the palm_plane by angle theta around vector K
            new_z_axis = (R1*z_axis')';
            
            % Determine the angle between the vector and the z-axis
            theta = acos(dot(new_z_axis, E_axis)/(norm(new_z_axis)*norm(E_axis)));
            % Determine the axis of rotation
            axis = cross(new_z_axis, E_axis)/norm(cross(new_z_axis, E_axis));
            % Construct the rotation matrix using Rodrigues' formula
            K = [0 -axis(3) axis(2); axis(3) 0 -axis(1); -axis(2) axis(1) 0];
            R2 = eye(3) + sin(theta)*K + (1-cos(theta))*K*K;
            % We are rotating the palm_plane by angle theta around vector K
            R_Shoulder = R2*R1;
            
            [x,y,z] = sph2cart(reproarr(grasp_time,26),reproarr(grasp_time,27),lengths(21));
            C = [x,y,z];
            C = (R_Shoulder*C')' + S;
            
            chest_axis = C-S;
            chest_axis = chest_axis/norm(chest_axis);
            shoulder_norm = cross(chest_axis,S-E,2);
            shoulder_normmag = sqrt(shoulder_norm(:,1).^2 + shoulder_norm(:,2).^2 + shoulder_norm(:,3).^2);
            unit_shoulder_norm = shoulder_norm./shoulder_normmag;
            K = [0 -chest_axis(3) chest_axis(2); chest_axis(3) 0 -chest_axis(1); -chest_axis(2) chest_axis(1) 0];
            R = eye(3) + sin(reproarr(grasp_time,28))*K + (1-cos(reproarr(grasp_time,28)))*K*K;
            N = (R*unit_shoulder_norm')'*-100+C;
            
            %% Drawing reprojection
            % New hand
            ThumbPlot = plot3([0,knuckles(1,1),TMCP(1),TIP(1),TT(1)], ...
                [0,knuckles(1,2),TMCP(2),TIP(2),TT(2)], ...
                [0,knuckles(1,3),TMCP(3),TIP(3),TT(3)], ...
                '-o', 'MarkerSize',3,'MarkerFaceColor',	'r', 'Color','r');
            hold on;
            IndexPlot = plot3([0,knuckles(2,1),PIP(1,1),DIP(1,1),Tip(1,1)], ...
                [0,knuckles(2,2),PIP(1,2),DIP(1,2),Tip(1,2)], ...
                [0,knuckles(2,3),PIP(1,3),DIP(1,3),Tip(1,3)], ...
                '-o', 'MarkerSize',3,'MarkerFaceColor',	'm', 'Color','m');
            MiddlePlot = plot3([0,knuckles(3,1),PIP(2,1),DIP(2,1),Tip(2,1)], ...
                [0,knuckles(3,2),PIP(2,2),DIP(2,2),Tip(2,2)], ...
                [0,knuckles(3,3),PIP(2,3),DIP(2,3),Tip(2,3)], ...
                '-o', 'MarkerSize',3,'MarkerFaceColor',	'b', 'Color','b');
            RingPlot = plot3([0,knuckles(4,1),PIP(3,1),DIP(3,1),Tip(3,1)], ...
                [0,knuckles(4,2),PIP(3,2),DIP(3,2),Tip(3,2)], ...
                [0,knuckles(4,3),PIP(3,3),DIP(3,3),Tip(3,3)], ...
                '-o', 'MarkerSize',3,'MarkerFaceColor',	'c', 'Color','c');
            LittlePlot = plot3([0,knuckles(5,1),PIP(4,1),DIP(4,1),Tip(4,1)], ...
                [0,knuckles(5,2),PIP(4,2),DIP(4,2),Tip(4,2)], ...
                [0,knuckles(5,3),PIP(4,3),DIP(4,3),Tip(4,3)], ...
                '-o', 'MarkerSize',3,'MarkerFaceColor',	'g', 'Color','g');
            RightArmPlot = plot3([0,E(1),S(1),C(1)], ...
                [0,E(2),S(2),C(2)], ...
                [0,E(3),S(3),C(3)], ...
                '-o', 'MarkerSize',3,'MarkerFaceColor',"#D95319", 'Color',"#D95319");
            ReferencePlot = plot3([C(1),N(1)], ...
                 [C(2),N(2)], ...
                 [C(3),N(3)], ...
                 '-o', 'MarkerSize',3,'MarkerFaceColor',"#7E2F8E", 'Color',"#7E2F8E");
            
            xlabel('x (mm)')
            ylabel('y (mm)')
            zlabel('z (mm)')
            
            %% Drawing original
            alignment_transformed = [];
            for l = 1:24
                temp = transpose(original_R*alignment(:,(l*3-2):l*3)');%.*[1,-1,1];
                alignment_transformed = cat(2,alignment_transformed,temp);
            end
            Nose = [T.N_x(grasp_time_3d),T.N_y(grasp_time_3d),T.N_z(grasp_time_3d)];
            Nose = Nose - original(grasp_time_3d,70:72);
            Nose = (original_R*Nose')';
            
            ThumbPlot = plot3([alignment_transformed(10),alignment_transformed(7),alignment_transformed(4),alignment_transformed(1),alignment_transformed(70)], ...
                [alignment_transformed(11),alignment_transformed(8),alignment_transformed(5),alignment_transformed(2),alignment_transformed(71)], ...
                [alignment_transformed(12),alignment_transformed(9),alignment_transformed(6),alignment_transformed(3),alignment_transformed(72)], ...
                '-o', 'MarkerSize',3,'MarkerFaceColor','r', 'Color','k');
            IndexPlot = plot3([alignment_transformed(22),alignment_transformed(19),alignment_transformed(16),alignment_transformed(13),alignment_transformed(70)], ...
                [alignment_transformed(23),alignment_transformed(20),alignment_transformed(17),alignment_transformed(14),alignment_transformed(71)], ...
                [alignment_transformed(24),alignment_transformed(21),alignment_transformed(18),alignment_transformed(15),alignment_transformed(72)], ...
                '-o', 'MarkerSize',3,'MarkerFaceColor','m', 'Color','k');
            MiddlePlot = plot3([alignment_transformed(34),alignment_transformed(31),alignment_transformed(28),alignment_transformed(25),alignment_transformed(70)], ...
                [alignment_transformed(35),alignment_transformed(32),alignment_transformed(29),alignment_transformed(26),alignment_transformed(71)], ...
                [alignment_transformed(36),alignment_transformed(33),alignment_transformed(30),alignment_transformed(27),alignment_transformed(72)], ...
                '-o', 'MarkerSize',3,'MarkerFaceColor','b', 'Color','k');
            RingPlot = plot3([alignment_transformed(46),alignment_transformed(43),alignment_transformed(40),alignment_transformed(37),alignment_transformed(70)], ...
                [alignment_transformed(47),alignment_transformed(44),alignment_transformed(41),alignment_transformed(38),alignment_transformed(71)], ...
                [alignment_transformed(48),alignment_transformed(45),alignment_transformed(42),alignment_transformed(39),alignment_transformed(72)], ...
                '-o', 'MarkerSize',3,'MarkerFaceColor','c', 'Color','k');
            LittlePlot = plot3([alignment_transformed(58),alignment_transformed(55),alignment_transformed(52),alignment_transformed(49),alignment_transformed(70)], ...
                [alignment_transformed(59),alignment_transformed(56),alignment_transformed(53),alignment_transformed(50),alignment_transformed(71)], ...
                [alignment_transformed(60),alignment_transformed(57),alignment_transformed(54),alignment_transformed(51),alignment_transformed(72)], ...
                '-o', 'MarkerSize',3,'MarkerFaceColor','g', 'Color','k');
            RightArmPlot = plot3([Nose(1),alignment_transformed(61),alignment_transformed(64),alignment_transformed(67),alignment_transformed(70)], ...
                 [Nose(2),alignment_transformed(62),alignment_transformed(65),alignment_transformed(68),alignment_transformed(71)], ...
                 [Nose(3),alignment_transformed(63),alignment_transformed(66),alignment_transformed(69),alignment_transformed(72)], ...
                 '-o', 'MarkerSize',3,'MarkerFaceColor',"#D95319", 'Color','k');
            ReferencePlot = plot3([Nose(1),alignment_transformed(61)], ...
                 [Nose(2),alignment_transformed(62)], ...
                 [Nose(3),alignment_transformed(63)], ...
                 '-o', 'MarkerSize',3,'MarkerFaceColor',"#7E2F8E", 'Color','k');
            hold off;
            xlabel('x (mm)')
            ylabel('y (mm)')
            zlabel('z (mm)')
            xlim([-500 200])
            ylim([-500 200])
            zlim([-200 500])
            titlestr = profolder{i};
            title(titlestr(1:end-4), 'Interpreter', 'none');
            % titlestr = strcat("Average hand position of grasp type: ", mainfolder{k});
            % title(titlestr, 'Interpreter', 'none')
            % set(gca, 'XDir','reverse')
            % set(gca, 'ZDir','reverse')
            % set(gca, 'fontsize',14)
            % set(gcf,'position',[200,200,700,600])
            view([190 20])
            % end
            saveas(gcf,'PCA.png')
            frame = getframe(gcf);
            writeVideo(myVideo, frame);
        end

        % meanreproarr = mean_pro(:,1:(num_pcs-j))*allpcarr(:,1:(num_pcs-j))' + repmat(meanarr,height(mean_pro),1);
        % meanreproarr(any(isnan(meanreproarr), 2), :) = [];
        % correlations = corrcoef(unproarr,meanreproarr);
        % mean_corr = [mean_corr, correlations(1,2)];

    end
    if i == 1
        close(myVideo)
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
    saveas(gcf,strcat('PCA_images\PC_min_max_test\coefficients_of_PC',int2str(j),'.png'))
    %saveas(gcf,strcat('ABEC_images_all\coefficients_of_PC',int2str(j),'.png'))
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