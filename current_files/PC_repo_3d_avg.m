% This file reprojects the PCs back into 3d position data and compares the
% performance of the principal components to the original 3d data at the
% moment of grasp.

close all; clear all; clc;

% src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\data\high_thresh\angles';
% src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\categorized_angles';
src = 'C:\Users\czhe0008\Documents\MATLAB\Openpose\archive_4cats';
src_3d = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\categorized_2';
pc_src = fullfile(src,'PCA_coeffs\');
mean_src = fullfile(src,'PCA_mean\');
pro_src = fullfile(src,'Projection\');
time_src = fullfile(src,'PCA_times\');
unpro_src = fullfile(src,'Unprojected\');
sig_src = fullfile(src,'Significance\');
comp_src = fullfile(src,'Comparison\');
joints = ["Thumb carpormetacarpal flexion","Thumb metacarpophalangeal flexion","Thumb interphalangeal flexion","Index metacarpophalangeal flexion",...
    "Index proximal interphalangeal flexion","Index distal interphalangeal flexion","Middle metacarpophalangeal flexion","Middle proximal interphalangeal flexion",...
    "Middle distal interphalangeal flexion","Ring metacarpophalangeal flexion","Ring proximal interphalangeal flexion","Ring distal interphalangeal flexion",...
    "Little metacarpophalangeal flexion","Little proximal interphalangeal flexion","Little distal interphalangeal flexion","Elbow flexion",...
    "Index abduction","Middle abduction","Ring abduction","Little abduction","Thumb abduction","Thumb rotation", "Wrist flexion","Wrist Abduction",...
    "Wrist rotation","Shoulder abduction","Shoulder flexion","Shoulder rotation"];
% Fixed joints by swapping shoulder abduction with flexion
colnames = categorical(joints);
colnames = reordercats(colnames,string(colnames));
colnames3d = {'TCMC','TMCP','TIP','TT','IMCP','IPIP','IDIP','IT','MMCP','MPIP','MDIP','MT','RMCP','RPIP','RDIP','RT','LMCP','LPIP','LDIP','LT','RE','RS','C','N'};

allpcarr = load("archive_4cats\All\pconcat_coeff.mat");
allpcarr = allpcarr.coeff;
allproarr = load("archive_4cats\All\pconcat_project.mat");
allproarr = allproarr.score;
allmeanarr = load("archive_4cats\All\pconcat_mean.mat");
allmeanarr = allmeanarr.mu;
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
comptemp = dir(fullfile(comp_src,'*'));
compfolder = setdiff({comptemp(~[comptemp.isdir]).name},{'.','..'});
sigtemp = dir(fullfile(sig_src,'*'));
sigfolder = setdiff({sigtemp(~[sigtemp.isdir]).name},{'.','..'});
sig_all = matfile('archive_4cats\All\pAll_sig_80.mat');
sig_all = sig_all.significant_80;
maintemp = dir(fullfile(src_3d,'*'));
mainfolder = setdiff({maintemp([maintemp.isdir]).name},{'.','..'});

focus_cat = 1;
avgerrcomp = [];
err_rangecomp = [];
avgerrcomp2 = [];
err_rangecomp2 = [];

opt = [2,4,22,25,27;2,4,10,25,27;2,14,19,20,25];

for a = 1:2
    avgerrcomp_t = [];
    err_rangecomp_t = [];
    repo_cp = 0;
    % Looping through categories
    for i = 1:numel(profolder)
        % Loading files
        category = profolder{i};
        file = fullfile(pro_src,profolder{i});
        % fprintf(1, 'Now reading %s\n', profolder{i});
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
        file = fullfile(comp_src,compfolder{i});
        compTable = readtable(file);                % 3d data generated from angle reprojection to compare PCs with

        if i < focus_cat
            reproarr(1:timeTable.end_time(end),:) = [];
            continue
        end
        if i > focus_cat
            break
        end
    
        %% Edit which PCs are being reprojected
        
        if a == 1
            reproarr = allmeanarr;
        elseif a == 2
            reproarr = (pi/2)*allpcarr(:,4)' + repmat(allmeanarr,height(allproarr),1);
        end

        total_error = [];
        subtemp = dir(fullfile(src_3d,mainfolder{i},'*.csv'));
        subfolder = {subtemp(~[subtemp.isdir]).name};
        % rand_order = 50;
        for j = 50:50
            %% Getting limb lengths
            grasp_time_3d = timeTable.grasp_time(j);
            grasp_time = 1;
            T = fullfile(src_3d,mainfolder{i},subfolder{j});
            T = readtable(T);
            original = table2array(T);
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

            %% Chest centred
            % Shoulder forward flexion: arm pointing down is 0, arm
            % pointing forward goes towards negative.
            % 2 plans, either match up the N values, or compare shoulder
            % flexion angle that's orthogonal to C-S.

            repo_c = C-S;
            x_axis = [compTable.C_x(j),compTable.C_y(j),compTable.C_z(j)]-[compTable.RS_x(j),compTable.RS_y(j),compTable.RS_z(j)];
            % Determine the angle between the vector and the x-axis
            theta = acos(dot(repo_c, x_axis)/(norm(repo_c)*norm(x_axis)));
            % Determine the axis of rotation
            axis = cross(repo_c, x_axis)/norm(cross(repo_c, x_axis));
            % Construct the rotation matrix using Rodrigues' formula
            K = [0 -axis(3) axis(2); axis(3) 0 -axis(1); -axis(2) axis(1) 0];
            R1 = eye(3) + sin(theta)*K + (1-cos(theta))*K*K;

            repo_n_temp = (R1*(N-C)')';
            comp_n = [compTable.N_x(j),compTable.N_y(j),compTable.N_z(j)]-[compTable.C_x(j),compTable.C_y(j),compTable.C_z(j)];
            theta = acos(dot(repo_n_temp, comp_n)/(norm(repo_n_temp)*norm(comp_n)));
            % Determine the axis of rotation
            % axis = x_axis/norm(x_axis);
            axis = (R1*(-repo_c)')'/norm((R1*(-repo_c)')');
            % Construct the rotation matrix using Rodrigues' formula
            K = [0 -axis(3) axis(2); axis(3) 0 -axis(1); -axis(2) axis(1) 0];
            R2 = eye(3) + sin(theta)*K + (1-cos(theta))*K*K;
            R = R2*R1;

            repo_s = (R*(-repo_c)')';
            repo_e = (R*(E-S)')'+repo_s;
            repo_w = (R*(-E)')'+repo_e;
            repo_n = (R*(N-C)')';

            repo_tcmc = (R*(knuckles(1,:))')'+repo_w;
            repo_tmcp = (R*(TMCP-knuckles(1,:))')'+repo_tcmc;
            repo_tip = (R*(TIP-TMCP)')'+repo_tmcp;
            repo_tt = (R*(TT-TIP)')'+repo_tip;

            repo_imcp = (R*(knuckles(2,:))')'+repo_w;
            repo_ipip = (R*(PIP(1,:)-knuckles(2,:))')'+repo_imcp;
            repo_idip = (R*(DIP(1,:)-PIP(1,:))')'+repo_ipip;
            repo_it = (R*(Tip(1,:)-DIP(1,:))')'+repo_idip;

            repo_mmcp = (R*(knuckles(3,:))')'+repo_w;
            repo_mpip = (R*(PIP(2,:)-knuckles(3,:))')'+repo_mmcp;
            repo_mdip = (R*(DIP(2,:)-PIP(2,:))')'+repo_mpip;
            repo_mt = (R*(Tip(2,:)-DIP(2,:))')'+repo_mdip;

            repo_rmcp = (R*(knuckles(4,:))')'+repo_w;
            repo_rpip = (R*(PIP(3,:)-knuckles(4,:))')'+repo_rmcp;
            repo_rdip = (R*(DIP(3,:)-PIP(3,:))')'+repo_rpip;
            repo_rt = (R*(Tip(3,:)-DIP(3,:))')'+repo_rdip;

            repo_lmcp = (R*(knuckles(5,:))')'+repo_w;
            repo_lpip = (R*(PIP(4,:)-knuckles(5,:))')'+repo_lmcp;
            repo_ldip = (R*(DIP(4,:)-PIP(4,:))')'+repo_lpip;
            repo_lt = (R*(Tip(4,:)-DIP(4,:))')'+repo_ldip;

            new_og = -[compTable.C_x(j),compTable.C_y(j),compTable.C_z(j)];
            tempcompTable = compTable(j,:) + repmat(new_og,1,24);
            tempcompTable.("RW_x") = new_og(1);
            tempcompTable.("RW_y") = new_og(2);
            tempcompTable.("RW_z") = new_og(3);
    
            %% Drawing reprojection
            if a == 1
                figure
                hold on;
                % New hand
                % ThumbPlotre = plot3([repo_w(1),repo_tcmc(1),repo_tmcp(1),repo_tip(1),repo_tt(1)], ...
                %     [repo_w(2),repo_tcmc(1,2),repo_tmcp(2),repo_tip(2),repo_tt(2)], ...
                %     [repo_w(3),repo_tcmc(1,3),repo_tmcp(3),repo_tip(3),repo_tt(3)], ...
                %     '-o', 'MarkerSize',3,'MarkerFaceColor',	[0.7, 0.7, 0.7], 'Color',[0.7, 0.7, 0.7]);
                % IndexPlotre = plot3([repo_w(1),repo_imcp(1),repo_ipip(1),repo_idip(1),repo_it(1)], ...
                %     [repo_w(2),repo_imcp(2),repo_ipip(2),repo_idip(2),repo_it(2)], ...
                %     [repo_w(3),repo_imcp(3),repo_ipip(3),repo_idip(3),repo_it(3)], ...
                %     '-o', 'MarkerSize',3,'MarkerFaceColor',	[0.7, 0.7, 0.7], 'Color',[0.7, 0.7, 0.7]);
                % MiddlePlotre = plot3([repo_w(1),repo_mmcp(1),repo_mpip(1),repo_mdip(1),repo_mt(1)], ...
                %     [repo_w(2),repo_mmcp(2),repo_mpip(2),repo_mdip(2),repo_mt(2)], ...
                %     [repo_w(3),repo_mmcp(3),repo_mpip(3),repo_mdip(3),repo_mt(3)], ...
                %     '-o', 'MarkerSize',3,'MarkerFaceColor',	[0.7, 0.7, 0.7], 'Color',[0.7, 0.7, 0.7]);
                % RingPlotre = plot3([repo_w(1),repo_rmcp(1),repo_rpip(1),repo_rdip(1),repo_rt(1)], ...
                %     [repo_w(2),repo_rmcp(2),repo_rpip(2),repo_rdip(2),repo_rt(2)], ...
                %     [repo_w(3),repo_rmcp(3),repo_rpip(3),repo_rdip(3),repo_rt(3)], ...
                %     '-o', 'MarkerSize',3,'MarkerFaceColor',	[0.7, 0.7, 0.7], 'Color',[0.7, 0.7, 0.7]);
                % LittlePlotre = plot3([repo_w(1),repo_lmcp(1),repo_lpip(1),repo_ldip(1),repo_lt(1)], ...
                %     [repo_w(2),repo_lmcp(2),repo_lpip(2),repo_ldip(2),repo_lt(2)], ...
                %     [repo_w(3),repo_lmcp(3),repo_lpip(3),repo_ldip(3),repo_lt(3)], ...
                %     '-o', 'MarkerSize',3,'MarkerFaceColor',	[0.7, 0.7, 0.7], 'Color',[0.7, 0.7, 0.7]);
                % RightArmPlotre = plot3([repo_n(1),0,repo_s(1),repo_e(1),repo_w(1)], ...
                %         [repo_n(2),0,repo_s(2),repo_e(2),repo_w(2)], ...
                %         [repo_n(3),0,repo_s(3),repo_e(3),repo_w(3)], ...
                %         '-o', 'MarkerSize',3,'MarkerFaceColor',[0.7, 0.7, 0.7], 'Color',[0.7, 0.7, 0.7]);

                ThumbPlotre = plot3([0,knuckles(1,1),TMCP(1),TIP(1),TT(1)], ...
                    [0,knuckles(1,2),TMCP(2),TIP(2),TT(2)], ...
                    [0,knuckles(1,3),TMCP(3),TIP(3),TT(3)], ...
                    '-o', 'MarkerSize',3,'MarkerFaceColor',	[0.7, 0.7, 0.7], 'Color',[0.7, 0.7, 0.7]);
                IndexPlotre = plot3([0,knuckles(2,1),PIP(1,1),DIP(1,1),Tip(1,1)], ...
                    [0,knuckles(2,2),PIP(1,2),DIP(1,2),Tip(1,2)], ...
                    [0,knuckles(2,3),PIP(1,3),DIP(1,3),Tip(1,3)], ...
                    '-o', 'MarkerSize',3,'MarkerFaceColor',	[0.7, 0.7, 0.7], 'Color',[0.7, 0.7, 0.7]);
                MiddlePlotre = plot3([0,knuckles(3,1),PIP(2,1),DIP(2,1),Tip(2,1)], ...
                    [0,knuckles(3,2),PIP(2,2),DIP(2,2),Tip(2,2)], ...
                    [0,knuckles(3,3),PIP(2,3),DIP(2,3),Tip(2,3)], ...
                    '-o', 'MarkerSize',3,'MarkerFaceColor', [0.7, 0.7, 0.7], 'Color',[0.7, 0.7, 0.7]);
                RingPlotre = plot3([0,knuckles(4,1),PIP(3,1),DIP(3,1),Tip(3,1)], ...
                    [0,knuckles(4,2),PIP(3,2),DIP(3,2),Tip(3,2)], ...
                    [0,knuckles(4,3),PIP(3,3),DIP(3,3),Tip(3,3)], ...
                    '-o', 'MarkerSize',3,'MarkerFaceColor',	[0.7, 0.7, 0.7], 'Color',[0.7, 0.7, 0.7]);
                LittlePlotre = plot3([0,knuckles(5,1),PIP(4,1),DIP(4,1),Tip(4,1)], ...
                    [0,knuckles(5,2),PIP(4,2),DIP(4,2),Tip(4,2)], ...
                    [0,knuckles(5,3),PIP(4,3),DIP(4,3),Tip(4,3)], ...
                    '-o', 'MarkerSize',3,'MarkerFaceColor',	[0.7, 0.7, 0.7], 'Color',[0.7, 0.7, 0.7]);
    
                xlabel('x (mm)')
                ylabel('y (mm)')
                zlabel('z (mm)')
                titlestr = ["Mean Hand Shape","Hand Shape with Synergy","Minimum"];
                title(titlestr(a))
                view([190 30])
            else
                % ThumbPlotre = plot3([repo_w(1),repo_tcmc(1),repo_tmcp(1),repo_tip(1),repo_tt(1)], ...
                %     [repo_w(2),repo_tcmc(1,2),repo_tmcp(2),repo_tip(2),repo_tt(2)], ...
                %     [repo_w(3),repo_tcmc(1,3),repo_tmcp(3),repo_tip(3),repo_tt(3)], ...
                %     '-o', 'MarkerSize',3,'MarkerFaceColor',	'k', 'Color','k');
                % IndexPlotre = plot3([repo_w(1),repo_imcp(1),repo_ipip(1),repo_idip(1),repo_it(1)], ...
                %     [repo_w(2),repo_imcp(2),repo_ipip(2),repo_idip(2),repo_it(2)], ...
                %     [repo_w(3),repo_imcp(3),repo_ipip(3),repo_idip(3),repo_it(3)], ...
                %     '-o', 'MarkerSize',3,'MarkerFaceColor',	'k', 'Color','k');
                % MiddlePlotre = plot3([repo_w(1),repo_mmcp(1),repo_mpip(1),repo_mdip(1),repo_mt(1)], ...
                %     [repo_w(2),repo_mmcp(2),repo_mpip(2),repo_mdip(2),repo_mt(2)], ...
                %     [repo_w(3),repo_mmcp(3),repo_mpip(3),repo_mdip(3),repo_mt(3)], ...
                %     '-o', 'MarkerSize',3,'MarkerFaceColor',	'k', 'Color','k');
                % RingPlotre = plot3([repo_w(1),repo_rmcp(1),repo_rpip(1),repo_rdip(1),repo_rt(1)], ...
                %     [repo_w(2),repo_rmcp(2),repo_rpip(2),repo_rdip(2),repo_rt(2)], ...
                %     [repo_w(3),repo_rmcp(3),repo_rpip(3),repo_rdip(3),repo_rt(3)], ...
                %     '-o', 'MarkerSize',3,'MarkerFaceColor',	'k', 'Color','k');
                % LittlePlotre = plot3([repo_w(1),repo_lmcp(1),repo_lpip(1),repo_ldip(1),repo_lt(1)], ...
                %     [repo_w(2),repo_lmcp(2),repo_lpip(2),repo_ldip(2),repo_lt(2)], ...
                %     [repo_w(3),repo_lmcp(3),repo_lpip(3),repo_ldip(3),repo_lt(3)], ...
                %     '-o', 'MarkerSize',3,'MarkerFaceColor',	'k', 'Color','k');
                % RightArmPlotre = plot3([repo_n(1),0,repo_s(1),repo_e(1),repo_w(1)], ...
                %         [repo_n(2),0,repo_s(2),repo_e(2),repo_w(2)], ...
                %         [repo_n(3),0,repo_s(3),repo_e(3),repo_w(3)], ...
                %         '-o', 'MarkerSize',10,'MarkerFaceColor','k', 'Color','k');

                ThumbPlotre = plot3([0,knuckles(1,1),TMCP(1),TIP(1),TT(1)], ...
                    [0,knuckles(1,2),TMCP(2),TIP(2),TT(2)], ...
                    [0,knuckles(1,3),TMCP(3),TIP(3),TT(3)], ...
                    '-o', 'MarkerSize',3,'MarkerFaceColor',	'k', 'Color','k');
                IndexPlotre = plot3([0,knuckles(2,1),PIP(1,1),DIP(1,1),Tip(1,1)], ...
                    [0,knuckles(2,2),PIP(1,2),DIP(1,2),Tip(1,2)], ...
                    [0,knuckles(2,3),PIP(1,3),DIP(1,3),Tip(1,3)], ...
                    '-o', 'MarkerSize',3,'MarkerFaceColor',	'k', 'Color','k');
                MiddlePlotre = plot3([0,knuckles(3,1),PIP(2,1),DIP(2,1),Tip(2,1)], ...
                    [0,knuckles(3,2),PIP(2,2),DIP(2,2),Tip(2,2)], ...
                    [0,knuckles(3,3),PIP(2,3),DIP(2,3),Tip(2,3)], ...
                    '-o', 'MarkerSize',3,'MarkerFaceColor',	'k', 'Color','k');
                RingPlotre = plot3([0,knuckles(4,1),PIP(3,1),DIP(3,1),Tip(3,1)], ...
                    [0,knuckles(4,2),PIP(3,2),DIP(3,2),Tip(3,2)], ...
                    [0,knuckles(4,3),PIP(3,3),DIP(3,3),Tip(3,3)], ...
                    '-o', 'MarkerSize',3,'MarkerFaceColor',	'k', 'Color','k');
                LittlePlotre = plot3([0,knuckles(5,1),PIP(4,1),DIP(4,1),Tip(4,1)], ...
                    [0,knuckles(5,2),PIP(4,2),DIP(4,2),Tip(4,2)], ...
                    [0,knuckles(5,3),PIP(4,3),DIP(4,3),Tip(4,3)], ...
                    '-o', 'MarkerSize',3,'MarkerFaceColor',	'k', 'Color','k');
            end
        end
    end
end