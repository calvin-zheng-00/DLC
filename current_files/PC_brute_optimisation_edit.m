% This file reprojects the PCs back into 3d position data and compares the
% performance of the principal components to the original 3d data at the
% moment of grasp.

close all; clear all; clc;

%src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\data\high_thresh\angles';
% src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\categorized_angles_2';
src_3d = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\categorized_2';
pc_src = 'C:\Users\czhe0008\Documents\MATLAB\Openpose\PCA_coeffs\';
mean_src = 'C:\Users\czhe0008\Documents\MATLAB\Openpose\PCA_mean\';
pro_src = 'C:\Users\czhe0008\Documents\MATLAB\Openpose\Projection\';
time_src = 'C:\Users\czhe0008\Documents\MATLAB\Openpose\PCA_times\';
unpro_src = 'C:\Users\czhe0008\Documents\MATLAB\Openpose\Unprojected\';
sig_src = 'C:\Users\czhe0008\Documents\MATLAB\Openpose\Significance\';
comp_src = 'C:\Users\czhe0008\Documents\MATLAB\Openpose\Comparison\';
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

allpcarr = load("All\pconcat_coeff.mat");
allpcarr = allpcarr.coeff;
allproarr = load("All\pconcat_project.mat");
allproarr = allproarr.score;
allmeanarr = load("All\pconcat_mean.mat");
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
sig_all = matfile('All\pAll_sig_80.mat');
sig_all = sig_all.significant_80;
maintemp = dir(fullfile(src_3d,'*'));
mainfolder = setdiff({maintemp([maintemp.isdir]).name},{'.','..'});

% cat1
% 1	NaN	NaN	NaN
% 1	4	NaN	NaN
% 1	2	9	NaN
% 1	2	4	9
%1 2 9 4 5
% 13.3244

% cat2
% 11	NaN	NaN	NaN	NaN
% 10	22	NaN	NaN	NaN
% 1	5	11	NaN	NaN
% 1	5	10	22	NaN
% 1	5	11	15	18

% cat2 alt
% 1 4 5 13 18
% Error: 14.69

% cat3
% 1	NaN	NaN	NaN	NaN
% 1	5	NaN	NaN	NaN
% 1	3	6	NaN	NaN
% 1	3	5	7	NaN
% 1	3	5	6	17

%% Methodology
% Start by brute forcing 2PCs and 3PCs.
% For every PC, Calculate previous optimum +1PC, opt-1 +2PCs, and opt-2 +3PCs.
% Select the best performance

focus_cat = 1;
allavgerr = [];
pcs = 1:28;
prechose = [2,9,15,25];
opt_num = 1;
pcs(prechose) = [];
opt = nchoosek(pcs,opt_num);
loopnum = height(opt);
avgerrcomp = [];
for r = 1:loopnum
    fprintf(1, 'Loops left %d\n', loopnum-r);
    repropc = [prechose,opt(r,:)];
    reproarr = allproarr(:,repropc)*allpcarr(:,repropc)' + repmat(allmeanarr,height(allproarr),1);
    avgerrcomp2 = [];
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
        file = fullfile(sig_src,sigfolder{i});
        sig = matfile(file);
        sig = sig.curr_sig_80;                      % Significant PCs for each participant
    
        %% Reprojection
        num_pcs = width(unproarr);
    
        %% Edit which PCs are being reprojected
        % reproarr = proarr(:,1:sig)*pcarr(:,1:sig)' + repmat(meanarr,height(proarr),1);
        % reproarr = proarr(:,:)*pcarr(:,:)' + repmat(meanarr,height(proarr),1);
        % reproarr(any(isnan(reproarr), 2), :) = [];
        subtemp = dir(fullfile(src_3d,mainfolder{i},'*.csv'));
        subfolder = {subtemp(~[subtemp.isdir]).name};
        total_error = [];
        if i < focus_cat
            reproarr(1:timeTable.end_time(end),:) = [];
            continue
        end
        if i > focus_cat
            break
        end
        rand_order = 1:numel(subfolder);
        rand_order = rand_order(randperm(length(rand_order)));
        rand_order = rand_order(1:min(100,length(rand_order)));
        for j_temp = 1:100
            j = rand_order(j_temp);
        % for j = 1:numel(subfolder)
            %% Getting limb lengths
            grasp_time_3d = timeTable.grasp_time(j);
            grasp_time = timeTable.reach_time(j);
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
    
            %% Calculating error
            % ComparisonArray = [knuckles(1,:),TMCP,TIP,TT,knuckles(2,:),PIP(1,:),DIP(1,:),Tip(1,:),...
            %     knuckles(3,:),PIP(2,:),DIP(2,:),Tip(2,:),knuckles(4,:),PIP(3,:),DIP(3,:),Tip(3,:),...
            %     knuckles(5,:),PIP(4,:),DIP(4,:),Tip(4,:),E,S,C,N];
            ComparisonArray = [repo_tcmc,repo_tmcp,repo_tip,repo_tt,repo_imcp,repo_ipip,repo_idip,repo_it,...
                    repo_mmcp,repo_mpip,repo_mdip,repo_mt,repo_rmcp,repo_rpip,repo_rdip,repo_rt,...
                    repo_lmcp,repo_lpip,repo_ldip,repo_lt,repo_e,repo_s,0,0,0,repo_n,repo_w];
            Error = table2array(tempcompTable)-ComparisonArray;
            Distance = []*(length(Error)/3);
            for k = 1:length(Error)/3
                Distance(k) = sqrt(Error((k-1)*3+1)^2+Error((k-1)*3+2)^2+Error((k-1)*3+3)^2);
            end
            total_error = [total_error;Distance];
        end
        % figure
        % boxchart(total_error(:,1:20))
        % xlabel('Category')
        % ylabel('Error (mm)')
        % avg_error = mean(median(total_error(:,1:20),"omitnan"),"omitnan");
        % titlestr = strcat("Average error: ", string(avg_error));
        % title(titlestr, 'Interpreter', 'none')
        % xticklabels(colnames3d)
        % total_error = array2table(total_error,"VariableNames",colnames3d);
        % writetable(total_error, strcat('Error_test\',mainfolder{i},'_error.csv'));
        %% optimisation criteria
        % if ismember(focus_cat, [1,6,7,8,9,12])
        %     avg_error = mean(median(total_error(:,[4,8,12,16,20]),"omitnan"),"omitnan"); % Finger tips
        % elseif focus_cat == 2
        %     avg_error = mean(median(total_error(:,[4,6,7,8]),"omitnan"),"omitnan"); % Thumb tip and index
        % elseif ismember(focus_cat, [4,11,13])
        %     avg_error = mean(median(total_error(:,[4,8,12]),"omitnan"),"omitnan"); % Thumb, index and middle tips
        % elseif focus_cat == 5
        %     avg_error = mean(median(total_error(:,[4,8,12,16]),"omitnan"),"omitnan"); % Thumb, index, middle and ring tips
        % elseif focus_cat == 10
        %     avg_error = mean(median(total_error(:,[4,8]),"omitnan"),"omitnan"); % Thumb and index, tips
        % elseif focus_cat == 3
        %     avg_error = mean(median(total_error(:,23),"omitnan"),"omitnan"); % Chest location
        % end
        if focus_cat == 1
            % avg_error = mean(median(total_error(:,[4,8,12,16,20]),"omitnan"),"omitnan"); % Finger tips
                avg_error = mean(median(total_error(:,21),"omitnan"),"omitnan"); % chest
        elseif focus_cat == 2
            avg_error = mean(median(total_error(:,[4,8]),"omitnan"),"omitnan"); % Thumb tip and index
        elseif focus_cat == 3
            avg_error = mean(median(total_error(:,[4,8,12]),"omitnan"),"omitnan"); % Thumb, index and middle tips
        end
        avgerrcomp2 = [avgerrcomp2;avg_error];
    end
    avgerrcomp = [avgerrcomp;mean(avgerrcomp2)];
end
[minimum,idx] = min(avgerrcomp);
good_pc = [prechose, opt(idx,:)];
allavgerr = [allavgerr,minimum];


figure
b = bar(allavgerr);
title('Average error of finger tips during grasp')
xlabel('PC number')
ylabel('Error (mm)')
box off