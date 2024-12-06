close all; clear all; clc;
% src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\categorized';
% dest = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\categorized_angles_2';
src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\filter_participant';
dest = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\angles_participant_2';

% cattemp = dir(fullfile(src,'*'));
% catfolder = setdiff({cattemp([cattemp.isdir]).name},{'.','..'});
% for catfolder_i = 1:numel(catfolder)
%     trialtemp = dir(fullfile(src,catfolder{catfolder_i},'*.csv'));
%     trialfolder = {trialtemp(~[trialtemp.isdir]).name};
%     for trialfolder_j = 1:numel(trialfolder)
%         data = fullfile(src,catfolder{catfolder_i},trialfolder{trialfolder_j});
%         outfolder = fullfile(dest,catfolder{catfolder_i});
%         outloc = fullfile(outfolder,strcat(trialfolder{trialfolder_j}(1:end-13),"_angles.csv"));

participanttemp = dir(fullfile(src,'*'));
participantfolder = setdiff({participanttemp([participanttemp.isdir]).name},{'.','..'});
for participantfolder_i = 1:numel(participantfolder)
    activitytemp = dir(fullfile(src,participantfolder{participantfolder_i},'*'));
    activityfolder = setdiff({activitytemp([activitytemp.isdir]).name},{'.','..'});  
    for activityfolder_i = 1:numel(activityfolder)
        trialtemp = dir(fullfile(src,participantfolder{participantfolder_i},activityfolder{activityfolder_i},'*.csv'));
        trialfolder = {trialtemp(~[trialtemp.isdir]).name};
        for trialfolder_j = 1:numel(trialfolder)
            data = fullfile(src,participantfolder{participantfolder_i},activityfolder{activityfolder_i},trialfolder{trialfolder_j});
            outfolder = fullfile(dest,participantfolder{participantfolder_i},activityfolder{activityfolder_i});
            outloc = fullfile(outfolder,strcat(trialfolder{trialfolder_j}(1:end-13),"_angles.csv"));

            if not(isfolder(outfolder))
                mkdir(outfolder)
            end
    
            %% Elbow
            finger = ["TCMC_x","TMCP_x","TIP_x","IMCP_x","IPIP_x","IDIP_x","MMCP_x","MPIP_x","MDIP_x","RMCP_x","RPIP_x","RDIP_x","LMCP_x","LPIP_x","LDIP_x"];
            abduct = ["IMCP_x","MMCP_x","RMCP_x","LMCP_x"];
            df = readtable(data);
            idx = find(strcmp(df.Properties.VariableNames, "RE_x"), 1);
            v1 = [df.(idx+3) - df.(idx), df.(idx+4) - df.(idx+1), df.(idx+5) - df.(idx+2)];
            v2 = [df.(idx-3) - df.(idx), df.(idx-2) - df.(idx+1), df.(idx-1) - df.(idx+2)];
            v1mag = sqrt(v1(:,1).^2 + v1(:,2).^2 + v1(:,3).^2);
            v1norm = v1./v1mag;
            v2mag = sqrt(v2(:,1).^2 + v2(:,2).^2 + v2(:,3).^2);
            v2norm = v2./v2mag;
            res = v1norm(:,1) .* v2norm(:,1) + v1norm(:,2) .* v2norm(:,2) + v1norm(:,3) .* v2norm(:,3);
            df2 = table(acos(res)-pi);
            df2.Properties.VariableNames = "RE_flex";
            %% Flexion
            for joint_idx = 1:length(finger)
                joint = char(finger(joint_idx));
                col_name = strcat(joint(1:end-2),"_flex");
                idx = find(strcmp(df.Properties.VariableNames, joint), 1);
                if mod(joint_idx-1,3) == 0
                    v1 = [df.(idx+3) - df.(idx), df.(idx+4) - df.(idx+1), df.(idx+5) - df.(idx+2)];
                    v2 = [df.RW_x - df.(idx), df.RW_y - df.(idx+1), df.RW_z - df.(idx+2)];
                else
                    v1 = [df.(idx+3) - df.(idx), df.(idx+4) - df.(idx+1), df.(idx+5) - df.(idx+2)];
                    v2 = [df.(idx-3) - df.(idx), df.(idx-2) - df.(idx+1), df.(idx-1) - df.(idx+2)];
                end
                v1mag = sqrt(v1(:,1).^2 + v1(:,2).^2 + v1(:,3).^2);
                v1norm = v1./v1mag;
                v2mag = sqrt(v2(:,1).^2 + v2(:,2).^2 + v2(:,3).^2);
                v2norm = v2./v2mag;
                res = v1norm(:,1) .* v2norm(:,1) + v1norm(:,2) .* v2norm(:,2) + v1norm(:,3) .* v2norm(:,3);
                df2.(col_name) = acos(res)-pi;
            end
        
            %% Palm
            % Calculating the normal vector to the plane that is the palm (using the wrist, index knuckle and little knuckle)
            v1 = [df.IMCP_x - df.RW_x, df.IMCP_y - df.RW_y, df.IMCP_z - df.RW_z];
            v2 = [df.LMCP_x - df.RW_x, df.LMCP_y - df.RW_y, df.LMCP_z - df.RW_z];
            normt = cross(v1,v2,2); % Normal to the palm plane
            normmag = sqrt(normt(:,1).^2 + normt(:,2).^2 + normt(:,3).^2);
            unitnorm = normt./normmag;
        
            %% Abduction
            % Calculating angle of abduction (angle between MCP-PIP and MCP-Wrist instead of between two different fingers)
            % Changed abduction angle to be between v1_p and right side vector
            for joint_idx = 1:length(abduct)
                joint = char(abduct(joint_idx));
                col_name = strcat(joint(1:end-2),"_abd");
                col_name2 = strcat(joint(1:end-2),"_flex");
                idx = find(strcmp(df.Properties.VariableNames, joint), 1);
                v1 = [df.(idx+3) - df.(idx), df.(idx+4) - df.(idx+1), df.(idx+5) - df.(idx+2)];
                v2 = [df.RW_x - df.(idx), df.RW_y - df.(idx+1), df.RW_z - df.(idx+2)];
                % v1_np = normt.*((v1(:,1).*normt(:,1) + v1(:,2).*normt(:,2) + v1(:,3).*normt(:,3))./(normmag.^2));  % Vector component normal to the palm plane
                % v1_p = [v1(:,1)-v1_np(:,1),v1(:,2)-v1_np(:,2),v1(:,3)-v1_np(:,3)];                                       % Vector component on the palm plane
                % v1mag = sqrt(v1_p(:,1).^2 + v1_p(:,2).^2 + v1_p(:,3).^2);
                % v1norm = v1_p./v1mag;
                v2_np = normt.*((v2(:,1).*normt(:,1) + v2(:,2).*normt(:,2) + v2(:,3).*normt(:,3))./(normmag.^2));   % Vector component normal to the palm plane
                v2_p = [v2(:,1)-v2_np(:,1),v2(:,2)-v2_np(:,2),v2(:,3)-v2_np(:,3)];                                       % Vector component on the palm plane
                v2mag = sqrt(v2_p(:,1).^2 + v2_p(:,2).^2 + v2_p(:,3).^2);
                v2norm = v2_p./v2mag;
                %Calculate knuckle_axis as the cross product between the palm plane and the component of V2 along the palm plane
                knuckle_axis = cross(v2norm,normt);
                knuckle_axismag = sqrt(knuckle_axis(:,1).^2 + knuckle_axis(:,2).^2 + knuckle_axis(:,3).^2);
                knuckle_axisunit = knuckle_axis./knuckle_axismag;
                
                z = dot(v1,unitnorm,2); % magnitude of the component of v1 along the palm normal
                x = dot(v1,v2norm,2); % magnitude of the component of v1 along the CMC vector
                y = dot(v1,knuckle_axisunit,2); % magnitude of the component of v1 along the knuckle axis
                df2.(col_name) = atan(y./x);
                df2.(col_name2) = acos(z./sqrt(x.^2 + y.^2 + z.^2)) - pi/2;
                
                % res = v1norm(:,1) .* rightsideunit(:,1) + v1norm(:,2) .* rightsideunit(:,2) + v1norm(:,3) .* rightsideunit(:,3);
                % df2.(col_name) = (pi/2) - acos(res);          % Minus 90 degrees to keep angle centred around 0
                % % Fixing up knuckle flex angle
                % % Component of V1 along the rightside plane
                % v1_np = rightside.*((v1(:,1).*rightside(:,1) + v1(:,2).*rightside(:,2) + v1(:,3).*rightside(:,3))./(rightsidemag.^2));  % Vector component normal to the palm plane
                % v1_p = [v1(:,1)-v1_np(:,1),v1(:,2)-v1_np(:,2),v1(:,3)-v1_np(:,3)];                                    % Vector component on the palm plane
                % v1mag = sqrt(v1_p(:,1).^2 + v1_p(:,2).^2 + v1_p(:,3).^2);
                % v1norm = v1_p./v1mag;
                % res = v1norm(:,1) .* unitnorm(:,1) + v1norm(:,2) .* unitnorm(:,2) + v1norm(:,3) .* unitnorm(:,3);
                % df2.(col_name2) = (pi/2) - acos(res);          % Minus 90 degrees to keep angle centred around 0
            end
        
            %% Thumb
            v1 = [df.TMCP_x - df.TCMC_x, df.TMCP_y - df.TCMC_y, df.TMCP_z - df.TCMC_z];
            v2 = [df.IMCP_x - df.TCMC_x, df.IMCP_y - df.TCMC_y, df.IMCP_z - df.TCMC_z];
            v3 = [df.RW_x - df.TCMC_x, df.RW_y - df.TCMC_y, df.RW_z - df.TCMC_z];
            
            v2mag = sqrt(v2(:,1).^2 + v2(:,2).^2 + v2(:,3).^2);
            tmb_norm = cross(v2,v3,2);  % Normal to the thumb plane
            tmb_normmag = sqrt(tmb_norm(:,1).^2 + tmb_norm(:,2).^2 + tmb_norm(:,3).^2);
            tmb_unitnorm = tmb_norm./tmb_normmag;
            tabd_norm = cross(v2,tmb_norm,2);  % Orthogonal to the thumb plane and index vector
            tabd_normmag = sqrt(tabd_norm(:,1).^2 + tabd_norm(:,2).^2 + tabd_norm(:,3).^2);
            
            z = dot(v1,tmb_norm,2)./tmb_normmag; % magnitude of the component of v1 along the thumb normal
            x = dot(v1,v2,2)./v2mag; % magnitude of the component of v1 along the index vector
            y = dot(v1,tabd_norm,2)./tabd_normmag; % magnitude of the component of v1 along the index vector
            df2.("TCMC_flex") = -atan(y./x);
            df2.("TMCP_abd") = -acos(z./sqrt(x.^2 + y.^2 + z.^2)) + pi/2;

            %TCMC_flex is the angle between TCMC-IMCP and TCMC-TMCP along the thumb plane (RW-TCMC-IMCP)
            %TMCP_abd is the angle between TCMC-IMCP and TCMC-TMCP out of the thumb plane (along the tabd plane described by the cross product between the thumb pane normal and TCMC-IMCP)
            %TCMC_rot is the angle between the component of the thumb rotation plane (TCMC-TMCP-TIP) that's orthogonal to TCMC-IMCP, and the tabd plane
            
            % tmb_norm = cross(v2,v3,2);  % Normal to the thumb plane
            % tmb_normmag = sqrt(tmb_norm(:,1).^2 + tmb_norm(:,2).^2 + tmb_norm(:,3).^2);
            % tmb_unitnorm = tmb_norm./tmb_normmag;
            % v1_np = tmb_norm.*((v1(:,1).*tmb_norm(:,1) + v1(:,2).*tmb_norm(:,2) + v1(:,3).*tmb_norm(:,3))./(tmb_normmag.^2));  % Vector component normal to the thumb plane
            % v1_p = [v1(:,1)-v1_np(:,1),v1(:,2)-v1_np(:,2),v1(:,3)-v1_np(:,3)];                                    % Vector component on the thumb plane
            % v1mag = sqrt(v1_p(:,1).^2 + v1_p(:,2).^2 + v1_p(:,3).^2);
            % v1norm = v1_p./v1mag;
            % v2mag = sqrt(v2(:,1).^2 + v2(:,2).^2 + v2(:,3).^2);
            % v2norm = v2./v2mag;
            % res = v1norm(:,1) .* v2norm(:,1) + v1norm(:,2) .* v2norm(:,2) + v1norm(:,3) .* v2norm(:,3);
            % df2.("TCMC_flex") = acos(res);
            % % TMCP_abd is the angle between TCMC-IMCP and TCMC-TMCP out of the thumb plane (along the tabd plane described by the cross product between the thumb plane normal and TCMC-IMCP)
            % tabd_norm = cross(v2,tmb_norm,2);  % Normal to the thumb plane
            % tabd_normmag = sqrt(tabd_norm(:,1).^2 + tabd_norm(:,2).^2 + tabd_norm(:,3).^2);
            % v1_np = tabd_norm.*((v1(:,1).*tabd_norm(:,1) + v1(:,2).*tabd_norm(:,2) + v1(:,3).*tabd_norm(:,3))./(tabd_normmag.^2));  % Vector component normal to the thumb abduction plane
            % v1_p = [v1(:,1)-v1_np(:,1),v1(:,2)-v1_np(:,2),v1(:,3)-v1_np(:,3)];                                    % Vector component on the thumb abduction plane
            % v1mag = sqrt(v1_p(:,1).^2 + v1_p(:,2).^2 + v1_p(:,3).^2);
            % v1norm = v1_p./v1mag;
            % % res = v1norm(:,1) .* v2norm(:,1) + v1norm(:,2) .* v2norm(:,2) + v1norm(:,3) .* v2norm(:,3);
            % res = v1norm(:,1) .* tmb_unitnorm(:,1) + v1norm(:,2) .* tmb_unitnorm(:,2) + v1norm(:,3) .* tmb_unitnorm(:,3);
            % df2.("TMCP_abd") = acos(res) - pi/2;
            % TCMC_rot is the angle between the thumb rotation plane (TCMC-TMCP-TIP) and the component of the thumb plane that's orthogonal to TCMC-TMCP
            v1mag = sqrt(v1(:,1).^2 + v1(:,2).^2 + v1(:,3).^2);
            tmb_unitnorm_np = v1.*((tmb_unitnorm(:,1).*v1(:,1) + tmb_unitnorm(:,2).*v1(:,2) + tmb_unitnorm(:,3).*v1(:,3))./(v1mag.^2));  % thumb plane component along v1
            tmb_unitnorm_p = tmb_unitnorm - tmb_unitnorm_np;                               % thumb plane component orthogonal to v1
            tmb_unitnorm_p_mag = sqrt(tmb_unitnorm_p(:,1).^2 + tmb_unitnorm_p(:,2).^2 + tmb_unitnorm_p(:,3).^2);
            tmb_unitnorm_p_norm = tmb_unitnorm_p./tmb_unitnorm_p_mag;
            v4 = [df.TIP_x - df.TMCP_x, df.TIP_y - df.TMCP_y, df.TIP_z - df.TMCP_z];
            trot_norm = cross(v1,v4,2);  % Normal to the thumb plane
            trotmag = sqrt(trot_norm(:,1).^2 + trot_norm(:,2).^2 + trot_norm(:,3).^2);
            trot_norm = trot_norm./trotmag;
            res = trot_norm(:,1) .* tmb_unitnorm_p_norm(:,1) + trot_norm(:,2) .* tmb_unitnorm_p_norm(:,2) + trot_norm(:,3) .* tmb_unitnorm_p_norm(:,3);
            df2.("TCMC_rot") = acos(res);
        
            %% Wrist
            v1 = [df.RW_x - df.RE_x, df.RW_y - df.RE_y, df.RW_z - df.RE_z];
            v2 = [df.MMCP_x - df.RW_x, df.MMCP_y - df.RW_y, df.MMCP_z - df.RW_z];
            v2_np = normt.*((v2(:,1).*normt(:,1) + v2(:,2).*normt(:,2) + v2(:,3).*normt(:,3))./(normmag.^2));   % Vector component normal to the palm plane
            v2_p = [v2(:,1)-v2_np(:,1),v2(:,2)-v2_np(:,2),v2(:,3)-v2_np(:,3)];                                       % Vector component on the palm plane
            v2mag = sqrt(v2_p(:,1).^2 + v2_p(:,2).^2 + v2_p(:,3).^2);
            wrist_forward = v2_p./v2mag;
            wrist_fe_axis = cross(normt,wrist_forward,2);  % Wrist flexion axis
            wrist_fe_axis_mag = sqrt(wrist_fe_axis(:,1).^2 + wrist_fe_axis(:,2).^2 + wrist_fe_axis(:,3).^2);
            
            z = dot(v1,unitnorm,2); % magnitude of the component of v1 along the thumb normal
            x = dot(v1,wrist_fe_axis,2)./wrist_fe_axis_mag; % magnitude of the component of v1 along the index vector
            y = dot(v1,wrist_forward,2); % magnitude of the component of v1 along the index vector
            df2.("W_abd") = atan(y./x);
            df2.("W_flex") = -acos(z./sqrt(x.^2 + y.^2 + z.^2)) + pi/2;

            % Find a reference vector orthogonal to the forearm
            % To calculate rotation of wrist, find the plane generated by the wrist, elbow and shoulder (arm plane).
            v1mag = sqrt(v1(:,1).^2 + v1(:,2).^2 + v1(:,3).^2);
            forearm_np = v1.*((normt(:,1).*v1(:,1) + normt(:,2).*v1(:,2) + normt(:,3).*v1(:,3))./(v1mag.^2));   % Vector component normal to the palm plane
            forearm_p = [normt(:,1)-forearm_np(:,1),normt(:,2)-forearm_np(:,2),normt(:,3)-forearm_np(:,3)];                                       % Vector component on the palm plane
            forearmmag = sqrt(forearm_p(:,1).^2 + forearm_p(:,2).^2 + forearm_p(:,3).^2);
            forearm_ref = forearm_p./forearmmag;
            v2 = [df.RS_x - df.RE_x, df.RS_y - df.RE_y, df.RS_z - df.RE_z];
            arm_norm = cross(v1,v2,2);
            arm_normmag = sqrt(arm_norm(:,1).^2 + arm_norm(:,2).^2 + arm_norm(:,3).^2);
            unit_arm_norm = arm_norm./arm_normmag;
            res = forearm_ref(:,1) .* unit_arm_norm(:,1) + forearm_ref(:,2) .* unit_arm_norm(:,2) + forearm_ref(:,3) .* unit_arm_norm(:,3);
            df2.("W_rot") = acos(res); 

            % Find wrist flex vs abduction
            % We have the normal of the wrist plane and two points for the orthogonal plane to pass through. (MMCP and Wrist)
            % Find the cross product between the normal to the palm plane and the vector (MMCP-Wrist) to get the plane
            % along which we will be finding the flex angle. Subtract flex angle from total angle to get abduction angle.
            % vp = [df.MMCP_x - df.RW_x, df.MMCP_y - df.RW_y, df.MMCP_z - df.RW_z];
            % vp_np = normt.*((vp(:,1).*normt(:,1) + vp(:,2).*normt(:,2) + vp(:,3).*normt(:,3))./(normmag.^2));  % Vector component normal to the palm plane
            % vp_p = [vp(:,1)-vp_np(:,1),vp(:,2)-vp_np(:,2),vp(:,3)-vp_np(:,3)];                                    % Vector component on the palm plane
            % vpmag = sqrt(vp_p(:,1).^2 + vp_p(:,2).^2 + vp_p(:,3).^2);
            % vp = vp_p./vpmag;   % finding RW-MMCP along palm plane first
            % flex_plane = cross(normt,vp,2);
            % flex_plane_mag = sqrt(flex_plane(:,1).^2 + flex_plane(:,2).^2 + flex_plane(:,3).^2);
            % flex_plane_unit = flex_plane./flex_plane_mag;
            % v1 = [df.RE_x - df.RW_x, df.RE_y - df.RW_y, df.RE_z - df.RW_z];
            % v1_np = flex_plane.*((v1(:,1).*flex_plane(:,1) + v1(:,2).*flex_plane(:,2) + v1(:,3).*flex_plane(:,3))./(flex_plane_mag.^2));  % Vector component normal to the flex plane
            % v1_p = [v1(:,1)-v1_np(:,1),v1(:,2)-v1_np(:,2),v1(:,3)-v1_np(:,3)];                                    % Vector component on the flex plane
            % v1mag = sqrt(v1_p(:,1).^2 + v1_p(:,2).^2 + v1_p(:,3).^2);
            % v1norm = v1_p./v1mag;
            % res = v1norm(:,1) .* unitnorm(:,1) + v1norm(:,2) .* unitnorm(:,2) + v1norm(:,3) .* unitnorm(:,3);
            % df2.("W_flex") = (pi/2) - acos(res);
            % % Abuct = angle between RE-RW and negative RW-MMCP along the palm plane (need to find RW-MMCP along palm plane)
            % v1_np = normt.*((v1(:,1).*normt(:,1) + v1(:,2).*normt(:,2) + v1(:,3).*normt(:,3))./(normmag.^2));   % Vector component normal to the palm plane
            % v1_p = [v1(:,1)-v1_np(:,1),v1(:,2)-v1_np(:,2),v1(:,3)-v1_np(:,3)];                                    % Vector component on the palm plane
            % v1mag = sqrt(v1_p(:,1).^2 + v1_p(:,2).^2 + v1_p(:,3).^2);
            % v1norm = v1_p./v1mag;
            % res = v1norm(:,1) .* flex_plane_unit(:,1) + v1norm(:,2) .* flex_plane_unit(:,2) + v1norm(:,3) .* flex_plane_unit(:,3);
            % df2.("W_abd") = (pi/2) - acos(res);    
            % To calculate rotation of wrist, find the plane generated by the wrist, elbow and shoulder (arm plane).
            % Find the angle between the palm plane and the arm plane.
            % v1 = [df.RW_x - df.RE_x, df.RW_y - df.RE_y, df.RW_z - df.RE_z];
            % v2 = [df.RS_x - df.RE_x, df.RS_y - df.RE_y, df.RS_z - df.RE_z];
            % arm_norm = cross(v1,v2,2);
            % arm_normmag = sqrt(arm_norm(:,1).^2 + arm_norm(:,2).^2 + arm_norm(:,3).^2);
            % unit_arm_norm = arm_norm./arm_normmag;
            % % Find plane perpendicular to RE_RW vector; RE_RW is the normal to the plane
            % v1mag = sqrt(v1(:,1).^2 + v1(:,2).^2 + v1(:,3).^2); % Forearm
            % % Find angle between palm normal and arm normal
            % palm_np = v1.*((unitnorm(:,1).*v1(:,1) + unitnorm(:,2).*v1(:,2) + unitnorm(:,3).*v1(:,3))./(v1mag.^2));   % Vector component along the forearm
            % palm_p = [unitnorm(:,1)-palm_np(:,1),unitnorm(:,2)-palm_np(:,2),unitnorm(:,3)-palm_np(:,3)];                                    % The vector component of the palm normal that's orthogonal to the forearm
            % palmmag = sqrt(palm_p(:,1).^2 + palm_p(:,2).^2 + palm_p(:,3).^2);
            % palmnorm = palm_p./palmmag;
            % res = palmnorm(:,1) .* unit_arm_norm(:,1) + palmnorm(:,2) .* unit_arm_norm(:,2) + palmnorm(:,3) .* unit_arm_norm(:,3);
            % df2.("W_rot") = acos(res);    
        
            %% Shoulder
            v1 = [df.RS_x - df.RE_x, df.RS_y - df.RE_y, df.RS_z - df.RE_z];
            v1mag = sqrt(v1(:,1).^2 + v1(:,2).^2 + v1(:,3).^2);
            v1norm = v1./v1mag;
            forearm_axis = cross(unit_arm_norm,v1norm,2);
            forearm_axismag = sqrt(forearm_axis(:,1).^2 + forearm_axis(:,2).^2 + forearm_axis(:,3).^2);
            forearm_axis = forearm_axis./forearm_axismag;
            v2 = [df.RS_x - df.C_x, df.RS_y - df.C_y, df.RS_z - df.C_z];
            z = dot(v2,unit_arm_norm,2); % magnitude of the component of v1 along the thumb normal
            x = dot(v2,v1norm,2); % magnitude of the component of v1 along the index vector
            y = dot(v2,forearm_axis,2); % magnitude of the component of v1 along the index vector
            df2.("RS_flex") = atan(y./x);
            df2.("RS_abd") = acos(z./sqrt(x.^2 + y.^2 + z.^2)) - pi/2;


            % Find shoulder abduction by using Shoulder, Chest, Nose to find the plane
            % v1 = [df.N_x - df.C_x, df.N_y - df.C_y, df.N_z - df.C_z];
            % v2 = [df.RS_x - df.C_x, df.RS_y - df.C_y, df.RS_z - df.C_z];
            % chest_norm = cross(v1,v2,2);
            % chest_normmag = sqrt(chest_norm(:,1).^2 + chest_norm(:,2).^2 + chest_norm(:,3).^2);
            % unit_chest_norm = chest_norm./chest_normmag;
            % v1 = [df.RE_x - df.RS_x, df.RE_y - df.RS_y, df.RE_z - df.RS_z];
            % v1_np = chest_norm.*((v1(:,1).*chest_norm(:,1) + v1(:,2).*chest_norm(:,2) + v1(:,3).*chest_norm(:,3))./(chest_normmag.^2)); % Vector component normal to the chest plane
            % v1_p = [v1(:,1)-v1_np(:,1),v1(:,2)-v1_np(:,2),v1(:,3)-v1_np(:,3)];                                    % Vector component on the chest plane
            % v1mag = sqrt(v1_p(:,1).^2 + v1_p(:,2).^2 + v1_p(:,3).^2);
            % v1norm = v1_p./v1mag;
            % v2 = [df.C_x - df.RS_x, df.C_y - df.RS_y, df.C_z - df.RS_z];
            % v2mag = sqrt(v2(:,1).^2 + v2(:,2).^2 + v2(:,3).^2);
            % v2norm = v2./v2mag;
            % res = v1norm(:,1) .* v2norm(:,1) + v1norm(:,2) .* v2norm(:,2) + v1norm(:,3) .* v2norm(:,3);
            % df2.("RS_abd") = acos(res); 
            % % Shoulder flexion (Forwards and backwards)
            % v1_np = v2.*((v1(:,1).*v2(:,1) + v1(:,2).*v2(:,2) + v1(:,3).*v2(:,3))./(v2mag.^2)); % Vector component along to C-RS
            % v1_p = [v1(:,1)-v1_np(:,1),v1(:,2)-v1_np(:,2),v1(:,3)-v1_np(:,3)];                                    % Vector component orthogonal to C-RS
            % v1mag = sqrt(v1_p(:,1).^2 + v1_p(:,2).^2 + v1_p(:,3).^2);
            % v1norm = v1_p./v1mag;
            % res = v1norm(:,1) .* unit_chest_norm(:,1) + v1norm(:,2) .* unit_chest_norm(:,2) + v1norm(:,3) .* unit_chest_norm(:,3);
            % df2.("RS_flex") = acos(res);

            % Shoulder rotation comparing arm plane (RW,RE,RS) to shoulder plane (C-RS-RE)
            % Finding shoulder plane (C-RS-RE)
            % v1 = [df.C_x - df.RS_x, df.C_y - df.RS_y, df.C_z - df.RS_z];
            % v2 = [df.RE_x - df.RS_x, df.RE_y - df.RS_y, df.RE_z - df.RS_z];
            % shoulder_norm = cross(v1,v2,2);
            % shoulder_normmag = sqrt(shoulder_norm(:,1).^2 + shoulder_norm(:,2).^2 + shoulder_norm(:,3).^2);
            % unit_shoulder_norm = shoulder_norm./shoulder_normmag;
            % res = unit_shoulder_norm(:,1) .* unit_arm_norm(:,1) + unit_shoulder_norm(:,2) .* unit_arm_norm(:,2) + unit_shoulder_norm(:,3) .* unit_arm_norm(:,3);
            % df2.("RS_rot") = acos(res);
            v1 = [df.C_x - df.RS_x, df.C_y - df.RS_y, df.C_z - df.RS_z];
            v2 = [df.RE_x - df.RS_x, df.RE_y - df.RS_y, df.RE_z - df.RS_z];
            v3 = [df.N_x - df.C_x, df.N_y - df.C_y, df.N_z - df.C_z];
            shoulder_norm = cross(v1,v2,2);
            shoulder_normmag = sqrt(shoulder_norm(:,1).^2 + shoulder_norm(:,2).^2 + shoulder_norm(:,3).^2);
            unit_shoulder_norm = shoulder_norm./shoulder_normmag;
            v1mag = sqrt(v1(:,1).^2 + v1(:,2).^2 + v1(:,3).^2);
            v3_np = v1.*((v3(:,1).*v1(:,1) + v3(:,2).*v1(:,2) + v3(:,3).*v1(:,3))./(v1mag.^2));   % Vector component normal to the palm plane
            v3_p = [v3(:,1)-v3_np(:,1),v3(:,2)-v3_np(:,2),v3(:,3)-v3_np(:,3)];                                       % Vector component on the palm plane
            v3mag = sqrt(v3_p(:,1).^2 + v3_p(:,2).^2 + v3_p(:,3).^2);
            v3unit = v3_p./v3mag;
            res = unit_shoulder_norm(:,1) .* v3unit(:,1) + unit_shoulder_norm(:,2) .* v3unit(:,2) + unit_shoulder_norm(:,3) .* v3unit(:,3);
            df2.("RS_rot") = acos(res);
    
            df2 = [df2(:,2:16) df2(:,21:22) df2(:,17:20) df2(:,23:25) df2(:,1) df2(:,26:28)];
        
            writetable(df2,outloc)
        end
    end
end