close all; clear all; clc;
%src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\angles_participant_2\p10\Instructions_BAT1\20240730T101748-101755_angles.csv';
%src_3d = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\filter_participant\p10\Instructions_BAT1\20240730T101748-101755_filtered.csv';
src_3d = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\filter_participant\p21\Instructions_DIG1\20240816T134146-134155_filtered.csv';
src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\angles_participant_2\p21\Instructions_DIG1\20240816T134146-134155_angles.csv';

unproTable = fullfile(src);
unproTable = readtable(unproTable);
%% Getting limb lengths
T = fullfile(src_3d);
T = readtable(T);
original = table2array(T);
original_table = T(:,[4:63,67:78]);
original = original(:,[4:63,67:78]);            
grasp_time = 164; %164
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

% Fix TCMC_flex, IMCP_flex, MMCP_flex, RMCP_flex, LMCP_flex,
% IMCP_abd, MMCP_abd, RMCP_abd, LMCP_abd, all W, all RS
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

%% Old knuckle code
% K = [0 -palm_plane(3) palm_plane(2); palm_plane(3) 0 -palm_plane(1); -palm_plane(2) palm_plane(1) 0];
% R = eye(3) + sin(-unproTable.IMCP_abd(grasp_time))*K + (1-cos(-unproTable.IMCP_abd(grasp_time)))*K*K;
% IPIP = R*knuckles(2,:)'.*(lengths(6)/lengths(5)) + knuckles(2,:)';
% new_plane = cross(knuckles(2,:),palm_plane);
% new_plane = new_plane/norm(new_plane);
% K = [0 -new_plane(3) new_plane(2); new_plane(3) 0 -new_plane(1); -new_plane(2) new_plane(1) 0];
% R = eye(3) + sin(-unproTable.IMCP_flex(grasp_time))*K + (1-cos(-unproTable.IMCP_flex(grasp_time)))*K*K;
% IPIP = R*(IPIP-knuckles(2,:)') + knuckles(2,:)';
% K = [0 -new_plane(3) new_plane(2); new_plane(3) 0 -new_plane(1); -new_plane(2) new_plane(1) 0];
% R = eye(3) + sin(unproTable.IPIP_flex(grasp_time)-pi)*K + (1-cos(unproTable.IPIP_flex(grasp_time)-pi))*K*K;
% IDIP = R*(IPIP-knuckles(2,:)').*(lengths(7)/lengths(6)) + IPIP;
% R = eye(3) + sin(unproTable.IDIP_flex(grasp_time)-pi)*K + (1-cos(unproTable.IDIP_flex(grasp_time)-pi))*K*K;
% IT = R*(IDIP-IPIP).*(lengths(8)/lengths(7)) + IDIP;
%% End of old code

new_plane = cross(knuckles(2,:),palm_plane);
new_plane = new_plane/norm(new_plane);
K = [0 -new_plane(3) new_plane(2); new_plane(3) 0 -new_plane(1); -new_plane(2) new_plane(1) 0];
R = eye(3) + sin(-unproTable.IMCP_flex(grasp_time) - pi/2)*K + (1-cos(-unproTable.IMCP_flex(grasp_time) - pi/2))*K*K;
IPIP = R*palm_plane'.*lengths(6);
K = [0 -palm_plane(3) palm_plane(2); palm_plane(3) 0 -palm_plane(1); -palm_plane(2) palm_plane(1) 0];
R = eye(3) + sin(-unproTable.IMCP_abd(grasp_time))*K + (1-cos(-unproTable.IMCP_abd(grasp_time)))*K*K;
IPIP = R*IPIP + knuckles(2,:)';
K = [0 -new_plane(3) new_plane(2); new_plane(3) 0 -new_plane(1); -new_plane(2) new_plane(1) 0];
R = eye(3) + sin(unproTable.IPIP_flex(grasp_time)-pi)*K + (1-cos(unproTable.IPIP_flex(grasp_time)-pi))*K*K;
IDIP = R*(IPIP-knuckles(2,:)').*(lengths(7)/lengths(6)) + IPIP;
R = eye(3) + sin(unproTable.IDIP_flex(grasp_time)-pi)*K + (1-cos(unproTable.IDIP_flex(grasp_time)-pi))*K*K;
IT = R*(IDIP-IPIP).*(lengths(8)/lengths(7)) + IDIP;

new_plane = cross(knuckles(3,:),palm_plane);
new_plane = new_plane/norm(new_plane);
K = [0 -new_plane(3) new_plane(2); new_plane(3) 0 -new_plane(1); -new_plane(2) new_plane(1) 0];
R = eye(3) + sin(-unproTable.MMCP_flex(grasp_time) - pi/2)*K + (1-cos(-unproTable.MMCP_flex(grasp_time) - pi/2))*K*K;
MPIP = R*palm_plane'.*lengths(10);
K = [0 -palm_plane(3) palm_plane(2); palm_plane(3) 0 -palm_plane(1); -palm_plane(2) palm_plane(1) 0];
R = eye(3) + sin(-unproTable.MMCP_abd(grasp_time))*K + (1-cos(-unproTable.MMCP_abd(grasp_time)))*K*K;
MPIP = R*MPIP + knuckles(3,:)';
K = [0 -new_plane(3) new_plane(2); new_plane(3) 0 -new_plane(1); -new_plane(2) new_plane(1) 0];
R = eye(3) + sin(unproTable.MPIP_flex(grasp_time)-pi)*K + (1-cos(unproTable.MPIP_flex(grasp_time)-pi))*K*K;
MDIP = R*(MPIP-knuckles(3,:)').*(lengths(11)/lengths(10)) + MPIP;
R = eye(3) + sin(unproTable.MDIP_flex(grasp_time)-pi)*K + (1-cos(unproTable.MDIP_flex(grasp_time)-pi))*K*K;
MT = R*(MDIP-MPIP).*(lengths(12)/lengths(11)) + MDIP;

new_plane = cross(knuckles(4,:),palm_plane);
new_plane = new_plane/norm(new_plane);
K = [0 -new_plane(3) new_plane(2); new_plane(3) 0 -new_plane(1); -new_plane(2) new_plane(1) 0];
R = eye(3) + sin(-unproTable.RMCP_flex(grasp_time) - pi/2)*K + (1-cos(-unproTable.RMCP_flex(grasp_time) - pi/2))*K*K;
RPIP = R*palm_plane'.*lengths(14);
K = [0 -palm_plane(3) palm_plane(2); palm_plane(3) 0 -palm_plane(1); -palm_plane(2) palm_plane(1) 0];
R = eye(3) + sin(-unproTable.RMCP_abd(grasp_time))*K + (1-cos(-unproTable.RMCP_abd(grasp_time)))*K*K;
RPIP = R*RPIP + knuckles(4,:)';
K = [0 -new_plane(3) new_plane(2); new_plane(3) 0 -new_plane(1); -new_plane(2) new_plane(1) 0];
R = eye(3) + sin(unproTable.RPIP_flex(grasp_time)-pi)*K + (1-cos(unproTable.RPIP_flex(grasp_time)-pi))*K*K;
RDIP = R*(RPIP-knuckles(4,:)').*(lengths(15)/lengths(14)) + RPIP;
R = eye(3) + sin(unproTable.RDIP_flex(grasp_time)-pi)*K + (1-cos(unproTable.RDIP_flex(grasp_time)-pi))*K*K;
RT = R*(RDIP-RPIP).*(lengths(16)/lengths(15)) + RDIP;

new_plane = cross(knuckles(5,:),palm_plane);
new_plane = new_plane/norm(new_plane);
K = [0 -new_plane(3) new_plane(2); new_plane(3) 0 -new_plane(1); -new_plane(2) new_plane(1) 0];
R = eye(3) + sin(-unproTable.LMCP_flex(grasp_time) - pi/2)*K + (1-cos(-unproTable.LMCP_flex(grasp_time) - pi/2))*K*K;
LPIP = R*palm_plane'.*lengths(18);
K = [0 -palm_plane(3) palm_plane(2); palm_plane(3) 0 -palm_plane(1); -palm_plane(2) palm_plane(1) 0];
R = eye(3) + sin(-unproTable.LMCP_abd(grasp_time))*K + (1-cos(-unproTable.LMCP_abd(grasp_time)))*K*K;
LPIP = R*LPIP + knuckles(5,:)';
K = [0 -new_plane(3) new_plane(2); new_plane(3) 0 -new_plane(1); -new_plane(2) new_plane(1) 0];
R = eye(3) + sin(unproTable.LPIP_flex(grasp_time)-pi)*K + (1-cos(unproTable.LPIP_flex(grasp_time)-pi))*K*K;
LDIP = R*(LPIP-knuckles(5,:)').*(lengths(19)/lengths(18)) + LPIP;
R = eye(3) + sin(unproTable.LDIP_flex(grasp_time)-pi)*K + (1-cos(unproTable.LDIP_flex(grasp_time)-pi))*K*K;
LT = R*(LDIP-LPIP).*(lengths(20)/lengths(19)) + LDIP;


%TCMC_flex is the angle between TCMC-IMCP and TCMC-TMCP along the thumb plane (RW-TCMC-IMCP)
%TMCP_abd is the angle between TCMC-IMCP and TCMC-TMCP out of the thumb plane (along the tabd plane described by the cross product between the thumb pane normal and TCMC-IMCP)
%TCMC_rot is the angle between the component of the thumb rotation plane (TCMC-TMCP-TIP) that's orthogonal to TCMC-IMCP, and the tabd plane


thumb_rot_axis = knuckles(2,:) - knuckles(1,:);
thumb_plane = cross(knuckles(1,:),knuckles(2,:));
thumb_plane = (thumb_plane./norm(thumb_plane));
K = [0 -thumb_plane(3) thumb_plane(2); thumb_plane(3) 0 -thumb_plane(1); -thumb_plane(2) thumb_plane(1) 0];
R = eye(3) + sin(-unproTable.TCMC_flex(grasp_time))*K + (1-cos(-unproTable.TCMC_flex(grasp_time)))*K*K;
TMCP = (R*thumb_rot_axis')'.*lengths(2)/norm(thumb_rot_axis);
thumb_abd_axis = cross(thumb_plane,thumb_rot_axis);
thumb_abd_axis = (thumb_abd_axis./norm(thumb_abd_axis));
K = [0 -thumb_abd_axis(3) thumb_abd_axis(2); thumb_abd_axis(3) 0 -thumb_abd_axis(1); -thumb_abd_axis(2) thumb_abd_axis(1) 0];
R = eye(3) + sin(-unproTable.TMCP_abd(grasp_time))*K + (1-cos(-unproTable.TMCP_abd(grasp_time)))*K*K;
TMCP = (R*TMCP')';

K = [0 -thumb_rot_axis(3) thumb_rot_axis(2); thumb_rot_axis(3) 0 -thumb_rot_axis(1); -thumb_rot_axis(2) thumb_rot_axis(1) 0];
R = eye(3) + sin(-unproTable.TCMC_rot(grasp_time))*K + (1-cos(-unproTable.TCMC_rot(grasp_time)))*K*K;
thumb_plane2 = (R*thumb_abd_axis')'./norm(thumb_abd_axis);
thumb_plane2_np = TMCP.*((thumb_plane2(:,1).*TMCP(:,1) + thumb_plane2(:,2).*TMCP(:,2) + thumb_plane2(:,3).*TMCP(:,3))./(norm(TMCP).^2));  % Vector component along TCMC-IMCP
thumb_plane2 = [thumb_plane2(:,1)-thumb_plane2_np(:,1),thumb_plane2(:,2)-thumb_plane2_np(:,2),thumb_plane2(:,3)-thumb_plane2_np(:,3)];  % Vector component orthogonal to TCMC-IMCP
thumb_plane2 = (thumb_plane2./norm(thumb_plane2));
TMCP = TMCP + knuckles(1,:)';

% vmag = dot(TMCP,thumb_rot_axis)/(norm(thumb_rot_axis)^2);
% v_inline = vmag*thumb_rot_axis;
% v_perp = TMCP - v_inline;
% thumb_rot_axis_unit = thumb_rot_axis/norm(thumb_rot_axis);
% K = [0 -thumb_rot_axis_unit(3) thumb_rot_axis_unit(2); thumb_rot_axis_unit(3) 0 -thumb_rot_axis_unit(1); -thumb_rot_axis_unit(2) thumb_rot_axis_unit(1) 0];
% R = eye(3) + sin(pi/2-unproTable.TCMC_rot(grasp_time))*K + (1-cos(pi/2-unproTable.TCMC_rot(grasp_time)))*K*K;
% TMCP = R*v_perp' + v_inline' + knuckles(1,:)';
% thumb_plane2 = R*thumb_plane';

K = [0 -thumb_plane2(3) thumb_plane2(2); thumb_plane2(3) 0 -thumb_plane2(1); -thumb_plane2(2) thumb_plane2(1) 0];
R = eye(3) + sin(unproTable.TMCP_flex(grasp_time)-pi)*K + (1-cos(unproTable.TMCP_flex(grasp_time)-pi))*K*K;
TIP = R*(TMCP-knuckles(1,:)').*(lengths(3)/lengths(2)) + TMCP;
R = eye(3) + sin(unproTable.TIP_flex(grasp_time)-pi)*K + (1-cos(unproTable.TIP_flex(grasp_time)-pi))*K*K;
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
arm_plane = R*palm_plane';
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
% RightArmPlot = plot3([0,E(1),S(1)], ...
%     [0,E(2),S(2)], ...
%     [0,E(3),S(3)], ...
%     '-o', 'MarkerSize',3,'MarkerFaceColor',"#D95319", 'Color',"#D95319");
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
% RightArmPlot = plot3([alignment_transformed(63),alignment_transformed(66),alignment_transformed(69),alignment_transformed(72)], ...
%      [alignment_transformed(61),alignment_transformed(64),alignment_transformed(67),alignment_transformed(70)], ...
%      [alignment_transformed(62),alignment_transformed(65),alignment_transformed(68),alignment_transformed(71)], ...
%      '-o', 'MarkerSize',3,'MarkerFaceColor',"#D95319", 'Color','k');
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