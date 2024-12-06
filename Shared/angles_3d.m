close all; clear all; clc;
%src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\angles_participant_2\p10\Instructions_BAT1\20240730T101748-101755_angles.csv';
%src_3d = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\filter_participant\p10\Instructions_BAT1\20240730T101748-101755_filtered.csv';
src_3d = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\filter_participant\p21\Instructions_DIG1\20240816T134146-134155_filtered.csv';
src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\angles_participant_2\p21\Instructions_DIG1\20240816T134146-134155_angles.csv';
%src_3d = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\filter_participant\p23\Instructions_COO6\20240902T110850-110901_filtered.csv';
%src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\angles_participant_2\p23\Instructions_COO6\20240902T110850-110901_angles.csv';
%src_3d = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\filter_participant\p19\Instructions_COO1\20240809T124742-124752_filtered.csv';
%src = 'C:\Users\czhe0008\Documents\DLCprojects\openpose\data_conversion\massive_3d\angles_participant_2\p19\Instructions_COO1\20240809T124742-124752_angles.csv';

unproTable = fullfile(src);
unproTable = readtable(unproTable);
%% Getting limb lengths
T = fullfile(src_3d);
T = readtable(T);
original = table2array(T);
original_table = T(:,[4:63,67:78]);
original = original(:,[4:63,67:78]);            
grasp_time = 164; %164, 120
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

for i = 2:5
    j = (i-1)*4+2;
    v2_np = palm_plane.*((knuckles(i,1).*palm_plane(1) + knuckles(i,2).*palm_plane(2) + knuckles(i,3).*palm_plane(3)));   % Vector component normal to the palm plane
    v2_p = [knuckles(i,1)-v2_np(1),knuckles(i,2)-v2_np(2),knuckles(i,3)-v2_np(3)];                                       % Vector component on the palm plane
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
    [x,y,z] = sph2cart(unproTable(grasp_time,i+16).(1),unproTable(grasp_time,(i-1)*3+1).(1),lengths(j));
    PIP_temp = [x,y,z];
    PIP_temp = (R_knuckle*PIP_temp')';

    finger_axis_np = palm_plane.*(PIP_temp(1).*palm_plane(1) + PIP_temp(2).*palm_plane(2) + PIP_temp(3).*palm_plane(3));  % Vector component normal to the thumb plane
    finger_axis_p = [PIP_temp(1)-finger_axis_np(1),PIP_temp(2)-finger_axis_np(2),PIP_temp(3)-finger_axis_np(3)];                                    % Vector component on the thumb plane
    finger_axismag = sqrt(finger_axis_p(1).^2 + finger_axis_p(2).^2 + finger_axis_p(3).^2);
    finger_axis = finger_axis_p./finger_axismag;
    K = [0 -palm_plane(3) palm_plane(2); palm_plane(3) 0 -palm_plane(1); -palm_plane(2) palm_plane(1) 0];
    R = eye(3) + sin(-pi/2)*K + (1-cos(-pi/2))*K*K;
    finger_axis = (R*finger_axis')';
    PIP = [PIP;PIP_temp + knuckles(i,:)];
    K = [0 -finger_axis(3) finger_axis(2); finger_axis(3) 0 -finger_axis(1); -finger_axis(2) finger_axis(1) 0];
    R = eye(3) + sin(unproTable(grasp_time,(i-1)*3+2).(1))*K + (1-cos(unproTable(grasp_time,(i-1)*3+2).(1)))*K*K;
    DIP = [DIP;(R*(PIP(i-1,:)-knuckles(i,:))')'.*(lengths(j+1)/lengths(j)) + PIP(i-1,:)];
    R = eye(3) + sin(unproTable(grasp_time,(i-1)*3+3).(1))*K + (1-cos(unproTable(grasp_time,(i-1)*3+3).(1)))*K*K;
    Tip = [Tip;(R*(DIP(i-1,:)-PIP(i-1,:))')'.*(lengths(j+2)/lengths(j+1)) + DIP(i-1,:)];
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

[x,y,z] = sph2cart(unproTable.TCMC_flex(grasp_time),unproTable.TMCP_abd(grasp_time),lengths(2));
TMCP = [x,y,z];
TMCP = (R_thumb*TMCP')';
TMCP_mag = sqrt(TMCP(:,1).^2 + TMCP(:,2).^2 + TMCP(:,3).^2);

thumb_plane_np = TMCP.*((thumb_plane(:,1).*TMCP(:,1) + thumb_plane(:,2).*TMCP(:,2) + thumb_plane(:,3).*TMCP(:,3))./TMCP_mag.^2);  % thumb plane component along v1
thumb_plane_p = thumb_plane - thumb_plane_np;                               % thumb plane component orthogonal to v1
thumb_plane2 = thumb_plane_p./norm(thumb_plane_p);
TMCP_norm = TMCP./TMCP_mag;
K = [0 -TMCP_norm(3) TMCP_norm(2); TMCP_norm(3) 0 -TMCP_norm(1); -TMCP_norm(2) TMCP_norm(1) 0];
R = eye(3) + sin(unproTable.TCMC_rot(grasp_time))*K + (1-cos(unproTable.TCMC_rot(grasp_time)))*K*K;
thumb_plane2 = (R*thumb_plane2')';
thumb_plane2 = -thumb_plane2./norm(thumb_plane2);
TMCP = TMCP + knuckles(1,:);

K = [0 -thumb_plane2(3) thumb_plane2(2); thumb_plane2(3) 0 -thumb_plane2(1); -thumb_plane2(2) thumb_plane2(1) 0];
R = eye(3) + sin(unproTable.TMCP_flex(grasp_time))*K + (1-cos(unproTable.TMCP_flex(grasp_time)))*K*K;
TIP = (R*(TMCP-knuckles(1,:))')'.*(lengths(3)/lengths(2)) + TMCP;
R = eye(3) + sin(unproTable.TIP_flex(grasp_time))*K + (1-cos(unproTable.TIP_flex(grasp_time)))*K*K;
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
[x,y,z] = sph2cart(unproTable.W_abd(grasp_time),unproTable.W_flex(grasp_time),lengths(23));
E = [x,y,z];
E = (R_Elbow*E')';
%E = E*lengths(23)/norm(E);

forearm_np = E.*((palm_plane(:,1).*E(:,1) + palm_plane(:,2).*E(:,2) + palm_plane(:,3).*E(:,3))./(lengths(23).^2));   % Vector component normal to the palm plane
forearm_p = [palm_plane(:,1)-forearm_np(:,1),palm_plane(:,2)-forearm_np(:,2),palm_plane(:,3)-forearm_np(:,3)];                                       % Vector component on the palm plane
forearmmag = sqrt(forearm_p(:,1).^2 + forearm_p(:,2).^2 + forearm_p(:,3).^2);
forearm_ref = forearm_p./forearmmag;
E_rot = E/lengths(23);
K = [0 -E_rot(3) E_rot(2); E_rot(3) 0 -E_rot(1); -E_rot(2) E_rot(1) 0];
R = eye(3) + sin(unproTable.W_rot(grasp_time))*K + (1-cos(unproTable.W_rot(grasp_time)))*K*K;
E_axis = -(R*forearm_ref')';
K = [0 -E_axis(3) E_axis(2); E_axis(3) 0 -E_axis(1); -E_axis(2) E_axis(1) 0];
% R = eye(3) + sin(unproTable.RE_flex(grasp_time)-pi/2)*K + (1-cos(unproTable.RE_flex(grasp_time)-pi/2))*K*K;
R = eye(3) + sin(unproTable.RE_flex(grasp_time))*K + (1-cos(unproTable.RE_flex(grasp_time)))*K*K;
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

[x,y,z] = sph2cart(unproTable.RS_flex(grasp_time),unproTable.RS_abd(grasp_time),lengths(21));
C = [x,y,z];
C = (R_Shoulder*C')' + S;

chest_axis = C-S;
chest_axis = chest_axis/norm(chest_axis);
shoulder_norm = cross(chest_axis,S-E,2);
shoulder_normmag = sqrt(shoulder_norm(:,1).^2 + shoulder_norm(:,2).^2 + shoulder_norm(:,3).^2);
unit_shoulder_norm = shoulder_norm./shoulder_normmag;
K = [0 -chest_axis(3) chest_axis(2); chest_axis(3) 0 -chest_axis(1); -chest_axis(2) chest_axis(1) 0];
R = eye(3) + sin(unproTable.RS_rot(grasp_time))*K + (1-cos(unproTable.RS_rot(grasp_time)))*K*K;
N = (R*unit_shoulder_norm')'*-100+C;

%% Drawing reprojection
figure
hold on;
% New hand
ThumbPlot = plot3([0,knuckles(1,1),TMCP(1),TIP(1),TT(1)], ...
    [0,knuckles(1,2),TMCP(2),TIP(2),TT(2)], ...
    [0,knuckles(1,3),TMCP(3),TIP(3),TT(3)], ...
    '-o', 'MarkerSize',3,'MarkerFaceColor',	'r', 'Color','r');
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

% xlim([-200 100])
% ylim([-100 200])
% zlim([-200 100])
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
flag = 1;


%% Drawing original
alignment_transformed = [];
for l = 1:24
    temp = transpose(original_R*alignment(:,(l*3-2):l*3)');%.*[1,-1,1];
    alignment_transformed = cat(2,alignment_transformed,temp);
end
Nose = [T.N_x(grasp_time),T.N_y(grasp_time),T.N_z(grasp_time)];
Nose = Nose - original(grasp_time,70:72);
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
% end