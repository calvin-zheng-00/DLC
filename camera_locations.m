close all; clear all; clc;

%x = cos(yaw)*cos(pitch)
%y = sin(yaw)*cos(pitch)
%z = sin(pitch)

%830E
%831F
%8337
%833C
%834C
%8380
%8389

%Err = 0.38
rotation = [ -0.07181233229369766, 0.015279590711491415, -0.0066928385808217895;
    0.3225769706715982, -2.0834137490616977, -0.9730116017647109;
    -0.0013849002559040363, -1.07860599191417, 0.054503854481861985;
    0.8410567495575895, -1.4559956658837143, -0.8365845847059495;
    0.408637822257164, -0.5464004987672559, -0.45273310516734977;
    -0.07680574381798971, 1.2517246162624618, -0.1683097202399387;
    0.4843403689317798, 0.4017221053302444, -0.12083787562634772];

translation = [ -3.8575871522589518, -15.37429290546911, 811.5948612587762;
    212.5285603702377, -140.64649189711326, 942.5865133043698;
    214.22876350401594, 24.697249356055252, 861.1415892116283;
    199.94068263217991, -116.04527934309331, 748.0548962088026;
    71.53752912995607, -104.85621839444416, 744.8280784800652;
    -148.6372340409785, -25.394184186920214, 820.6781364816705;
    -170.13837788799313, -14.161518304213828, 921.1877443187309];

figure

% cos(yaw)*cos(pitch)
% sin(yaw)*cos(pitch)
% sin(pitch)

% pitch
% roll
% yaw

r_x = rotation(:,1);
r_y = rotation(:,2);
r_z = rotation(:,3);

r_x2 = cos(rotation(:,3)).*cos(rotation(:,1));
r_y2 = sin(rotation(:,3)).*cos(rotation(:,1));
r_z2 = sin(rotation(:,1));

r_x3 = cos(rotation(:,1)).*cos(rotation(:,2));
r_y3 = sin(rotation(:,1)).*cos(rotation(:,2));
r_z3 = sin(rotation(:,2));

r_x4 = cos(rotation(:,2)).*cos(rotation(:,1));
r_y4 = sin(rotation(:,2)).*cos(rotation(:,1));
r_z4 = sin(rotation(:,1));

r_x5 = cos(rotation(:,1)).*cos(rotation(:,3));
r_y5 = sin(rotation(:,1)).*cos(rotation(:,3));
r_z5 = sin(rotation(:,3));

r_x6 = cos(rotation(:,2)).*cos(rotation(:,3));
r_y6 = sin(rotation(:,2)).*cos(rotation(:,3));
r_z6 = sin(rotation(:,3));

r_x7 = cos(rotation(:,3)).*cos(rotation(:,2));
r_y7 = sin(rotation(:,3)).*cos(rotation(:,2));
r_z7 = sin(rotation(:,2));

%quiver3(translation(:,1),translation(:,2),translation(:,3),rotation(:,1),rotation(:,2),rotation(:,3))
quiver3(translation(:,1),translation(:,2),translation(:,3),r_x,r_y,r_z)
hold on
plot3(translation(1,1),translation(1,2),translation(1,3),'go')%830E
plot3(translation(2,1),translation(2,2),translation(2,3),'mo')%831F
plot3(translation(3,1),translation(3,2),translation(3,3),'yo')%8337
plot3(translation(4,1),translation(4,2),translation(4,3),'co')%833C
plot3(translation(5,1),translation(5,2),translation(5,3),'bo')%834C
plot3(translation(6,1),translation(6,2),translation(6,3),'ko')%8380
plot3(translation(7,1),translation(7,2),translation(7,3),'ro')%8389
xlabel('x')
ylabel('y')
zlabel('z')

figure
quiver3(translation(:,1),translation(:,2),translation(:,3),r_x2,r_y2,r_z2)
hold on
plot3(translation(1,1),translation(1,2),translation(1,3),'go')%830E
plot3(translation(2,1),translation(2,2),translation(2,3),'mo')%831F
plot3(translation(3,1),translation(3,2),translation(3,3),'yo')%8337
plot3(translation(4,1),translation(4,2),translation(4,3),'co')%833C
plot3(translation(5,1),translation(5,2),translation(5,3),'bo')%834C
plot3(translation(6,1),translation(6,2),translation(6,3),'ko')%8380
plot3(translation(7,1),translation(7,2),translation(7,3),'ro')%8389
xlabel('x')
ylabel('y')
zlabel('z')

figure
quiver3(translation(:,1),translation(:,2),translation(:,3),r_x3,r_y3,r_z3)
hold on
plot3(translation(1,1),translation(1,2),translation(1,3),'go')%830E
plot3(translation(2,1),translation(2,2),translation(2,3),'mo')%831F
plot3(translation(3,1),translation(3,2),translation(3,3),'yo')%8337
plot3(translation(4,1),translation(4,2),translation(4,3),'co')%833C
plot3(translation(5,1),translation(5,2),translation(5,3),'bo')%834C
plot3(translation(6,1),translation(6,2),translation(6,3),'ko')%8380
plot3(translation(7,1),translation(7,2),translation(7,3),'ro')%8389
xlabel('x')
ylabel('y')
zlabel('z')

figure
quiver3(translation(:,1),translation(:,2),translation(:,3),r_x4,r_y4,r_z4)
hold on
plot3(translation(1,1),translation(1,2),translation(1,3),'go')%830E
plot3(translation(2,1),translation(2,2),translation(2,3),'mo')%831F
plot3(translation(3,1),translation(3,2),translation(3,3),'yo')%8337
plot3(translation(4,1),translation(4,2),translation(4,3),'co')%833C
plot3(translation(5,1),translation(5,2),translation(5,3),'bo')%834C
plot3(translation(6,1),translation(6,2),translation(6,3),'ko')%8380
plot3(translation(7,1),translation(7,2),translation(7,3),'ro')%8389
xlabel('x')
ylabel('y')
zlabel('z')

figure
quiver3(translation(:,1),translation(:,2),translation(:,3),r_x5,r_y5,r_z5)
hold on
plot3(translation(1,1),translation(1,2),translation(1,3),'go')%830E
plot3(translation(2,1),translation(2,2),translation(2,3),'mo')%831F
plot3(translation(3,1),translation(3,2),translation(3,3),'yo')%8337
plot3(translation(4,1),translation(4,2),translation(4,3),'co')%833C
plot3(translation(5,1),translation(5,2),translation(5,3),'bo')%834C
plot3(translation(6,1),translation(6,2),translation(6,3),'ko')%8380
plot3(translation(7,1),translation(7,2),translation(7,3),'ro')%8389
xlabel('x')
ylabel('y')
zlabel('z')

figure
quiver3(translation(:,1),translation(:,2),translation(:,3),r_x6,r_y6,r_z6)
hold on
plot3(translation(1,1),translation(1,2),translation(1,3),'go')%830E
plot3(translation(2,1),translation(2,2),translation(2,3),'mo')%831F
plot3(translation(3,1),translation(3,2),translation(3,3),'yo')%8337
plot3(translation(4,1),translation(4,2),translation(4,3),'co')%833C
plot3(translation(5,1),translation(5,2),translation(5,3),'bo')%834C
plot3(translation(6,1),translation(6,2),translation(6,3),'ko')%8380
plot3(translation(7,1),translation(7,2),translation(7,3),'ro')%8389
xlabel('x')
ylabel('y')
zlabel('z')

figure
quiver3(translation(:,1),translation(:,2),translation(:,3),r_x7,r_y7,r_z7)
hold on
plot3(translation(1,1),translation(1,2),translation(1,3),'go')%830E
plot3(translation(2,1),translation(2,2),translation(2,3),'mo')%831F
plot3(translation(3,1),translation(3,2),translation(3,3),'yo')%8337
plot3(translation(4,1),translation(4,2),translation(4,3),'co')%833C
plot3(translation(5,1),translation(5,2),translation(5,3),'bo')%834C
plot3(translation(6,1),translation(6,2),translation(6,3),'ko')%8380
plot3(translation(7,1),translation(7,2),translation(7,3),'ro')%8389
xlabel('x')
ylabel('y')
zlabel('z')

%Err=0.99
rotation2 = [ -0.01767508690520194, -0.005808548209518449, -0.005277101014853847;
    0.7270164201049648, -2.175524208436703, -0.6996642896099586;
    0.039409732877158536, -1.0368140742229455, 0.10514421444860234;
    0.996215717912918, -1.5296030612168305, -0.7622067528912126;
    0.44063311587646675, -0.5525613195300977, -0.39831514830192766;
    -0.03871406807882289, 1.4510428204315244, -0.14447125945780814;
    0.6659293035171362, 0.4734423752242747, -0.11582258408131459];
translation2 = [ 4.831297770649868, -14.211952915337841, 46.71096704318427;
    547.7746150382152, -309.6671806417661, 1430.0270648293138;
    866.5103317371418, 101.16215754095909, 532.3143097108809;
    820.9953232753862, -178.18152666700504, 929.4231161417605;
    519.2085607231639, 116.4363620450521, 208.44631588894129;
    -720.979391648478, 39.075566864831714, 694.1264921145132;
    -392.6697699282378, 351.61885353585797, 467.2416315131919];

figure
quiver3(translation2(:,1),translation2(:,2),translation2(:,3),rotation2(:,1),rotation2(:,2),rotation2(:,3))
hold on
plot3(translation2(1,1),translation2(1,2),translation2(1,3),'go')%830E
plot3(translation2(2,1),translation2(2,2),translation2(2,3),'mo')%831F
plot3(translation2(3,1),translation2(3,2),translation2(3,3),'yo')%8337
plot3(translation2(4,1),translation2(4,2),translation2(4,3),'co')%833C
plot3(translation2(5,1),translation2(5,2),translation2(5,3),'bo')%834C
plot3(translation2(6,1),translation2(6,2),translation2(6,3),'ko')%8380
plot3(translation2(7,1),translation2(7,2),translation2(7,3),'ro')%8389
xlabel('x')
ylabel('y')
zlabel('z')