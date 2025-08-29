close all; clear all; clc;
% y = [488, 512, 358, 597, 517, 513, 504, 573]; %GP
y = [551, 750, 449, 734, 568, 801, 568, 780, 528, 749, 512, 751, 579, 803, 527, 790]; %PP
% x = ["P1 Grasp" "P1 Pinch" "P1 Grasp PCA" "P1 Pinch PCA" "P2 Grasp" "P2 Pinch" "P2 Grasp PCA" "P2 Pinch PCA"]; %GP
x = ["P1 Grasp Pre" "P1 Grasp Post" "P1 Grasp Pre PCA" "P1 Grasp Post PCA" "P1 Pinch Pre" "P1 Pinch Post" ...
    "P1 Pinch Pre PCA" "P1 Pinch Post PCA" "P2 Grasp Pre" "P2 Grasp Post" "P2 Grasp Pre PCA" "P2 Grasp Post PCA" ...
    "P2 Pinch Pre" "P2 Pinch Post" "P2 Pinch Pre PCA" "P2 Pinch Post PCA"]; %PP
b = bar(x,y);
b.FaceColor = 'flat';
b.CData(2:2:16,:) = repmat([.9 .3 .1],8,1);
% ylim([0 930])
% legend(["Grasp Correct", "Pinch Correct"])

hold on

err = 930 * (y./930) .* (1-(y./930));
errlow = err/-2;
errhigh = err/2;
er = errorbar(categorical(x),y,errlow,errhigh);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';

figstr = 'PP_LDA_Score';
titlestr = 'Pre/Post LDA Classification Score';
title(titlestr)
ylabel('Number of successful guesses out of 930')
xlabel('Category')
ylim([0 930])
fontsize(14,"points")
save_folder = 'C:\Users\czhe0008\Documents\EEG\LDA\';
saveas(gcf,strcat(save_folder,figstr,'.png'))